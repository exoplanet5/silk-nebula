"""verify_shockmodel_validate: validate the shockmodel network/cooling/shock code against published models.

 1. CIE ion fractions vs Gnat & Sternberg 2007 (GS07, VizieR J/ApJS/168/213 ion_CIE; Asplund05 abundances).
 2. Isobaric non-equilibrium cooling from 5e6 K vs GS07 ion_bar (Z = 1) and cef_bar.
 3. Complete steady shocks vs Allen et al. 2008 Table 6 (MAPPINGS III, solar, n = 1, B = 3.23 uG), 100 & 200 km/s.
Output: data/derived/verify_shockmodel_validation.json, figures/verify_shockmodel_validation.png
"""
import sys, json, numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import verify_shockmodel_atomic as A
import verify_shockmodel_flow as F

ROOT = '/Users/mickey/astronomy/silk-nebula'
LIT = ROOT + '/data/raw/verify_shockmodel_lit'


def read_gs07(fn):
    lines = [l.rstrip('\n') for l in open(f'{LIT}/{fn}') if not l.startswith('#') and l.strip()]
    head = lines[0].split('\t')
    rows = []
    for l in lines[3:]:
        s = l.split('\t')
        try:
            rows.append([float(x) if x.strip() not in ('', 'Sun') else np.nan for x in s])
        except ValueError:
            continue
    return head, np.array(rows)


def cie_compare(net):
    head, a = read_gs07('gs07_ion_CIE.tsv')
    col = {h: i for i, h in enumerate(head)}
    out = {}
    for el, Z in (('He', 2), ('C', 6), ('N', 7), ('O', 8)):
        rows = []
        for T in (2e4, 3e4, 5e4, 1e5, 2e5, 3e5, 5e5, 1e6):
            k = np.argmin(np.abs(a[:, col['T']] - T))
            f = net.cie_state(a[k, col['T']])
            mine = f[net.slices[el]]
            names = [f'{el}0/{el}'] + [f'{el}+/{el}'] + [f'{el}{q}+/{el}' for q in range(2, Z + 1)]
            gs = np.array([a[k, col[n]] for n in names])
            rows.append({'T': float(a[k, col['T']]), 'mine': mine.round(4).tolist(), 'gs07': gs.tolist()})
        out[el] = rows
    return out


def isobaric(net, T0=5e6, nH0=1e-3, Tstop=1.2e4):
    f = net.cie_state(T0)
    P = nH0 * net.ntot_per_H(f) * F.K_B * T0
    h = 2.5 * P / (nH0 * net.mass_per_H)
    T, nH = T0, nH0
    rec = []
    dt = 1e8
    while T > Tstop:
        L = net.cooling(f, T, nH)
        for _ in range(40):
            rho = nH * net.mass_per_H
            hn = h - L / rho * dt
            rhon = 2.5 * P / hn
            nHn = rhon / net.mass_per_H
            fn = net.implicit_ion_step(f, T, 0.5 * (nH + nHn), dt)
            Tn = P / (nHn * net.ntot_per_H(fn) * F.K_B)
            if abs(Tn - T) / T > 0.003 or np.max(np.abs(fn - f)) > 0.02:
                dt *= 0.5; continue
            break
        rec.append((T, f.copy(), L / (nH * nH * net.ne_per_H(f))))
        f, T, nH, h = fn, Tn, nHn, hn
        dt *= 1.3
    return rec


def bar_compare(net):
    head, a = read_gs07('gs07_ion_bar.tsv')
    col = {h: i for i, h in enumerate(head)}
    # The VizieR Z column is blank/inconsistent for most rows; blocks are delimited by T resets to 5e6 K.
    # Block 0 (rows 1-130) is identified as Z = 1 by its own GF12-weighted cooling efficiency matching the
    # GS07 cef_bar Z=1 column (7.4e-23 vs 7.8e-23 at 3e4 K; 1.9e-22 vs 2.2e-22 at 1e5 K); blocks 1-3 are
    # near-CIE (low Z) and block 4 is the most over-ionized (Z = 2).  See notes.
    a1 = a[0:130]
    hc, c = read_gs07('gs07_cef_bar.tsv')
    colc = {h: i for i, h in enumerate(hc)}
    rec = isobaric(net)
    Tm = np.array([r[0] for r in rec])
    Lm = np.array([r[2] for r in rec])
    out = {'T': [], 'O_mine': [], 'O_gs07': [], 'C_mine': [], 'C_gs07': [], 'Lam_mine': [], 'Lam_gs07': []}
    for T in (2e4, 3e4, 5e4, 7e4, 1e5, 1.5e5, 2e5, 3e5, 5e5, 1e6):
        k = np.argmin(np.abs(a1[:, col['T']] - T))
        Tg = a1[k, col['T']]
        km = np.argmin(np.abs(Tm - Tg))
        fo = rec[km][1][net.slices['O']]
        fc = rec[km][1][net.slices['C']]
        gO = [a1[k, col[n]] for n in ['O0/O', 'O+/O'] + [f'O{q}+/O' for q in range(2, 9)]]
        gC = [a1[k, col[n]] for n in ['C0/C', 'C+/C'] + [f'C{q}+/C' for q in range(2, 7)]]
        kc = np.argmin(np.abs(c[:, colc['T']] - Tg))
        out['T'].append(float(Tg))
        out['O_mine'].append(fo.round(4).tolist()); out['O_gs07'].append(gO)
        out['C_mine'].append(fc.round(4).tolist()); out['C_gs07'].append(gC)
        out['Lam_mine'].append(float(Lm[km])); out['Lam_gs07'].append(float(c[kc, colc['Lambda(Z=1)']]))
    out['curve_T'] = Tm[::5].tolist(); out['curve_L'] = Lm[::5].tolist()
    out['curve_O2'] = [float(r[1][net.slices['O']][2]) for r in rec[::5]]
    out['gs07_T'] = a1[:, col['T']].tolist(); out['gs07_O2'] = a1[:, col['O2+/O']].tolist()
    out['gs07_cefT'] = c[:, colc['T']].tolist(); out['gs07_cef'] = c[:, colc['Lambda(Z=1)']].tolist()
    return out


ALLEN6 = {  # Allen et al. 2008 Table 6 (shock component, solar, n=1, B=3.23 uG), relative to Hbeta
    100: {'CIV1549': 15.032, 'HeII1640': 0.469, 'OIII1664': 1.735, 'OII3727': 10.652, 'OIII4363': 0.290,
          'HeII4686': 0.039, 'OIII5007': 3.190, 'Ha': 3.214, 'NII6583': 1.641, 'SII6716_31': 1.445,
          'NeIII3869': 0.612, 'logHb': -5.317},
    200: {'CIV1549': 4.593, 'HeII1640': 0.479, 'OIII1664': 0.711, 'OII3727': 5.964, 'OIII4363': 0.115,
          'HeII4686': 0.060, 'OIII5007': 1.188, 'Ha': 3.010, 'NII6583': 2.024, 'SII6716_31': 1.975,
          'NeIII3869': 0.309, 'logHb': -4.444},
}


def ratios_from(F_):
    hb = F_['Hb'][-1]
    return {'CIV1549': (F_['CIV1548'][-1] + F_['CIV1551'][-1]) / hb, 'HeII1640': F_['HeII1640'][-1] / hb,
            'OIII1664': (F_['OIII1666'][-1] + F_['OIII1661'][-1]) / hb,
            'OII3727': (F_['OII3726'][-1] + F_['OII3729'][-1]) / hb, 'OIII4363': F_['OIII4363'][-1] / hb,
            'HeII4686': F_['HeII4686'][-1] / hb, 'OIII5007': F_['OIII5007'][-1] / hb, 'Ha': F_['Ha'][-1] / hb,
            'NII6583': F_['NII6583'][-1] / hb, 'SII6716_31': (F_['SII6716'][-1] + F_['SII6731'][-1]) / hb,
            'NeIII3869': F_['NeIII3869'][-1] / hb, 'logHb': float(np.log10(hb))}


def allen_compare():
    net = F.Network('allen08')
    res = {}
    for v in (100, 200):
        for pre, Tf in (('precursor', 1e4), ('wim', 1e4), ('precursor', 7e3)):
            sh = F.Shock(net, v, 1.0, B0_uG=3.23, pre=pre, complete=True, Tfloor=Tf).run()
            j = F.emissivities(sh, thin=2)
            Fc = F.cumulative(sh, j)
            r = ratios_from(Fc)
            res[f'{v}_{pre}_{int(Tf)}'] = {'mine': {k: float(x) for k, x in r.items()}, 'allen08': ALLEN6[v],
                                 'Ts': float(sh.Ts), 'Nfinal': float(sh.N[-1]), 'nstep': int(sh.nstep)}
            print(v, pre, Tf, {k: (round(r[k], 2), ALLEN6[v][k]) for k in r})
    return res


if __name__ == '__main__':
    net05 = F.Network('asplund05')
    out = {'cie': cie_compare(net05)}
    for el in ('O', 'C'):
        for row in out['cie'][el]:
            print('CIE', el, f"{row['T']:.2e}", 'mine', np.round(row['mine'], 3), 'gs07', np.round(row['gs07'], 3))
    out['isobaric'] = bar_compare(net05)
    ib = out['isobaric']
    for i, T in enumerate(ib['T']):
        print('ISOBAR', f'{T:.2e}', 'O mine', np.round(ib['O_mine'][i], 3), 'gs07', np.round(ib['O_gs07'][i], 3),
              'Lam', f"{ib['Lam_mine'][i]:.2e} {ib['Lam_gs07'][i]:.2e}")
    out['allen'] = allen_compare()
    json.dump(out, open(ROOT + '/data/derived/verify_shockmodel_validation.json', 'w'), indent=1)
