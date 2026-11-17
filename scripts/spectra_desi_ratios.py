"""spectra_desi_ratios.py -- differential line ratios, velocities and widths of the [O III] filament component.

Builds on spectra_desi_regspec.py (regression spectrum S = sum w P c / sum w P^2, c = continuum-subtracted
residual in in-fibre R/A, P = instrument = expected [O III]5007 residual of the fibre in R):
  * fixed-velocity Gaussian fits of every line in S at the DESI [O III] velocity (and the LAMOST red-line
    velocity for the low-ionisation lines), sigma_int = 20 km/s, DESI LSF from spectra_desi_lsf.json;
  * healpix-block bootstrap (DESI nside-64 coadd pixels, ~0.9 deg; 400 resamples) for all amplitudes and ratios,
    which absorbs spatially correlated (per-petal) sky-subtraction systematics that ivar errors miss;
  * ratios in photon units (R/R, as measured) and energy units (x lambda_ref/lambda_line);
  * [O III] velocity and width from S; velocity in Dec bands; sign split (positive / negative residual fibres);
  * electron temperature limit from [O III] (4959+5007)/4363.
Outputs: data/derived/spectra_desi_ratios.json, data/derived/spectra_desi_ratios_boot.npz,
         data/derived/spectra_desi_regspec_full.npz (S, eS, bootstrap percentile band, stacks)
Usage  : python -I scripts/spectra_desi_ratios.py [--nboot 400]
"""
import csv
import json
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_desi_lines import LINES, WIN_WAVE, C_KMS  # noqa: E402
import spectra_desi_measure as meas  # noqa: E402
from spectra_desi_regspec import continuum_subtract, fit_lines_1d, fit_width, GROUPS_1D  # noqa: E402
from spectra_desi_analyze import rayleigh_per_unit  # noqa: E402

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
D = ROOT / 'data/derived'
RNG = np.random.default_rng(42)
V_LAMOST_RED = 15.0      # km/s helio, mean of LAMOST Halpha/[NII]/[SII] (Wu+2026 Table 4: 11.7-19.7)
V_LAMOST_OIII = 25.0     # km/s helio, LAMOST [O III] (23.8-26.8)


def design_fixed(cols, v, S, eS, sig_int=20.0, half=30.0):
    """Design matrix and selection for a fixed-velocity fit (Gaussians + constant)."""
    lam_c = np.array([c[0][0] for c in cols])
    sel = (WIN_WAVE >= lam_c.min() - half) & (WIN_WAVE <= lam_c.max() + half) & np.isfinite(S) & (eS > 0) & \
        np.isfinite(eS)
    lam = WIN_WAVE[sel]
    G = np.zeros((len(cols) + 1, lam.size))
    for k, c in enumerate(cols):
        for l0, amp in c:
            s = np.hypot(meas.lsf_sigma_A(l0), l0 * sig_int / C_KMS)
            mu = l0 * (1 + v / C_KMS)
            G[k] += amp * np.exp(-0.5 * ((lam - mu) / s) ** 2) / (np.sqrt(2 * np.pi) * s)
    G[-1] = 1.0
    w = 1 / eS[sel] ** 2
    A = (G * w) @ G.T
    cov = np.linalg.pinv(A)
    proj = cov @ (G * w)                    # amplitude = proj @ S[sel]
    return sel, proj, cov, G


def fit_fixed(S, eS, cols, v, **kw):
    sel, proj, cov, G = design_fixed(cols, v, S, eS, **kw)
    p = proj @ S[sel]
    chi2 = float(np.sum(((S[sel] - p @ G) / eS[sel]) ** 2))
    dof = int(sel.sum() - len(cols) - 1)
    return dict(amp=p[:-1].tolist(), err=np.sqrt(np.diag(cov))[:-1].tolist(), red_chi2=chi2 / max(dof, 1),
                v=float(v)), (sel, proj)


def block_sums(C, W, P, okP, blocks):
    """Per-block weighted sums for the per-pixel regression c(lambda) = a(lambda) + S(lambda) P."""
    ub, inv = np.unique(blocks, return_inverse=True)
    B = {k: np.zeros((ub.size, C.shape[1])) for k in ('s0', 's1', 's2', 't0', 't1')}
    p = np.where(okP, P, 0.0)
    for i0 in range(0, C.shape[0], 4000):
        sl = slice(i0, i0 + 4000)
        w = W[sl] * okP[sl, None]
        wp = w * p[sl, None]
        np.add.at(B['s0'], inv[sl], w)
        np.add.at(B['s1'], inv[sl], wp)
        np.add.at(B['s2'], inv[sl], wp * p[sl, None])
        np.add.at(B['t0'], inv[sl], w * C[sl])
        np.add.at(B['t1'], inv[sl], wp * C[sl])
    return ub, B


INTERCEPT = True    # fit c = a + S P per pixel (absorbs any uniform additive residual offset, e.g. at Halpha)


def reg(B, cnt=None, intercept=None):
    """Regression spectrum S(lambda) and its formal error from block sums (optionally block-resampled)."""
    intercept = INTERCEPT if intercept is None else intercept
    if cnt is None:
        sm = {k: v.sum(0) for k, v in B.items()}
    else:
        sm = {k: cnt @ v for k, v in B.items()}
    with np.errstate(invalid='ignore', divide='ignore'):
        if intercept:
            det = sm['s0'] * sm['s2'] - sm['s1'] ** 2
            S = (sm['s0'] * sm['t1'] - sm['s1'] * sm['t0']) / det
            eS = np.sqrt(sm['s0'] / det)
        else:
            S = sm['t1'] / sm['s2']
            eS = 1 / np.sqrt(sm['s2'])
    return S, eS


def te_from_oiii(R):
    """T_e (K) from R = (I4959+I5007)/I4363 (energy), low-density limit (Osterbrock & Ferland 2006 eq. 5.4)."""
    if not np.isfinite(R) or R <= 7.90:
        return np.nan
    # 7.90 exp(32900/T) / (1 + 4.5e-4 ne/T^0.5) with ne -> 0
    return 32900.0 / np.log(R / 7.90)


RATIO_NAMES = None


def ratio_table(fits, boots, L, vlab='vO', ref='OIII5007'):
    def b_(g, j=0):
        return boots[f'{g}@{vlab}'][:, j]

    def f_(g, j=0):
        return fits[f'{g}@{vlab}']['amp'][j]

    defs = {
        'OIII5007/Hb': (b_('OIII5007'), b_('Hb'), f_('OIII5007'), f_('Hb'), L['Hb'] / L['OIII5007']),
        'OIII5007/Ha': (b_('OIII5007'), b_('HaNII', 0), f_('OIII5007'), f_('HaNII', 0), L['Ha'] / L['OIII5007']),
        'Ha/Hb': (b_('HaNII', 0), b_('Hb'), f_('HaNII', 0), f_('Hb'), L['Hb'] / L['Ha']),
        'NII6583/Ha': (b_('HaNII', 1), b_('HaNII', 0), f_('HaNII', 1), f_('HaNII', 0), L['Ha'] / L['NII6583']),
        'SII/Ha': (b_('SII', 0) + b_('SII', 1), b_('HaNII', 0), f_('SII', 0) + f_('SII', 1), f_('HaNII', 0),
                   L['Ha'] / L['SII6716']),
        'SII6716/6731': (b_('SII', 0), b_('SII', 1), f_('SII', 0), f_('SII', 1), L['SII6731'] / L['SII6716']),
        'OI6300/Ha': (b_('OI6300'), b_('HaNII', 0), f_('OI6300'), f_('HaNII', 0), L['Ha'] / L['OI6300']),
        'HeII4686/Hb': (b_('HeII4686'), b_('Hb'), f_('HeII4686'), f_('Hb'), L['Hb'] / L['HeII4686']),
        'NI5200/Hb': (b_('NI5200'), b_('Hb'), f_('NI5200'), f_('Hb'), L['Hb'] / L['NI5198']),
        'HeI5876/Hb': (b_('HeI5876'), b_('Hb'), f_('HeI5876'), f_('Hb'), L['Hb'] / L['HeI5876']),
        'OII3729/3726': (b_('OII3727', 1), b_('OII3727', 0), f_('OII3727', 1), f_('OII3727', 0),
                         L['OII3726'] / L['OII3729']),
        'OIII5007/OIII4959': (b_('OIII5007_single'), b_('OIII4959'), f_('OIII5007_single'), f_('OIII4959'),
                              L['OIII4959'] / L['OIII5007']),
    }
    # everything relative to the reference line (OIII5007 or Ha)
    refb = b_('OIII5007') if ref == 'OIII5007' else b_('HaNII', 0)
    reff = f_('OIII5007') if ref == 'OIII5007' else f_('HaNII', 0)
    lref = L['OIII5007'] if ref == 'OIII5007' else L['Ha']
    for nm, g, j, lam in (('Hb', 'Hb', 0, L['Hb']), ('Ha', 'HaNII', 0, L['Ha']), ('NII6583', 'HaNII', 1, L['NII6583']),
                          ('SII6716+31', 'SII', 0, L['SII6716']), ('OI6300', 'OI6300', 0, L['OI6300']),
                          ('HeII4686', 'HeII4686', 0, L['HeII4686']), ('NI5200', 'NI5200', 0, L['NI5198']),
                          ('OII3727', 'OII3727', 0, L['OII3726']), ('OIII4363', 'OIII4363', 0, L['OIII4363']),
                          ('NeIII3869', 'NeIII3869', 0, L['NeIII3869']), ('SIII9069', 'SIII9069', 0, L['SIII9069']),
                          ('OII7325', 'OII7325', 0, L['OII7320']), ('HeI5876', 'HeI5876', 0, L['HeI5876']),
                          ('NII5755', 'NII5755', 0, L['NII5755']), ('ArIII7136', 'ArIII7136', 0, L['ArIII7136']),
                          ('OIII5007', 'OIII5007', 0, L['OIII5007']), ('sky5577', 'sky5577', 0, L['sky5577'])):
        if nm == ('OIII5007' if ref == 'OIII5007' else 'Ha'):
            continue
        nb_, nf_ = b_(g, j), f_(g, j)
        if g in ('SII', 'OII3727'):
            nb_, nf_ = b_(g, 0) + b_(g, 1), f_(g, 0) + f_(g, 1)
        defs[f'{nm}/{ref}'] = (nb_, refb, nf_, reff, lref / lam if nm != 'sky5577' else 1.0)
    out = {}
    for name, (nb, db, nf, df, efac) in defs.items():
        q = nb / db
        out[name] = dict(photon=nf / df, energy=nf / df * efac,
                         p16=float(np.nanpercentile(q, 16) * efac), p50=float(np.nanpercentile(q, 50) * efac),
                         p84=float(np.nanpercentile(q, 84) * efac), p2p5=float(np.nanpercentile(q, 2.5) * efac),
                         p97p5=float(np.nanpercentile(q, 97.5) * efac),
                         num_snr_boot=float(nf / np.nanstd(nb)), den_snr_boot=float(df / np.nanstd(db)))
    return out


def regress(C, W, P, okP, hp, nboot, ref='OIII5007', groups=None, vgrid_ref=np.arange(-60, 80.1, 0.5),
            extra_v=True, intercept=None):
    """Regression spectrum on instrument P (units: R of ref line) with healpix-block bootstrap and line fits."""
    ub, Bk = block_sums(C, W, P, okP, hp)
    S, eS = reg(Bk, intercept=intercept)
    Sb = np.empty((nboot, S.size))
    for b in range(nboot):
        cnt = np.bincount(RNG.integers(0, ub.size, ub.size), minlength=ub.size).astype(float)
        Sb[b] = reg(Bk, cnt, intercept=intercept)[0]
    gref = 'OIII5007' if ref == 'OIII5007' else 'HaNII'
    fR = fit_lines_1d(S, eS, GROUPS_1D[gref], vgrid=vgrid_ref)
    vR = fR['v']
    vb = np.array([fit_lines_1d(Sb[b], eS, GROUPS_1D[gref], vgrid=np.arange(-40, 70.1, 1.0))['v']
                   for b in range(min(nboot, 150))])
    r = {'n_fibres': int(okP.sum()), 'n_blocks': int(ub.size), 'ref': ref, 'v_ref_free': fR,
         'v_ref_boot': dict(median=float(np.median(vb)), p16=float(np.percentile(vb, 16)),
                            p84=float(np.percentile(vb, 84)))}
    if ref == 'OIII5007':
        wid = fit_width(S, eS, [LINES['OIII5007'], LINES['OIII4959']], [1, 1 / 2.98], vR)
    else:
        wid = fit_width(S, eS, [LINES['Ha']], [1], vR)
    l0 = LINES['OIII5007'] if ref == 'OIII5007' else LINES['Ha']
    lsf = meas.lsf_sigma_A(l0) / l0 * C_KMS
    wid['lsf_sigma_kms'] = lsf
    for q, qn in (('sigma_tot_kms', 'sigma_int_kms'), ('lo', 'sigma_int_lo'), ('hi', 'sigma_int_hi')):
        s2 = wid[q] ** 2 - lsf ** 2
        wid[qn] = float(np.sqrt(s2)) if s2 > 0 else 0.0
    wid['fwhm_int_kms'] = 2.3548 * wid['sigma_int_kms']
    for f_, lab in ((0.98, 'lsf_minus2pct'), (1.02, 'lsf_plus2pct')):
        s2 = wid['sigma_tot_kms'] ** 2 - (f_ * lsf) ** 2
        wid['fwhm_int_' + lab] = float(2.3548 * np.sqrt(s2)) if s2 > 0 else 0.0
    r['width'] = wid
    fits, boots = {}, {}
    vlist = [('vO', vR), ('vLAMOST_red', V_LAMOST_RED), ('vLAMOST_OIII', V_LAMOST_OIII)] if extra_v else [('vO', vR)]
    for vlab, vv in vlist:
        for g, cols in GROUPS_1D.items():
            f, (s_, proj) = fit_fixed(S, eS, cols, vv)
            ab = (Sb[:, s_] @ proj.T)[:, :-1]
            f['boot_err'] = np.nanstd(ab, 0).tolist()
            fits[f'{g}@{vlab}'] = f
            boots[f'{g}@{vlab}'] = ab
    r['fixed'] = fits
    r['ratios'] = ratio_table(fits, boots, LINES, 'vO', ref)
    if ref == 'OIII5007':
        L = LINES
        q4363 = boots['OIII4363@vO'][:, 0] / (boots['OIII5007@vO'][:, 0] * (1 + 1 / 2.98)) * (L['OIII5007'] / L['OIII4363'])
        R_full = (fits['OIII5007@vO']['amp'][0] * (1 + 1 / 2.98)) / fits['OIII4363@vO']['amp'][0] * \
            (L['OIII4363'] / L['OIII5007'])
        r['Te_OIII'] = dict(R_energy=float(R_full), Te_K=te_from_oiii(R_full),
                            Te_upper_97p5=te_from_oiii(1 / np.nanpercentile(q4363, 97.5)),
                            q4363_p2p5=float(np.nanpercentile(q4363, 2.5)), q4363_p50=float(np.nanpercentile(q4363, 50)),
                            q4363_p97p5=float(np.nanpercentile(q4363, 97.5)),
                            note='low-density limit; q4363 = I4363/(I4959+I5007), energy; bootstrap percentiles')
    band = dict(S=S, eS=eS, p16=np.nanpercentile(Sb, 16, 0), p84=np.nanpercentile(Sb, 84, 0))
    return r, band, boots


def main(nboot=400):
    from spectra_desi_analyze import neighbour_map
    ids, F, I, M, Sg, files = meas.load_windows()
    arr = np.load(D / 'spectra_desi_model_arrays.npz', allow_pickle=False)
    fitz = np.load(D / 'spectra_desi_fits.npz', allow_pickle=False)
    pos = {s: k for k, s in enumerate(arr['sparcl_id'])}
    assert (fitz['sparcl_id'] == arr['sparcl_id']).all()
    sel = np.array([i in pos for i in ids])
    ids, F, I, M = ids[sel], F[sel], I[sel], M[sel]
    k = np.array([pos[i] for i in ids])
    cat = {r['sparcl_id']: r for r in csv.DictReader(open(D / 'spectra_desi_catalog.csv'))}
    z = np.array([float(cat[i]['redshift']) for i in ids])
    ra = np.array([float(cat[i]['ra']) for i in ids])
    dec = np.array([float(cat[i]['dec']) for i in ids])
    hp = np.array([int(cat[i]['healpix']) for i in ids])
    dany = np.array([float(cat[i]['d_any']) for i in ids])
    ffrac = np.clip(np.array([float(cat[i]['mean_psf_to_fiber_specflux']) for i in ids]), 0.3, 1.0)
    C, W = continuum_subtract(F - M, I, z)
    del F, M
    kR = rayleigh_per_unit(WIN_WAVE)
    conv = ffrac[:, None] * kR[None, :]
    C = C * conv
    W = W / conv ** 2
    # H-alpha neighbour instrument (in R), from per-fibre fixed-velocity fits
    Ha = fitz['HaNII__amp_fix'][k, 0] * ffrac * rayleigh_per_unit(LINES['Ha'])
    eHa = fitz['HaNII__err_fix'][k, 0] * ffrac * rayleigh_per_unit(LINES['Ha'])
    okHa = np.isfinite(Ha) & np.isfinite(eHa) & (eHa > 0) & ~fitz['HaNII__flag_own'][k] & \
        (fitz['HaNII__npix_good'][k] > 40)
    nbHa, _ = neighbour_map(ra, dec, Ha, eHa, okHa, radius=0.12)
    trace = json.loads((ROOT / 'data/derived/silk_trace.json').read_text())
    patch = np.array(trace['halpha_patch_info']['outer_polygon_30pct'])
    from matplotlib.path import Path as MPath
    in_patch = MPath(patch).contains_points(np.c_[ra, dec])
    near_patch = np.hypot((ra - 115.6) * np.cos(np.radians(20.5)), dec - 20.5) < 1.2
    P_model, P_nbr = arr['pred'][k], arr['nbO'][k]
    out = {'n_fibres_total': int(len(ids)), 'nboot': nboot, 'n_in_patch30': int(in_patch.sum())}
    full = {'wave': WIN_WAVE}
    runs = [
        ('model_all', P_model, np.isfinite(P_model), 'OIII5007'),
        ('nbr_all', P_nbr, np.isfinite(P_nbr), 'OIII5007'),
        ('model_south', P_model, np.isfinite(P_model) & (dec < 19.0), 'OIII5007'),
        ('model_mid', P_model, np.isfinite(P_model) & (dec >= 19.0) & (dec < 22.0) & ~near_patch, 'OIII5007'),
        ('model_north', P_model, np.isfinite(P_model) & (dec >= 22.0), 'OIII5007'),
        ('model_nopatch', P_model, np.isfinite(P_model) & ~near_patch, 'OIII5007'),
        ('model_posP', P_model, np.isfinite(P_model) & (P_model > 0.3), 'OIII5007'),
        ('model_negP', P_model, np.isfinite(P_model) & (P_model < -0.3), 'OIII5007'),
        ('nbrHa_all', nbHa, np.isfinite(nbHa), 'Ha'),
        ('nbrHa_nopatch', nbHa, np.isfinite(nbHa) & ~near_patch, 'Ha'),
        ('nbrHa_patch', nbHa, np.isfinite(nbHa) & near_patch, 'Ha'),
    ]
    boots_all = {}
    for name, P, okP, ref in runs:
        if okP.sum() < 40:
            continue
        r, band, boots = regress(C[okP], W[okP], P[okP], np.ones(okP.sum(), bool), hp[okP],
                                 nboot if name in ('model_all', 'nbr_all', 'nbrHa_all') else max(nboot // 2, 100),
                                 ref=ref, extra_v=name in ('model_all', 'model_south'),
                                 intercept=False if name in ('model_posP', 'model_negP') else None)
        r['intercept'] = name not in ('model_posP', 'model_negP')
        out[name] = r
        for kk, vv in band.items():
            full[f'{kk}_{name}'] = vv
        if name == 'model_all':
            boots_all = boots
        print(f"== {name}: n={r['n_fibres']} blocks={r['n_blocks']} v_ref={r['v_ref_free']['v']:.1f} "
              f"[{r['v_ref_free']['v_lo']:.1f},{r['v_ref_free']['v_hi']:.1f}] boot {r['v_ref_boot']}", flush=True)
        print('   width', {kk: round(vv, 1) for kk, vv in r['width'].items()})
        for nm, v in r['ratios'].items():
            print(f"   {nm:20s} E={v['energy']:8.3f}  [{v['p16']:.3f},{v['p84']:.3f}] 95%[{v['p2p5']:.3f},"
                  f"{v['p97p5']:.3f}] snr_num {v['num_snr_boot']:.1f} snr_den {v['den_snr_boot']:.1f}")
        if 'Te_OIII' in r:
            print('   Te', r['Te_OIII'])
    # velocity by Dec band (model instrument)
    vb_ = []
    okP = np.isfinite(P_model)
    for lo, hi in ((12.0, 15.5), (15.5, 16.5), (16.5, 17.5), (17.5, 19.5), (19.5, 22.0), (22.0, 25.0), (25.0, 28.6)):
        m = okP & (dec >= lo) & (dec < hi)
        if m.sum() < 30:
            continue
        ubm, Bm = block_sums(C[m], W[m], P_model[m], okP[m], hp[m])
        Sm, eSm = reg(Bm)
        fo = fit_lines_1d(Sm, eSm, GROUPS_1D['OIII5007'], vgrid=np.arange(-60, 80.1, 0.5))
        vbb = []
        for b in range(150):
            cnt = np.bincount(RNG.integers(0, ubm.size, ubm.size), minlength=ubm.size).astype(float)
            Sbb = reg(Bm, cnt)[0]
            vbb.append(fit_lines_1d(Sbb, eSm, GROUPS_1D['OIII5007'], vgrid=np.arange(-40, 70.1, 1.0))['v'])
        vb_.append(dict(dec_lo=lo, dec_hi=hi, n=int(m.sum()), n_blocks=int(ubm.size), v=fo['v'], v_lo=fo['v_lo'],
                        v_hi=fo['v_hi'], amp=fo['amp'][0], amp_err=fo['err'][0],
                        v_boot_p16=float(np.percentile(vbb, 16)), v_boot_p84=float(np.percentile(vbb, 84)),
                        red_chi2=fo['red_chi2']))
    out['v_by_dec'] = vb_
    # velocity by survey/program and by nearest strand, all Dec and in the P1 band (Dec 15.5-17.5)
    prog = np.array([cat[i]['survey'] + '/' + cat[i]['program'] for i in ids])
    tr = json.loads((ROOT / 'data/derived/silk_trace.json').read_text())
    from spectra_desi_select import dist_to_polyline
    dE = dist_to_polyline(ra, dec, np.array([q['ra'] for q in tr['spine']]), np.array([q['dec'] for q in tr['spine']]))
    wpts = tr['strands']['W']['points']
    dW = dist_to_polyline(ra, dec, np.array([q['ra'] for q in wpts]), np.array([q['dec'] for q in wpts]))
    subsets = {}
    for lab, m0 in (('sv1', prog == 'sv1/dark'), ('sv1b', prog == 'sv1/bright'), ('main/dark', prog == 'main/dark'),
                    ('main/bright', prog == 'main/bright'), ('nearer_E', dE < dW), ('nearer_W', dW <= dE)):
        for band, mb in (('all', np.ones(len(ids), bool)), ('dec15.5-17.5', (dec >= 15.5) & (dec < 17.5))):
            m = okP & m0 & mb
            if m.sum() < 40:
                continue
            ubm, Bm = block_sums(C[m], W[m], P_model[m], okP[m], hp[m])
            Sm, eSm = reg(Bm)
            fo = fit_lines_1d(Sm, eSm, GROUPS_1D['OIII5007'], vgrid=np.arange(-60, 80.1, 0.5))
            subsets[f'{lab}|{band}'] = dict(n=int(m.sum()), n_blocks=int(ubm.size), v=fo['v'], v_lo=fo['v_lo'],
                                            v_hi=fo['v_hi'], amp=fo['amp'][0], amp_err=fo['err'][0])
    out['v_by_subset'] = subsets
    # simple stacks: on (pred > 90th pct), neg (< 10th pct), far (> 3 deg), ctrl (> 5 deg)
    for lab, m in (('on', okP & (P_model > np.nanpercentile(P_model, 90))),
                   ('neg', okP & (P_model < np.nanpercentile(P_model, 10))),
                   ('far', dany > 3), ('ctrl', dany > 5), ('patch', in_patch)):
        if m.sum() < 15:
            continue
        w = W[m]
        with np.errstate(invalid='ignore', divide='ignore'):
            s = np.sum(w * C[m], 0) / np.sum(w, 0)
            e = 1 / np.sqrt(np.sum(w, 0))
        full[f'stack_{lab}'], full[f'estack_{lab}'] = s, e
        st = {'n': int(m.sum()), 'mean_pred_R': float(np.nanmean(P_model[m])) if okP[m].any() else None}
        for g in ('OIII5007', 'Hb', 'HaNII', 'SII', 'sky5577', 'OI6300', 'OII3727', 'HeII4686'):
            st[g] = fit_fixed(s, e, GROUPS_1D[g], out['model_all']['v_ref_free']['v'])[0]
        out.setdefault('stacks', {})[lab] = st
    out['notes'] = ('Amplitudes: integrated line intensity per unit reference-line residual (photon units, R/R), from a '
                    'per-pixel weighted regression c = a + S P with intercept (INTERCEPT=True). '
                    'Energy ratio = photon ratio x lambda_den/lambda_num. Percentiles from healpix-block bootstrap. '
                    'Velocities barycentric (DESI) = heliocentric within 0.015 km/s. near_patch = within 1.2 deg of '
                    '(115.6, +20.5).')
    np.savez_compressed(D / 'spectra_desi_regspec_full.npz', **full)
    np.savez_compressed(D / 'spectra_desi_ratios_boot.npz',
                        **{k.replace('@', '_at_'): v for k, v in boots_all.items()})
    (D / 'spectra_desi_ratios.json').write_text(json.dumps(out, indent=1, default=float))
    print('v_by_dec', json.dumps(vb_, indent=0))
    print('v_by_subset', json.dumps(subsets, indent=0))
    print('stacks', json.dumps({k: {g: (v[g]['amp'] if isinstance(v[g], dict) else v[g]) for g in v}
                                for k, v in out['stacks'].items()}, indent=0, default=float))


if __name__ == '__main__':
    nb = int(sys.argv[sys.argv.index('--nboot') + 1]) if '--nboot' in sys.argv else 400
    main(nb)
