"""verify_xrayfuv: predicted GALEX FUV/NUV surface brightness (CU = photons s-1 cm-2 sr-1 A-1)
of shock / interface UV line emission for I([O III]5007) = 1, 2, 4 R.

Line tables were transcribed BY HAND from the ADS scans in data/raw/verify_xrayfuv_lit/ (energy
fluxes relative to H-beta; see each block's source). GALEX effective areas: SVO FPS
(GALEX/GALEX.FUV, GALEX/GALEX.NUV; GALEX team curves, Morrissey+2007). Extinction: Cardelli,
Clayton & Mathis (1989) R_V = 3.1, E(B-V) = 0.03 (SFD total column, an upper bound for a
~300 pc object) and 0.01. Line -> CU: a line of N photons s-1 cm-2 sr-1 produces the same count
rate as a flat-F_lambda continuum of  N * A(lam_line) * lam_eff / Int(A lam dlam)  CU (the GALEX
cps->F_lambda calibration is C = hc/Int(A lam dlam) for a flat F_lambda source; the curve-based C is
within 2% (FUV) / 8% (NUV) of the published 1.40e-15 / 2.06e-16).
H 2-photon continuum: spectral shape A(y) of Nussbaumer & Schmutz (1984) fit (y = nu/nu_Lya),
total energy from the model tables where given (R88, CR85: 14-21 x H-beta), otherwise 14-21 x H-beta.
Output: data/derived/verify_xrayfuv_uv_predictions.json
"""
import json, sys
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_xrayfuv_common import ROOT, DER

LIT = ROOT / 'data/raw/verify_xrayfuv_lit'
R2PH = 1e6 / (4 * np.pi)

# ---------------------------------------------------------------- GALEX curves
curves = {}
for b in ['FUV', 'NUV']:
    lam, A = np.loadtxt(LIT / f'svo_GALEX.{b}.dat').T
    curves[b] = (lam, A)
LAM_EFF = {'FUV': 1528.0, 'NUV': 2271.0}
C_PUB = {'FUV': 1.40e-15, 'NUV': 2.06e-16}
HC = 1.98645e-8


def Aeff(b, lam):
    l, a = curves[b]
    return np.interp(lam, l, a, left=0.0, right=0.0)


def cu_per_photon(b, lam):
    """CU produced by 1 photon s-1 cm-2 sr-1 of a line at lam."""
    l, a = curves[b]
    return Aeff(b, lam) * LAM_EFF[b] / np.trapezoid(a * l, l)


def ccm89(lam_A, Rv=3.1):
    x = 1e4 / np.atleast_1d(lam_A)
    a = np.zeros_like(x); b = np.zeros_like(x)
    o = (x >= 1.1) & (x < 3.3)
    y = x[o] - 1.82
    a[o] = 1 + 0.17699 * y - 0.50447 * y**2 - 0.02427 * y**3 + 0.72085 * y**4 + 0.01979 * y**5 - 0.77530 * y**6 + 0.32999 * y**7
    b[o] = 1.41338 * y + 2.28305 * y**2 + 1.07233 * y**3 - 5.38434 * y**4 - 0.62251 * y**5 + 5.30260 * y**6 - 2.09002 * y**7
    u_ = (x >= 3.3) & (x <= 8)
    xx = x[u_]
    Fa = np.where(xx >= 5.9, -0.04473 * (xx - 5.9)**2 - 0.009779 * (xx - 5.9)**3, 0)
    Fb = np.where(xx >= 5.9, 0.2130 * (xx - 5.9)**2 + 0.1207 * (xx - 5.9)**3, 0)
    a[u_] = 1.752 - 0.316 * xx - 0.104 / ((xx - 4.67)**2 + 0.341) + Fa
    b[u_] = -3.090 + 1.825 * xx + 1.206 / ((xx - 4.62)**2 + 0.263) + Fb
    return (a + b / Rv) * Rv   # A_lambda / E(B-V)


def twophoton_shape():
    """Photons per unit wavelength per decay for H 2s->1s (Nussbaumer & Schmutz 1984 fit)."""
    lam12 = 1215.67
    y = np.linspace(1e-4, 1 - 1e-4, 20001)
    w = y * (1 - y)
    Ay = 202.0 * (w * (1 - (4 * w)**0.8) + 0.88 * w**1.53 * (4 * w)**0.8)
    norm = np.trapezoid(Ay, y)          # = 2 x decay rate (two photons per decay)
    lam = lam12 / y
    n_lam = (Ay / (norm / 2)) * (lam12 / lam**2)   # photons per decay per A
    return lam[::-1], n_lam[::-1], norm


# ---------------------------------------------------------------- model line tables
# energy relative to H-beta = 100 (R79 Table 3/4; R88 Table 5 'Model'; CR85 Table 3) or = 1 (DS96 Tables 8/9).
# keys: wavelength (A) of the line or blend used for the bandpass; '5007' means 5007 alone,
# '4959+5007' means the doublet sum; '3727' the [O II] doublet.
MODELS = {}
# Raymond 1979 (ApJS 39, 1), models E-I: n0 = 10 cm^-3, abundance set A, complete steady flow.
r79 = {
    'E': dict(vs=100, O3=547, O2=260 + 322, HeII4686=3.10, C2_1336=289, C2_2324=324, C3_1910=657, C4=756,
              O2_2471=55, O3_1664=176, O4_1401=27, Mg2=34, Si4=107),
    'F': dict(vs=120, O3=572, O2=273 + 320, HeII4686=5.60, C2_1336=267, C2_2324=304, C3_1910=638, C4=1610,
              O2_2471=56, O3_1664=180, O4_1401=58, Mg2=37, Si4=100),
    'G': dict(vs=141, O3=596, O2=287 + 320, HeII4686=7.30, C2_1336=169, C2_2324=200, C3_1910=513, C4=2580,
              O2_2471=61, O3_1664=229, O4_1401=130, Mg2=34, Si4=51),
    'H': dict(vs=160, O3=495, O2=284 + 294, HeII4686=5.20, C2_1336=87, C2_2324=116, C3_1910=202, C4=1080,
              O2_2471=57, O3_1664=221, O4_1401=197, Mg2=26, Si4=37),
    'I': dict(vs=200, O3=373, O2=273 + 246, HeII4686=3.41, C2_1336=73, C2_2324=112, C3_1910=195, C4=612,
              O2_2471=49, O3_1664=157, O4_1401=159, Mg2=31, Si4=53),
}
for k, m in r79.items():
    MODELS[f'R79_{k}'] = dict(source='Raymond 1979 ApJS 39, 1, Tables 1/3/4 (complete, n0=10)', kind='complete',
                              vs=m['vs'], Hb=100.0, lines={
        '5007': m['O3'] * 0.749, '3727': m['O2'], 1336: m['C2_1336'], 1399: m['O4_1401'] + m['Si4'],
        1549: m['C4'], 1640: ('HeII', m['HeII4686']), 1664: m['O3_1664'], 1909: m['C3_1910'],
        2326: m['C2_2324'], 2471: m['O2_2471'], 2798: m['Mg2']}, twophoton=None)
# Dopita & Sutherland 1996 (ApJS 102, 161), Tables 8 (optical) and 9 (UV), shock only (no precursor),
# n = 1 cm^-3; suffix 10/11/12/14 = magnetic parameter B/sqrt(n) = 0/1/2/4.
ds = {  # vs, 5007, 3727+9, HeII4686, CII1335, OIV1402, NIV1486, CIV1550, OIII1666, NIII1750, SiIII1886, CIII1909, CII2326, NeIV2423, OII2470, MgII2800
    's15010': (150, 1.840, 5.556, 0.048, 1.761, 3.445, 0.370, 7.265, 1.158, 0.356, 0.367, 0.722, 0.516, 0.427, 0.208, 0.259),
    's20010': (200, 1.210, 3.739, 0.039, 1.791, 2.221, 0.125, 4.808, 0.739, 0.171, 0.321, 0.759, 0.468, 0.678, 0.137, 0.226),
    's15011': (150, 1.764, 5.707, 0.046, 1.689, 3.264, 0.353, 6.937, 1.103, 0.341, 0.355, 0.696, 0.505, 0.405, 0.198, 0.292),
    's20011': (200, 1.178, 4.841, 0.042, 1.738, 2.131, 0.120, 4.619, 0.623, 0.165, 0.269, 0.746, 0.490, 0.650, 0.136, 0.388),
    's15012': (150, 1.791, 6.174, 0.047, 1.740, 3.296, 0.351, 6.864, 1.113, 0.340, 0.363, 0.707, 0.534, 0.407, 0.201, 0.345),
    's20012': (200, 1.200, 5.656, 0.046, 1.738, 1.875, 0.120, 4.619, 0.716, 0.168, 0.275, 0.770, 0.513, 0.648, 0.140, 0.500),
    's15014': (150, 1.905, 6.807, 0.049, 1.796, 3.341, 0.364, 7.166, 1.152, 0.357, 0.385, 0.747, 0.564, 0.413, 0.206, 0.418),
    's20014': (200, 1.301, 6.715, 0.051, 1.772, 2.186, 0.124, 4.738, 0.746, 0.176, 0.295, 0.844, 0.545, 0.663, 0.149, 0.620),
}
for k, v in ds.items():
    vs, o3, o2, he2, c2, o4, n4, c4, o3s, n3, si3, c3, c2s, ne4, o2s, mg2 = v
    MODELS[f'DS96_{k}'] = dict(source='Dopita & Sutherland 1996 ApJS 102, 161, Tables 8/9 (shock only, n=1)',
                               kind='complete', vs=vs, Hb=1.0, lines={
        '5007': o3, '3727': o2, 1335: c2, 1402: o4, 1486: n4, 1549: c4, 1640: ('HeII', he2), 1666: o3s,
        1750: n3, 1886: si3, 1909: c3, 2326: c2s, 2423: ne4, 2470: o2s, 2800: mg2}, twophoton=None)
# Cox & Raymond 1985 (ApJ 298, 651) Table 3, runs 1-5 (complete flows, preionization families),
# normalised to recombination H-beta = 100; total H-beta given separately. He II 1640 case A (B).
cr = {  # run: vs, NV(out), OIV1400, SiIV1400, NIV1486, CIV, HeII1640(A), HeII1640(B), OIII]1666, NIII1750, SiIII1892, CIII1909, CII2325, MgII2800, OII3727, Hb, OIII5007+4959, 2g
    1: (200, 265, 35, 40, 573, 28, 100, 152, 34, 88, 171, 127, 126, 510, 129, 333, 2550),
    2: (186, 260, 33, 33, 560, 28, 99, 144, 31, 85, 169, 136, 137, 519, 102, 322, 1800),
    3: (140, 89, 150, 31, 932, 22, 97, 158, 60, 200, 577, 317, 152, 778, 157, 445, 3340),
    4: (130, 93, 146, 32, 981, 24, 106, 160, 61, 195, 582, 323, 160, 789, 130, 450, 2600),
    5: (119, 99, 121, 33, 1040, 27, 113, 155, 59, 168, 533, 308, 174, 758, 103, 440, 1850),
}
for k, v in cr.items():
    vs, o4, si4, n4, c4, heA, heB, o3s, n3, si3, c3, c2s, mg2, o2, hb, o3, tg = v
    MODELS[f'CR85_run{k}'] = dict(source='Cox & Raymond 1985 ApJ 298, 651, Table 3 (complete, preionization family)',
                                  kind='complete', vs=vs, Hb=float(hb), lines={
        '5007': o3 * 0.749, '3727': o2, 1400: o4 + si4, 1486: n4, 1549: c4, 1640: heA, 1666: o3s, 1750: n3,
        1892: si3, 1909: c3, 2325: c2s, 2800: mg2}, twophoton=float(tg), heII_caseB=float(heB))
# Raymond et al. 1988 (ApJ 324, 869) Table 5 'Model' columns: incomplete (truncated) shocks fitted to
# the Cygnus Loop spur; Table 7: SE vs=140, n0=4.4, log N_H(trunc)=18.20; M vs=120, n0=6.0, 18.25.
r88 = {  # CII1335, OIV]1400, NIV]1486, CIV1550, HeII1640, OIII]1664, NIII]1750, SiIII]1890, CIII]1909, CII]2325, NeIV2420, OII2470, OII3727, Hb, 4959, 5007, 2g
    'SE': (140, 486, 675, 222, 6800, 275, 537, 206, 185, 1420, 454, 219, 43, 1730, 100, 433, 1290, 1420),
    'M': (120, 717, 401, 129, 4330, 166, 427, 167, 261, 1570, 596, 140, 44, 1720, 100, 353, 1050, 1550),
}
for k, v in r88.items():
    vs, c2, o4, n4, c4, he2, o3s, n3, si3, c3, c2s, ne4, o2s, o2, hb, o3b, o3a, tg = v
    MODELS[f'R88_{k}'] = dict(source='Raymond et al. 1988 ApJ 324, 869, Tables 5 and 7 (incomplete/truncated)',
                              kind='incomplete', vs=vs, Hb=float(hb), lines={
        '5007': o3a, '3727': o2, 1335: c2, 1400: o4, 1486: n4, 1549: c4, 1640: he2, 1664: o3s, 1750: n3,
        1890: si3, 1909: c3, 2325: c2s, 2420: ne4, 2470: o2s}, twophoton=float(tg))

RESONANT = {1335, 1336, 1549, 1399, 1400, 1402, 2798, 2800}   # C II, C IV, Si IV (in 1400 blends), Mg II
HEII_1640_OVER_4686 = (5.6, 12.5)   # CR85 Table 3 case A / case B ratios (energy)
TWOG_OVER_HB = (14.2, 21.3)         # range of total 2-photon / total H-beta in R88 and CR85 tables


def model_cu(m, i5007_R, ebv, civ_supp=1.0, include_2g=True, he_ratio=None, tg_ratio=None):
    """CU in FUV and NUV for I(5007) = i5007_R Rayleigh; returns dict of totals and per-line."""
    L = m['lines']; e5007 = L['5007']
    lam_ph, n_lam, _ = twophoton_shape()
    out = {'FUV': {}, 'NUV': {}}
    for lam, val in L.items():
        if isinstance(lam, str):
            continue
        if isinstance(val, tuple):   # He II 1640 from 4686
            r = HEII_1640_OVER_4686[1] if he_ratio == 'B' else HEII_1640_OVER_4686[0] if he_ratio == 'A' else np.sqrt(np.prod(HEII_1640_OVER_4686))
            val = val[1] * r
        phot_ratio = (val / e5007) * (lam / 5007.0)        # photons(line)/photons(5007)
        nph = i5007_R * R2PH * phot_ratio
        if lam in RESONANT and abs(lam - 1549) < 2:
            nph /= civ_supp
        att = 10 ** (-0.4 * ccm89(lam)[0] * ebv)
        for b in ['FUV', 'NUV']:
            out[b][lam] = nph * cu_per_photon(b, lam) * att
    if include_2g:
        tg = m['twophoton'] if m.get('twophoton') else (tg_ratio if tg_ratio else np.sqrt(np.prod(TWOG_OVER_HB))) * m['Hb']
        # energy ratio 2g/5007 -> number of decays: E_2g / (h nu_Lya)
        e_ratio = tg / e5007
        n_decay = i5007_R * R2PH * e_ratio * (1215.67 / 5007.0)   # decays s-1 cm-2 sr-1 (photon-energy scaling)
        att = 10 ** (-0.4 * ccm89(lam_ph) * ebv)
        for b in ['FUV', 'NUV']:
            out[b]['2g'] = float(np.trapezoid(n_decay * n_lam * cu_per_photon(b, lam_ph) * att, lam_ph))
    res = {}
    for b in ['FUV', 'NUV']:
        tot = sum(out[b].values())
        nores = sum(v for k, v in out[b].items() if not (k in RESONANT))
        res[b] = {'total': float(tot), 'non_resonant': float(nores),
                  'lines': {str(k): float(v) for k, v in out[b].items()}}
    return res


# ---------------------------------------------------------------- conductive interfaces
# Borkowski, Balbus & Fristrom 1990 (ApJ 355, 501) Table 1: photon number fluxes (cm^-2 s^-1) of the
# standard non-magnetized evaporating/condensing front (T_h0 = 7.5e5 K, n_h0 = 2.3e-3 cm^-3) at
# t = 2.8e5, 2.2e6, 7.4e6 yr. [O III] 5007 is not tabulated; it is bounded from O III] 1663 with the
# smallest photon ratio N(1663)/N(5007) found in the shock tables above (lower-T gas gives a lower
# ratio, so this gives an UPPER bound on 5007 and a LOWER bound on FUV per R of 5007).
# Face-on intensity I = F / (4 pi) (F taken as emission into 4 pi per unit front area; x2 if one-sided).
BBF90 = {'2.8e5yr': {1640: 101, 2326: 237, 1909: 414, 1549: 1110, 1486: 35, 1750: 17, 1664: 46, 1407: 95, 1892: 118, 1398: 104, 2141: 39},
         '2.2e6yr': {1640: 13, 2326: 478, 1909: 680, 1549: 525, 1486: 53, 1750: 41, 1664: 123, 1407: 226, 1892: 182, 1398: 100, 2141: 71},
         '7.4e6yr': {1640: 29, 2326: 602, 1909: 554, 1549: 621, 1486: 50, 1750: 53, 1664: 210, 1407: 269, 1892: 236, 1398: 64, 2141: 100}}
SLAVIN89_CIV = 278.0   # photons cm-2 s-1 sr-1, Slavin 1989 ApJ 346, 718, Table 2 (Local Cloud boundary model)


def interfaces(ebv=0.03):
    # smallest N(1663)/N(5007) photon ratio among the shock models
    rmin = min((m['lines'].get(1664, m['lines'].get(1666, np.nan)) / m['lines']['5007']) * (1664 / 5007.0)
               for m in MODELS.values() if (1664 in m['lines'] or 1666 in m['lines']))
    out = {'min_photon_ratio_1663_over_5007_in_shock_tables': float(rmin)}
    for ep, L in BBF90.items():
        I = {k: v / (4 * np.pi) for k, v in L.items()}
        i5007_max = I[1664] / rmin
        rec = {'I5007_faceon_max_ph': float(i5007_max), 'I5007_faceon_max_R': float(i5007_max / R2PH),
               'limb_brightening_needed_for_1R': float(R2PH / i5007_max),
               'CIV_faceon_ph': float(I[1549]), 'CIV_over_5007_photons_min': float(I[1549] / i5007_max)}
        for b in ['FUV', 'NUV']:
            tot = 0.0; nores = 0.0
            for lam, val in I.items():
                c = val * cu_per_photon(b, lam) * 10 ** (-0.4 * ccm89(lam)[0] * ebv) / i5007_max * R2PH
                tot += c
                if lam not in (1549, 1398):
                    nores += c
            rec[f'{b}_CU_per_R5007_min'] = float(tot); rec[f'{b}_nonres_CU_per_R5007_min'] = float(nores)
        out[f'BBF90_{ep}'] = rec
    out['Slavin89_CIV_faceon_R'] = SLAVIN89_CIV / R2PH
    return out


def main():
    lam_ph, n_lam, norm = twophoton_shape()
    chk = {'twophoton_norm_int_A_dy': float(norm), 'twophoton_photons_per_decay': float(np.trapezoid(n_lam, lam_ph)),
           'twophoton_energy_per_decay_over_hnuLya': float(np.trapezoid(n_lam * 1215.67 / lam_ph, lam_ph))}
    for b in ['FUV', 'NUV']:
        l, a = curves[b]
        chk[f'{b}_C_curve'] = float(HC / np.trapezoid(a * l, l)); chk[f'{b}_C_pub'] = C_PUB[b]
        chk[f'{b}_Aext_per_EBV_at_lameff'] = float(ccm89(LAM_EFF[b])[0])
        for lam in [1335, 1400, 1486, 1549, 1640, 1664, 1750, 1892, 1909, 2326, 2423, 2470, 2800]:
            chk[f'{b}_A({lam})'] = float(Aeff(b, lam))
            chk[f'{b}_Weff({lam})_A'] = float(1 / cu_per_photon(b, lam)) if Aeff(b, lam) > 0 else None
    print(json.dumps(chk, indent=0)[:2500])
    out = {'checks': chk, 'models': {}, 'notes': __doc__}
    for name, m in MODELS.items():
        L = m['lines']
        o3hb = L['5007'] / m['Hb']; o2o3 = L['3727'] / L['5007']
        rec = {'source': m['source'], 'kind': m['kind'], 'vs': m['vs'], '5007/Hb': o3hb, '3727/5007': o2o3,
               '4959+5007/Hb': L['5007'] / 0.749 / m['Hb']}
        pred = {}
        for ebv in [0.0, 0.01, 0.03]:
            for supp in [1.0, 3.5, 20.0]:
                p = model_cu(m, 1.0, ebv, civ_supp=supp)
                pred[f'ebv{ebv}_civsupp{supp}'] = p
        # ranges for He II and 2-photon assumptions (per 1 R, E(B-V)=0.03, no suppression)
        lo = model_cu(m, 1.0, 0.03, he_ratio='A', tg_ratio=TWOG_OVER_HB[0])
        hi = model_cu(m, 1.0, 0.03, he_ratio='B', tg_ratio=TWOG_OVER_HB[1])
        rec['per_R'] = pred
        rec['per_R_assumption_range_ebv0.03'] = {b: [lo[b]['non_resonant'], hi[b]['non_resonant']] for b in ['FUV', 'NUV']}
        out['models'][name] = rec
        p = pred['ebv0.03_civsupp1.0']
        print(f"{name:12s} vs={m['vs']:3d} 5007/Hb={o3hb:5.1f} 3727/5007={o2o3:4.2f} | per 1 R(5007), E(B-V)=0.03: "
              f"FUV tot {p['FUV']['total']:6.0f} non-res {p['FUV']['non_resonant']:5.0f} (CIV {p['FUV']['lines'].get('1549',0):5.0f}, 2g {p['FUV']['lines']['2g']:4.0f}) | "
              f"NUV tot {p['NUV']['total']:5.0f} non-res {p['NUV']['non_resonant']:5.0f} (2g {p['NUV']['lines']['2g']:4.0f})")
    out['interfaces'] = interfaces()
    print(json.dumps(out['interfaces'], indent=1))
    (DER / 'verify_xrayfuv_uv_predictions.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
