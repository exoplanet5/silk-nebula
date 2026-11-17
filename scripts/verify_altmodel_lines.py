"""verify_altmodel_lines: line-emissivity tables, spectrum synthesis from ion fractions, GALEX band
conversion and extinction, shared by the verify_altmodel_* model scripts.

Atomic data: CHIANTI v7.1.4 (see verify_altmodel_atomic.py) for every collisionally excited line,
summed over all transitions of the ion that fall in a vacuum-wavelength window (multiplets).
Hydrogen: case-B recombination Hbeta/Halpha (Draine 2011 fits) + collisional excitation of H0 from
the CHIANTI H I model (n <= 5).  He II 4686 recombination = 12 x Hbeta emissivity per He++ (+-20 %);
He II 1640 recombination = 7.0 x 4686 (energy; case-B order, +-30 %) plus CHIANTI collisional 1640.
H two-photon continuum: collisional 2s population x A(2s-1s) from CHIANTI + recombination
(0.32 decays per case-B recombination, Osterbrock & Ferland 2006 order); spectral shape of
Nussbaumer & Schmutz (1984).
GALEX FUV/NUV effective-area curves: SVO FPS files already on disk (data/raw/verify_xrayfuv_lit/,
downloaded by the xrayfuv agent; read here as plain text).  CU = photons s^-1 cm^-2 sr^-1 A^-1; a line
of photon intensity N gives N * A(lambda) / int A dlambda CU (flat-photon-spectrum equivalent).
Extinction: CCM89 R_V = 3.1 (own implementation), E(B-V) = 0.03 default, applied differentially to
5007 (observed [O III] is also attenuated).
"""
import os
import numpy as np
import verify_altmodel_atomic as A

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
TAB = os.path.join(ROOT, "data/derived/verify_altmodel_linetab.npz")
SVO = os.path.join(ROOT, "data/raw/verify_xrayfuv_lit")

# name: (Z, ion stage (1 = neutral), vacuum wavelength window [A])
WINDOWS = {
    "OIII5007": (8, 3, 5007.5, 5009.0), "OIII4959": (8, 3, 4959.5, 4961.5), "OIII4363": (8, 3, 4363.5, 4365.5),
    "OII3726": (8, 2, 3726.5, 3727.6), "OII3729": (8, 2, 3729.3, 3730.5), "OII7325": (8, 2, 7318.0, 7335.0),
    "OII2470": (8, 2, 2469.0, 2472.0),
    "NII6583": (7, 2, 6584.5, 6586.0), "NII6548": (7, 2, 6549.0, 6551.0), "NII5755": (7, 2, 5755.5, 5757.0),
    "SII6716": (16, 2, 6717.5, 6719.5), "SII6731": (16, 2, 6732.0, 6734.0),
    "NeIII3869": (10, 3, 3869.0, 3870.5), "SIII9069": (16, 3, 9068.0, 9073.0), "ArIII7136": (18, 3, 7136.5, 7139.0),
    "CIV1549": (6, 4, 1547.0, 1552.0), "OIII1663": (8, 3, 1659.0, 1668.0), "NIII1750": (7, 3, 1745.0, 1756.0),
    "NIV1486": (7, 4, 1482.0, 1488.0), "OIV1400": (8, 4, 1392.0, 1409.0), "SiIV1400": (14, 4, 1392.0, 1404.0),
    "CIII1909": (6, 3, 1906.0, 1910.0), "CII2326": (6, 2, 2323.0, 2330.0), "NeIV2423": (10, 4, 2421.0, 2426.0),
    "MgII2800": (12, 2, 2795.0, 2805.0), "SiIII1892": (14, 3, 1882.0, 1893.0), "NV1240": (7, 5, 1238.0, 1244.0),
    "OVI1035": (8, 6, 1031.0, 1038.0), "CIII977": (6, 3, 976.0, 978.0), "HeII1640c": (2, 2, 1639.0, 1642.0),
}
LAM = {"OIII5007": 5008.2, "OIII4959": 4960.3, "OIII4363": 4364.4, "OII3726": 3727.1, "OII3729": 3729.9,
       "OII7325": 7322.0, "OII2470": 2471.0, "NII6583": 6585.3, "NII6548": 6549.9, "NII5755": 5756.2,
       "SII6716": 6718.3, "SII6731": 6732.7, "NeIII3869": 3869.9, "SIII9069": 9071.1, "ArIII7136": 7137.8,
       "CIV1549": 1549.5, "OIII1663": 1664.0, "NIII1750": 1750.0, "NIV1486": 1486.5, "OIV1400": 1401.0,
       "SiIV1400": 1397.0, "CIII1909": 1908.0, "CII2326": 2326.0, "NeIV2423": 2423.0, "MgII2800": 2800.0,
       "SiIII1892": 1892.0, "NV1240": 1240.0, "OVI1035": 1034.0, "CIII977": 977.0, "HeII1640": 1640.4,
       "HeII4686": 4687.0, "Hb": 4862.7, "Ha": 6564.6}
TGRID = np.logspace(np.log10(5e3), np.log10(3e6), 150)
NEGRID = np.array([1.0, 100.0])
LYA_A = 1215.67


def _window_emis(I, T, ne, lo, hi):
    p = I.populations(T, ne)
    m = (I.wl >= lo) & (I.wl <= hi) & (I.A > 0)
    dE = (I.E[:, None] - I.E[None, :]) * A.C_L * A.H_PL
    return float((p[:, None] * I.A * dE * m).sum() / ne)


def build_table(force=False):
    if os.path.exists(TAB) and not force:
        d = np.load(TAB)
        return {k: d[k] for k in d.files}
    out = {"T": TGRID, "ne": NEGRID}
    byion = {}
    for name, (Z, ion, lo, hi) in WINDOWS.items():
        byion.setdefault((Z, ion), []).append((name, lo, hi))
    for (Z, ion), wins in byion.items():
        I = A.get_ion(Z, ion)
        dE = (I.E[:, None] - I.E[None, :]) * A.C_L * A.H_PL
        arrs = {nm: np.zeros((len(NEGRID), len(TGRID))) for nm, _, _ in wins}
        masks = {nm: (I.wl >= lo) & (I.wl <= hi) & (I.A > 0) for nm, lo, hi in wins}
        for a, ne in enumerate(NEGRID):
            for b, T in enumerate(TGRID):
                p = I.populations(T, ne)
                P = p[:, None] * I.A * dE
                for nm in arrs:
                    arrs[nm][a, b] = float((P * masks[nm]).sum() / ne)
        for nm, arr in arrs.items():
            out["e_" + nm] = arr
            print("table", nm, I.n, "levels", "e(1e4, 3e4, 1e5) =",
                  np.interp(np.log([1e4, 3e4, 1e5]), np.log(TGRID), arr[0]), flush=True)
    # H I collisional Balmer + 2s two-photon decays (per H0 per n_e)
    H = A.get_ion(1, 1)
    hb = np.zeros((len(NEGRID), len(TGRID)))
    ha = np.zeros_like(hb)
    n2s = np.zeros_like(hb)
    i2s = 1
    for a, ne in enumerate(NEGRID):
        for b, T in enumerate(TGRID):
            c = A.h_collisional_balmer(T, ne)
            hb[a, b], ha[a, b] = c["Hb"], c["Ha"]
            p = H.populations(T, ne)
            n2s[a, b] = p[i2s] * H.A[i2s, 0] / ne     # two-photon decays per H0 per n_e (s^-1 cm^3)
    out["e_Hb_coll"], out["e_Ha_coll"], out["r_2g_coll"] = hb, ha, n2s
    np.savez(TAB, **out)
    return out


_T = None


def tab():
    global _T
    if _T is None:
        _T = build_table()
    return _T


def e_line(name, T, ne=1.0):
    """per-ion emissivity (erg cm^3 s^-1) of window line `name` at T (array ok); log-interp in T, and
    linear interpolation in log n_e between the two tabulated densities (clipped)."""
    t = tab()
    T = np.atleast_1d(np.asarray(T, float))
    lt = np.log(np.clip(T, TGRID[0], TGRID[-1]))
    arr = t["e_" + name]
    v = [np.exp(np.interp(lt, np.log(TGRID), np.log(np.maximum(arr[k], 1e-300)))) for k in range(2)]
    w = np.clip(np.log10(ne) / 2.0, 0, 1)
    return np.exp((1 - w) * np.log(v[0]) + w * np.log(v[1]))


# ------------------------------------------------------------------------------------------------
def synth(T, fr, ne_nH=1.2, ne=1.0, ab=None, caseB=True):
    """line emissivities per n_H^2 (erg cm^3 s^-1) for temperature T and ion fractions fr:
    dict Z -> array of stage fractions (index 0 = neutral). Returns dict name -> value (arrays ok) and
    the two-photon decay rate per n_H^2 ('2g', s^-1 cm^3)."""
    ab = ab or A.ABUND_A09
    T = np.atleast_1d(np.asarray(T, float))
    out = {}
    for name, (Z, ion, lo, hi) in WINDOWS.items():
        if Z not in fr:
            continue
        out[name] = ne_nH * ab[Z] * fr[Z][ion - 1] * e_line(name, T, ne)
    xH0, xHp = fr[1][0], fr[1][1]
    hbr = A.j_hbeta(T) * (1.0 if caseB else 0.67)
    out["Hb_rec"] = ne_nH * xHp * hbr
    out["Ha_rec"] = ne_nH * xHp * A.j_halpha(T) * (1.0 if caseB else 0.67 * 1.04)
    out["Hb_coll"] = ne_nH * xH0 * _e_h("e_Hb_coll", T, ne)
    out["Ha_coll"] = ne_nH * xH0 * _e_h("e_Ha_coll", T, ne)
    out["Hb"], out["Ha"] = out["Hb_rec"] + out["Hb_coll"], out["Ha_rec"] + out["Ha_coll"]
    if 2 in fr:
        he3 = fr[2][2]
        out["HeII4686"] = ne_nH * ab[2] * he3 * 12.0 * A.j_hbeta(T)
        out["HeII1640"] = 7.0 * out["HeII4686"] + out.get("HeII1640c", 0.0)
    rec = ne_nH * xHp * A.alpha_B(T)
    out["2g"] = 0.32 * rec + ne_nH * xH0 * _e_h("r_2g_coll", T, ne)
    return out


def _e_h(key, T, ne):
    t = tab()
    lt = np.log(np.clip(np.atleast_1d(T), TGRID[0], TGRID[-1]))
    v = [np.exp(np.interp(lt, np.log(TGRID), np.log(np.maximum(t[key][k], 1e-300)))) for k in range(2)]
    w = np.clip(np.log10(ne) / 2.0, 0, 1)
    return np.exp((1 - w) * np.log(v[0]) + w * np.log(v[1]))


# ------------------------------------------------------------------------------------------------
def ccm89(lam_A, Rv=3.1):
    x = 1e4 / np.asarray(lam_A, float)
    a = np.zeros_like(x)
    b = np.zeros_like(x)
    ir = x < 1.1
    a[ir], b[ir] = 0.574 * x[ir] ** 1.61, -0.527 * x[ir] ** 1.61
    op = (x >= 1.1) & (x < 3.3)
    y = x[op] - 1.82
    a[op] = 1 + 0.17699 * y - 0.50447 * y**2 - 0.02427 * y**3 + 0.72085 * y**4 + 0.01979 * y**5 - 0.77530 * y**6 + 0.32999 * y**7
    b[op] = 1.41338 * y + 2.28305 * y**2 + 1.07233 * y**3 - 5.38434 * y**4 - 0.62251 * y**5 + 5.30260 * y**6 - 2.09002 * y**7
    uv = (x >= 3.3) & (x <= 8)
    xx = x[uv]
    Fa = np.where(xx >= 5.9, -0.04473 * (xx - 5.9) ** 2 - 0.009779 * (xx - 5.9) ** 3, 0)
    Fb = np.where(xx >= 5.9, 0.2130 * (xx - 5.9) ** 2 + 0.1207 * (xx - 5.9) ** 3, 0)
    a[uv] = 1.752 - 0.316 * xx - 0.104 / ((xx - 4.67) ** 2 + 0.341) + Fa
    b[uv] = -3.090 + 1.825 * xx + 1.206 / ((xx - 4.62) ** 2 + 0.263) + Fb
    return (a + b / Rv) * Rv   # A_lambda / E(B-V)


def _read_svo(band):
    w, r = [], []
    with open(os.path.join(SVO, f"svo_GALEX.{band}.dat")) as f:
        for l in f:
            p = l.split()
            if len(p) >= 2:
                try:
                    w.append(float(p[0])); r.append(float(p[1]))
                except ValueError:
                    pass
    return np.array(w), np.array(r)


_BANDS = {}


def band(bname):
    if bname not in _BANDS:
        w, r = _read_svo(bname)
        _BANDS[bname] = (w, r, np.trapezoid(r, w))
    return _BANDS[bname]


def W_eff(bname, lam):
    w, r, tot = band(bname)
    a = np.interp(lam, w, r, left=0, right=0)
    return tot / a if a > 0 else np.inf


def twophoton_shape():
    """photons per decay per unit y (y = nu/nu_Lya), Nussbaumer & Schmutz (1984) fit, normalised."""
    y = np.linspace(1e-4, 1 - 1e-4, 4000)
    C, al, be, ga = 202.0, 0.88, 1.53, 0.8
    z = y * (1 - y)
    Ay = C * (z * (1 - (4 * z) ** ga) + al * z ** be * (4 * z) ** ga)
    Ay = 2.0 * Ay / np.trapezoid(Ay, y)
    return y, Ay


def uv_CU_per_R5007(L, ebv=0.03, lines_fuv=("CIV1549", "OIII1663", "NIII1750", "NIV1486", "OIV1400", "SiIV1400",
                                              "HeII1640", "NV1240"),
                    lines_nuv=("CIII1909", "CII2326", "NeIV2423", "MgII2800", "SiIII1892", "OII2470", "HeII1640",
                               "OIII1663", "NIII1750")):
    """GALEX FUV and NUV surface brightness (CU) per 1 R of observed [O III]5007 for a spectrum L
    (energy emissivities incl. 5007 and the '2g' decay rate). Returns dict with totals and parts."""
    E5007 = A.H_PL * A.C_L / (LAM["OIII5007"] * 1e-8)
    ph5007 = L["OIII5007"] / E5007
    R2ph = 1e6 / (4 * np.pi)
    k5007 = ccm89([LAM["OIII5007"]])[0]
    res = {}
    for bname, lines in (("FUV", lines_fuv), ("NUV", lines_nuv)):
        parts = {}
        w, r, tot = band(bname)
        for ln in lines:
            if ln not in L:
                continue
            lam = LAM[ln]
            a = np.interp(lam, w, r, left=0, right=0)
            if a <= 0:
                continue
            ph = (L[ln] / (A.H_PL * A.C_L / (lam * 1e-8))) / ph5007        # photons per 5007 photon
            att = 10 ** (-0.4 * ebv * (ccm89([lam])[0] - k5007))
            parts[ln] = float(ph * R2ph * a / tot * att)
        # two-photon continuum
        y, Ay = twophoton_shape()
        lam2 = LYA_A / y
        a2 = np.interp(lam2, w, r, left=0, right=0)
        att2 = 10 ** (-0.4 * ebv * (ccm89(np.clip(lam2, 1000, 3.3e4)) - k5007))
        frac = np.trapezoid(Ay * a2 * att2, y) / tot     # (photons per decay) weighted by response / int A
        parts["2g"] = float(L["2g"] / ph5007 * R2ph * frac)
        res[bname] = dict(total=float(sum(parts.values())), parts=parts,
                          non_resonant=float(sum(v for k, v in parts.items() if k not in ("CIV1549", "SiIV1400", "NV1240", "MgII2800"))))
    return res


def ratios(L):
    o3 = L["OIII5007"]
    r = {"OIII_Hb": o3 / L["Hb"], "OII_OIII": (L["OII3726"] + L["OII3729"]) / o3,
         "O4363_5007": L["OIII4363"] / o3, "HeII_OIII": L.get("HeII4686", 0 * o3) / o3,
         "NII_Ha": L["NII6583"] / L["Ha"], "SII_Ha": (L["SII6716"] + L["SII6731"]) / L["Ha"],
         "OIII_Ha_phot": o3 / L["Ha"] * LAM["OIII5007"] / LAM["Ha"],
         "NeIII_OIII": L["NeIII3869"] / o3, "OII7325_OIII": L["OII7325"] / o3, "OII_7325_3727": L["OII7325"] / (L["OII3726"] + L["OII3729"]),
         "SIII9069_OIII": L["SIII9069"] / o3, "ArIII_OIII": L.get("ArIII7136", 0 * o3) / o3,
         "NII5755_6583": L["NII5755"] / L["NII6583"], "OII3729_3726": L["OII3729"] / L["OII3726"],
         "CIV_OIII_phot": L["CIV1549"] / o3 * LAM["CIV1549"] / LAM["OIII5007"],
         "OVI_OIII_energy": L["OVI1035"] / o3, "Hb_coll_frac": L["Hb_coll"] / L["Hb"]}
    return r
