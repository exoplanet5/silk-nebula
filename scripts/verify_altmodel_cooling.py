"""verify_altmodel_cooling: collisional (non-photoionised) alternatives to a shock for the Silk [O III].

Model families (all my own calculations; ion fractions from the published Gnat & Sternberg 2007 tables
or from my own CIE/NEI solver; emissivities from verify_altmodel_lines / CHIANTI v7.1.4):
  1. NEI radiative cooling of hot gas (GS07 isochoric and isobaric tracks, Z = 1 and 2, start 5e6 K):
     integrated line energy per H nucleus  E_l = (c_p/k) k * int L_l(T)/(n_H^2) * (n_tot/n_H)/((n_e/n_H) Lambda) dT,
     c_p = 5/2 (isobaric) or 3/2 (isochoric); Lambda = GS07 cooling efficiency (rate = n_e n_H Lambda).
     The same track integrated only below T_start (e.g. 1.6e6 K = Monogem interior) and down to
     T_stop = 1e4 (full), 2e4, 3e4, 5e4 K (truncated).  A recombination 'tail' below 1e4 K is bracketed
     by adding x(H+,1e4 K) case-B recombinations (0.119 Hbeta photons each).
     The isobaric track is also the proxy for radiative turbulent mixing layers (Ji, Oh & Masterson 2019
     compare their TML columns with GS07 isobaric cooling; they find TML emission weighted to
     intermediate T even more strongly) and for radiative condensation fronts.
  2. Warm CIE layer at a single T (gas held at T by an unspecified local heating): GS07 CIE fractions,
     checked against my own CIE balance (Dere 2007 CI + CHIANTI RR/DR).
  3. Static radiative conduction front (conduction balancing radiation, no flow, CIE ions):
     F dF/dT = kappa0 T^5/2 n_e n_H Lambda,  dx = kappa0 T^5/2 dT / F, F(T_c = 1e4 K) = 0,
     kappa0 = 6e-7 cgs (Spitzer); thermal pressure P/k = 1.5e4 (Monogem interior) and 2e6 K cm^-3
     (the [O II] 3729/3726 density x T).  Face-on intensities scale as P.
  4. Steady planar evaporation front (Cowie & McKee 1977 classical, no radiation in the energy balance):
     dT/dx = (5/2) J k T^-3/2 / (mu m_H kappa0), v = J k T/(mu m_H P); ions integrated with my own
     NEI solver from T = 1e4 K (WIM-like start: H+, O+ ...) to 1.5e6 K; J = classical mass flux for
     cloud radius 1-60 pc, capped by the saturated flux.
Outputs: data/derived/verify_altmodel_cooling.json, figures/verify_altmodel_cooling.png
Run: ~/.venvs/astro313/bin/python -I scripts/verify_altmodel_cooling.py
"""
import json
import os
import sys

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, os.path.dirname(__file__))
import verify_altmodel_atomic as A  # noqa: E402
import verify_altmodel_lines as LN  # noqa: E402
from verify_altmodel_common import sheet_G  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
GS = os.path.join(ROOT, "data/raw/verify_altmodel/gs07")
OUTJ = os.path.join(ROOT, "data/derived/verify_altmodel_cooling.json")
FIG = os.path.join(ROOT, "figures/verify_altmodel_cooling.png")

K_B, M_H = A.K_B, 1.6726e-24
# GS07 Table 1 (Asplund+2005, enhanced Ne); He from Asplund+2005
AB_GS = {1: 1.0, 2: 0.085, 6: 2.45e-4, 7: 6.03e-5, 8: 4.57e-4, 10: 1.95e-4, 12: 3.39e-5, 14: 3.24e-5, 16: 1.38e-5}
ZS = {"H": 1, "He": 2, "C": 6, "N": 7, "O": 8, "Ne": 10, "Mg": 12, "Si": 14, "S": 16, "Fe": 26}
KAPPA0 = 6.0e-7
I_OBS = dict(fid=2.0, lo=1.4, hi=4.0)          # R, [O III]5007 (verify_desi / verify_sdss / NSNS)


# ------------------------------------------------------------------------------------------------
def read_gs_ion(fn):
    lab = []
    rows = []
    with open(os.path.join(GS, fn)) as f:
        L = f.readlines()
    import re
    pat = re.compile(r"^\s*\d+-\s*\d+\s+A\d+\s+\S+\s+(\S+)\s+Fractional")
    for l in L:
        mm = pat.match(l)
        if mm:
            lab.append(mm.group(1))
    for l in L:
        p = l.split()
        if len(p) == len(lab) + 1:
            try:
                rows.append([float(v) for v in p])
            except ValueError:
                pass
    arr = np.array(rows)
    T = arr[:, 0]
    fr = {}
    for el, Z in ZS.items():
        cols = [k for k, s in enumerate(lab) if s.split("{")[0] == el]
        if cols:
            fr[Z] = arr[:, [c + 1 for c in cols]].T     # [stage, T]
    o = np.argsort(T)
    return T[o], {Z: v[:, o] for Z, v in fr.items()}


def read_gs_cool(fn, col):
    rows = []
    with open(os.path.join(GS, fn)) as f:
        for l in f:
            p = l.split()
            if len(p) == 6:
                try:
                    rows.append([float(v) for v in p])
                except ValueError:
                    pass
    a = np.array(rows)
    o = np.argsort(a[:, 0])
    return a[o, 0], a[o, col]


def ne_nH_of(fr, ab):
    return sum(ab[Z] * (np.arange(fr[Z].shape[0])[:, None] * fr[Z]).sum(0) for Z in fr if Z in ab)


def ntot_nH_of(fr, ab):
    return sum(ab[Z] for Z in fr if Z in ab) + ne_nH_of(fr, ab)


def summarize(L, extra=None):
    r = {k: float(v) for k, v in LN.ratios(L).items()}
    uv = LN.uv_CU_per_R5007(L)
    r["FUV_CU_per_R"] = uv["FUV"]["total"]
    r["FUV_nonres_CU_per_R"] = uv["FUV"]["non_resonant"]
    r["NUV_CU_per_R"] = uv["NUV"]["total"]
    r["FUV_parts"] = uv["FUV"]["parts"]
    r["NUV_parts"] = uv["NUV"]["parts"]
    if extra:
        r.update(extra)
    return r


def integrate_track(T, fr, Lam, cp, ab, T_hi, T_lo, tail=False, ne_dens=1.0):
    """integrated line energies per H (erg) for a cooling track between T_hi and T_lo."""
    m = (T <= T_hi * 1.0001) & (T >= T_lo * 0.9999)
    Tm = T[m]
    frm = {Z: v[:, m] for Z, v in fr.items()}
    nenh = ne_nH_of(frm, ab)
    ntnh = ntot_nH_of(frm, ab)
    Lm = LN.synth(Tm, frm, ne_nH=nenh, ne=ne_dens, ab=ab)
    lam = np.exp(np.interp(np.log(Tm), np.log(Lam[0]), np.log(Lam[1])))
    w = cp * K_B * ntnh / (nenh * lam)          # n_H dt / dT   (s cm^-3 K^-1)
    E = {k: float(np.trapezoid(v * w, Tm)) for k, v in Lm.items() if np.ndim(v)}
    E["total_radiated"] = float(np.trapezoid(nenh * lam * w, Tm))
    if tail:
        xHp = frm[1][1][np.argmin(Tm)]
        ehb = A.H_PL * A.C_L / 4.8627e-5
        E["Hb"] += xHp * 0.119 * ehb
        E["Ha"] += xHp * 0.119 * 2.86 * ehb      # case-B Halpha/Hbeta energy ratio 2.86
        E["2g"] += xHp * 0.32
        E["tail_xHp_1e4"] = float(xHp)
    return E


def required_flux(E, I_R=2.0, G=None):
    """energy flux and H-nucleus flux through the cooling layer needed for I_R of 5007 viewed with
    limb-brightening G (path / layer thickness)."""
    if G is None:
        G, _ = sheet_G()
    ph5007 = 4 * np.pi * I_R * 1e6 / (4 * np.pi) / G          # face-on 5007 photon flux (cm^-2 s^-1)
    e5007 = A.H_PL * A.C_L / 5.0082e-5
    jH = ph5007 * e5007 / E["OIII5007"]                        # H nuclei cm^-2 s^-1
    return dict(G=float(G), face_on_5007_erg=float(ph5007 * e5007), H_flux=float(jH),
                energy_flux=float(jH * E["total_radiated"]), f5007=float(E["OIII5007"] / E["total_radiated"]))


# ------------------------------------------------------------------------------------------------
# own CIE / NEI solver
ELEM = [1, 2, 6, 7, 8, 10, 14, 16]


def rate_matrices(T):
    out = {}
    for Z in ELEM:
        n = Z + 1
        M = np.zeros((n, n))
        for i in range(Z):            # stage i (0-based, i = 0 neutral) -> i+1
            c = float(A.ci_rate(Z, i + 1, T))
            a = float(A.rec_rate(Z, i + 2, T)[0])
            M[i + 1, i] += c
            M[i, i] -= c
            M[i, i + 1] += a
            M[i + 1, i + 1] -= a
        out[Z] = M
    return out


def own_cie(T):
    fr = {Z: np.zeros((Z + 1, len(T))) for Z in ELEM}
    for k, t in enumerate(T):
        for Z in ELEM:
            f = A.cie_fractions(Z, [t])[0]
            fr[Z][:, k] = f
    return fr


def evaporation_front(J, P_k, T0=1e4, T1=1.5e6, start="wim", mu=0.6, nstep=700):
    """steady classical evaporation; returns face-on line intensities (R) and columns."""
    P = P_k * K_B
    lT = np.linspace(np.log(T0), np.log(T1), nstep)
    T = np.exp(lT)
    # initial ionisation
    x = {}
    for Z in ELEM:
        f = np.zeros(Z + 1)
        if start == "wim":
            f[1] = 1.0            # H+, He+, C+, N+, O+, Ne+, Si+, S+
            if Z in (14, 16, 6):
                f[:] = 0; f[2 if Z != 6 else 1] = 1.0
        else:
            f = A.cie_fractions(Z, [T0])[0]
        x[Z] = f
    Ltot = {}
    cols = {}
    ab = A.ABUND_A09
    prev = None
    for k in range(nstep - 1):
        Ta, Tb = T[k], T[k + 1]
        Tm = np.sqrt(Ta * Tb)
        dT = Tb - Ta
        dtdT = (mu * M_H) ** 2 * KAPPA0 * P * Tm ** 0.5 / (2.5 * J ** 2 * K_B ** 2)
        dt = dtdT * dT
        fr = {Z: x[Z][:, None] for Z in ELEM}
        nenh = float(ne_nH_of(fr, ab)[0])
        ntnh = float(ntot_nH_of(fr, ab)[0])
        nH = P / (K_B * Tm * ntnh)
        ne = nenh * nH
        Ms = rate_matrices(Tm)
        for Z in ELEM:
            x[Z] = np.clip(expm(Ms[Z] * ne * dt) @ x[Z], 0, None)
            x[Z] /= x[Z].sum()
        fr = {Z: x[Z][:, None] for Z in ELEM}
        Lm = LN.synth(np.array([Tm]), fr, ne_nH=nenh, ne=max(ne, 1e-3), ab=ab)
        v = J * K_B * Tm / (mu * M_H * P)
        dx = v * dt
        for key, val in Lm.items():
            Ltot[key] = Ltot.get(key, 0.0) + float(val[0]) * nH ** 2 * dx
        for Z, st, nm in ((6, 3, "CIV"), (8, 5, "OVI"), (7, 4, "NV"), (8, 2, "OIII"), (14, 3, "SiIV")):
            cols[nm] = cols.get(nm, 0.0) + ab[Z] * x[Z][st] * nH * dx
    E5007 = A.H_PL * A.C_L / 5.0082e-5
    I5007_R = Ltot["OIII5007"] / E5007 / 1e6          # face-on Rayleigh (photon rate / 4pi / (1e6/4pi))
    return Ltot, cols, float(I5007_R)


def static_front(P_k, Tcie, frcie, LamCIE, T_c=1e4, T_h=1.6e6, ab=None):
    ab = ab or AB_GS
    m = (Tcie >= T_c) & (Tcie <= T_h)
    T = Tcie[m]
    fr = {Z: v[:, m] for Z, v in frcie.items()}
    nenh = ne_nH_of(fr, ab)
    ntnh = ntot_nH_of(fr, ab)
    P = P_k * K_B
    nH = P / (K_B * T * ntnh)
    lam = np.exp(np.interp(np.log(T), np.log(LamCIE[0]), np.log(LamCIE[1])))
    src = KAPPA0 * T ** 2.5 * nenh * nH ** 2 * lam            # F dF/dT
    F2 = np.concatenate([[0.0], np.cumsum(0.5 * (src[1:] + src[:-1]) * np.diff(T)) * 2])
    F = np.sqrt(np.maximum(F2, 1e-300))
    F[0] = F[1] * 0.5
    dxdT = KAPPA0 * T ** 2.5 / F
    Lm = LN.synth(T, fr, ne_nH=nenh, ne=max(1e-3, float(np.median(nenh * nH))), ab=ab)
    Lint = {k: float(np.trapezoid(v * nH ** 2 * dxdT, T)) for k, v in Lm.items() if np.ndim(v)}
    Lint["total_radiated"] = float(np.trapezoid(nenh * nH ** 2 * lam * dxdT, T))
    cols = {"CIV": float(np.trapezoid(ab[6] * fr[6][3] * nH * dxdT, T)),
            "OVI": float(np.trapezoid(ab[8] * fr[8][5] * nH * dxdT, T)),
            "OIII": float(np.trapezoid(ab[8] * fr[8][2] * nH * dxdT, T))}
    thick = float(np.trapezoid(dxdT, T))
    E5007 = A.H_PL * A.C_L / 5.0082e-5
    return Lint, cols, float(Lint["OIII5007"] / E5007 / 1e6), float(F[-1]), thick


# ------------------------------------------------------------------------------------------------
def main():
    res = {"notes": __doc__, "I_obs_R": I_OBS}
    G, delta = sheet_G()
    res["geometry"] = {"G_sheet": float(G), "thickness_pc": float(delta),
                       "formula": "G = 2 sqrt(2 R_c/Delta), R_c = 60 pc, Delta = 0.28 deg x 280 pc"}
    Tc, frc = read_gs_ion("tab2.txt")
    lamCIE = read_gs_cool("tab13.txt", 4)
    lamCIE2 = read_gs_cool("tab13.txt", 5)
    tracks = {"isochoric_Z1": ("tab3.txt", "tab14.txt", 4, 1.5), "isobaric_Z1": ("tab4.txt", "tab15.txt", 4, 2.5),
              "isochoric_Z2": ("tab11.txt", "tab14.txt", 5, 1.5), "isobaric_Z2": ("tab12.txt", "tab15.txt", 5, 2.5)}
    res["cooling_flows"] = {}
    for name, (fi, fc, col, cp) in tracks.items():
        T, fr = read_gs_ion(fi)
        lam = read_gs_cool(fc, col)
        ab = dict(AB_GS)
        if name.endswith("Z2"):
            ab = {Z: (v * 2 if Z > 2 else v) for Z, v in AB_GS.items()}
        out = {}
        for Thi in (1.6e6, 3e5):
            for Tlo in (1e4, 2e4, 3e4, 5e4):
                for tail in ((False, True) if Tlo == 1e4 else (False,)):
                    E = integrate_track(T, fr, (lam[0], lam[1]), cp, ab, Thi, Tlo, tail=tail)
                    key = f"T{Thi:.1e}-{Tlo:.0e}" + ("+tail" if tail else "")
                    s = summarize(E)
                    s.update(required_flux(E))
                    E100 = integrate_track(T, fr, (lam[0], lam[1]), cp, ab, Thi, Tlo, tail=tail, ne_dens=100.0)
                    s["OII3729_3726_at_ne1"] = s["OII3729_3726"]
                    s["OII3729_3726_at_ne100"] = float(E100["OII3729"] / E100["OII3726"])
                    out[key] = s
        # T-resolved emission weighting (for the figure)
        m = T >= 1e4
        frm = {Z: v[:, m] for Z, v in fr.items()}
        nenh, ntnh = ne_nH_of(frm, ab), ntot_nH_of(frm, ab)
        Lm = LN.synth(T[m], frm, ne_nH=nenh, ab=ab)
        lamm = np.exp(np.interp(np.log(T[m]), np.log(lam[0]), np.log(lam[1])))
        w = cp * K_B * ntnh / (nenh * lamm) * T[m]
        out["_dE_dlnT"] = {"T": T[m].tolist(), "OIII5007": (Lm["OIII5007"] * w).tolist(),
                           "CIV1549": (Lm["CIV1549"] * w).tolist(), "OII3727": ((Lm["OII3726"] + Lm["OII3729"]) * w).tolist(),
                           "Hb": (Lm["Hb"] * w).tolist(), "total": (nenh * lamm * w).tolist(),
                           "xOpp": frm[8][2].tolist(), "xOp": frm[8][1].tolist(), "xC3": frm[6][3].tolist()}
        res["cooling_flows"][name] = out
        f = out["T1.6e+06-1e+04"]
        ft = out["T1.6e+06-1e+04+tail"]
        print(f"{name:14s} 1.6e6->1e4: OIII/Hb={f['OIII_Hb']:.1f} (tail {ft['OIII_Hb']:.1f}) OII/OIII={f['OII_OIII']:.3f}"
              f" 4363={f['O4363_5007']:.3f} CIV/OIII(ph)={f['CIV_OIII_phot']:.2f} FUV={f['FUV_CU_per_R']:.0f} CU/R"
              f" f5007={f['f5007']:.2e} Fneed={f['energy_flux']:.2e}")
    # ---- warm CIE layer
    fr_own = own_cie(Tc[(Tc >= 1e4) & (Tc <= 3e5)])
    Tsel = Tc[(Tc >= 1e4) & (Tc <= 3e5)]
    cie = {"T": Tsel.tolist(), "GS07": [], "own": []}
    for src, frx in (("GS07", {Z: v[:, (Tc >= 1e4) & (Tc <= 3e5)] for Z, v in frc.items()}), ("own", fr_own)):
        abx = AB_GS if src == "GS07" else A.ABUND_A09
        nenh = ne_nH_of(frx, abx)
        for ne_d in (1.0, 100.0):
            Lw = LN.synth(Tsel, frx, ne_nH=nenh, ne=ne_d, ab=abx)
            lam = np.exp(np.interp(np.log(Tsel), np.log(lamCIE[0]), np.log(lamCIE[1])))
            for k, T in enumerate(Tsel):
                Lk = {kk: (vv[k] if np.ndim(vv) else vv) for kk, vv in Lw.items()}
                s = summarize(Lk)
                s["T"] = float(T)
                s["ne_dens"] = ne_d
                s["f5007"] = float(Lk["OIII5007"] / (nenh[k] * lam[k]))
                s["xOpp"] = float(frx[8][2][k]); s["xOp"] = float(frx[8][1][k]); s["xC3"] = float(frx[6][3][k])
                s["xH0"] = float(frx[1][0][k])
                cie[src].append(s)
    res["warm_CIE"] = cie
    # consistency window (GS07 CIE, ne=1)
    win = [s for s in cie["GS07"] if s["ne_dens"] == 1.0 and 1.0 <= s["OII_OIII"] <= 2.0]
    res["warm_CIE_window_OII_OIII_1_2"] = win
    for s in win:
        print(f"CIE T={s['T']:.0f} OII/OIII={s['OII_OIII']:.2f} OIII/Hb={s['OIII_Hb']:.1f} 4363={s['O4363_5007']:.3f}"
              f" CIV/OIII(ph)={s['CIV_OIII_phot']:.2f} FUV={s['FUV_CU_per_R']:.0f} NUV={s['NUV_CU_per_R']:.0f} NII/Ha={s['NII_Ha']:.2f}"
              f" SII/Ha={s['SII_Ha']:.2f} Hbcoll={s['Hb_coll_frac']:.2f} f5007={s['f5007']:.3f} xH0={s['xH0']:.1e}")
    # own-vs-GS07 CIE check for O++ / C3+
    chk = []
    for T in (2e4, 3e4, 4e4, 5e4, 7e4, 1e5, 2e5):
        i = int(np.argmin(abs(Tsel - T)))
        chk.append(dict(T=float(Tsel[i]), GS07_Opp=float(frc[8][2][(Tc >= 1e4) & (Tc <= 3e5)][i]), own_Opp=float(fr_own[8][2][i]),
                        GS07_Op=float(frc[8][1][(Tc >= 1e4) & (Tc <= 3e5)][i]), own_Op=float(fr_own[8][1][i]),
                        GS07_C3=float(frc[6][3][(Tc >= 1e4) & (Tc <= 3e5)][i]), own_C3=float(fr_own[6][3][i])))
    res["CIE_check_own_vs_GS07"] = chk
    for c in chk:
        print("CIE check", {k: round(v, 4) if v < 1e3 else round(v) for k, v in c.items()})
    # ---- static radiative conduction fronts
    res["static_fronts"] = {}
    for P_k in (1.5e4, 2e6):
        L, cols, I5007, Fh, thick = static_front(P_k, Tc, frc, lamCIE)
        s = summarize(L)
        need = 2.0 / G
        s.update(dict(P_k=P_k, I5007_faceon_R=I5007, fronts_needed_for_2R_at_G=float(need / I5007),
                      heat_flux_hot_side=Fh, total_radiated=L["total_radiated"], thickness_pc=thick / 3.0857e18,
                      cols=cols, f5007=float(L["OIII5007"] / L["total_radiated"])))
        res["static_fronts"][f"P{P_k:.1e}"] = s
        print(f"static front P/k={P_k:.1e}: I5007 face-on={I5007:.2e} R; need {need / I5007:.0f} fronts (G={G:.1f});"
              f" CIV/OIII(ph)={s['CIV_OIII_phot']:.1f} 4363={s['O4363_5007']:.3f} OII/OIII={s['OII_OIII']:.2f}"
              f" OIII/Hb={s['OIII_Hb']:.1f} FUV={s['FUV_CU_per_R']:.0f} N(OVI)={cols['OVI']:.2e} N(CIV)={cols['CIV']:.2e} thick={thick/3.0857e18:.2e} pc")
    # ---- evaporation fronts
    Th, neh = 1.6e6, 5e-3
    P_k = 1.92 * neh * Th
    ckT = K_B * Th
    q_sat = 0.4 * np.sqrt(2 * ckT / (np.pi * 9.109e-28)) * neh * ckT
    J_sat = q_sat / (2.5 * ckT / (0.6 * M_H))
    res["evaporation"] = {"P_k": P_k, "q_sat": float(q_sat), "J_sat": float(J_sat), "runs": {}}
    for Rpc in (1.0, 3.0, 10.0, 60.0):
        Jcl = 4 * 0.6 * M_H * KAPPA0 * Th ** 2.5 / (25 * K_B * Rpc * 3.0857e18)
        J = min(Jcl, J_sat)
        L, cols, I5007 = evaporation_front(J, P_k)
        s = summarize(L)
        s.update(dict(R_cloud_pc=Rpc, J_classical=float(Jcl), J_used=float(J), I5007_faceon_R=I5007,
                      fronts_needed_for_2R_at_G=float(2.0 / G / max(I5007, 1e-30)), cols=cols))
        res["evaporation"]["runs"][f"R{Rpc:g}pc"] = s
        print(f"evaporation R={Rpc}pc J={J:.2e} (cl {Jcl:.2e}, sat {J_sat:.2e}): I5007={I5007:.2e} R face-on; "
              f"4363={s['O4363_5007']:.3f} OII/OIII={s['OII_OIII']:.3f} OIII/Hb={s['OIII_Hb']:.1f} CIV/OIII(ph)={s['CIV_OIII_phot']:.2f}"
              f" N(OVI)={cols['OVI']:.2e} N(CIV)={cols['CIV']:.2e}")
    with open(OUTJ, "w") as f:
        json.dump(res, f, indent=1)
    plot(res)


def plot(res):
    from silk_style import new_sky_figure, label, finish, OIII_C, WARN, OK, FGDIM, FG, ACCENT, DANGER, SII_C, HA_C, NII_C
    import matplotlib.pyplot as plt
    fig, axs = new_sky_figure(13.5, 5.2, "Collisional alternatives: cooling flows, warm layers, conduction fronts",
                              "line-emission weighting of GS07 isobaric cooling from 1.6 MK, and warm-CIE line ratios vs T",
                              1, 3)
    ax = axs[0]
    d = res["cooling_flows"]["isobaric_Z1"]["_dE_dlnT"]
    T = np.array(d["T"])
    tot = np.trapezoid(np.array(d["total"]), np.log(T))
    for k, c, lab in (("total", FGDIM, "all cooling"), ("CIV1549", WARN, "C IV 1549"), ("OIII5007", OIII_C, "[O III] 5007"),
                      ("OII3727", NII_C, "[O II] 3727"), ("Hb", OK, "Hβ")):
        y = np.array(d[k]) / tot
        ax.plot(T, np.maximum(y, 1e-8), color=c, lw=1.6, label=lab)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-6, 2); ax.set_xlim(1e4, 2e6)
    ax.axvspan(3.1e4, 5.3e4, color=DANGER, alpha=0.12, lw=0)
    label(ax, 3.2e4, 0.6, "T_e([O III])", color=DANGER, size=7.5)
    ax.set_xlabel("T (K)"); ax.set_ylabel("dE/dlnT / E_total  (isobaric, Z = 1)")
    ax.legend(loc="lower right", fontsize=7.5); ax.grid(alpha=0.3, which="both")
    ax = axs[1]
    cie = [s for s in res["warm_CIE"]["GS07"] if s["ne_dens"] == 1.0]
    Tw = np.array([s["T"] for s in cie])
    for k, c, lab in (("OIII_Hb", OIII_C, "[O III]/Hβ"), ("OII_OIII", NII_C, "[O II]/[O III]"),
                      ("O4363_5007", DANGER, "4363/5007 ×100"), ("CIV_OIII_phot", WARN, "C IV/[O III] (photons)")):
        y = np.array([s[k] for s in cie]) * (100 if k == "O4363_5007" else 1)
        ax.plot(Tw, y, color=c, lw=1.6, label=lab)
    ax.axhspan(15, 60, color=OIII_C, alpha=0.08, lw=0)
    ax.axhspan(1.1, 1.9, color=NII_C, alpha=0.18, lw=0)
    ax.axhspan(5.9, 9.1, color=DANGER, alpha=0.12, lw=0)
    ax.axhline(0.45, color=WARN, ls="--", lw=0.9)
    label(ax, 1.05e4, 0.5, "C IV/[O III] FUV ceiling", color=WARN, size=7)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1e4, 3e5); ax.set_ylim(1e-3, 1e4)
    ax.set_xlabel("T (K) of a warm CIE layer"); ax.set_ylabel("line ratio (energy unless noted)")
    ax.legend(loc="upper left", fontsize=7.5); ax.grid(alpha=0.3, which="both")
    ax = axs[2]
    cfc = {"isochoric_Z1": "#2979ff", "isobaric_Z1": OIII_C, "isobaric_Z2": NII_C}
    for name, col in cfc.items():
        for key, mk in (("T1.6e+06-1e+04", "o"), ("T1.6e+06-2e+04", "^"), ("T1.6e+06-3e+04", "v"), ("T3.0e+05-1e+04", "D")):
            s = res["cooling_flows"][name][key]
            ax.scatter(s["CIV_OIII_phot"], s["O4363_5007"], color=col, marker=mk, s=30, zorder=3,
                       label=name.replace("_", " ") if key == "T1.6e+06-1e+04" else None)
    for s in cie:
        if 4e4 < s["T"] < 1.5e5:
            ax.scatter(s["CIV_OIII_phot"], s["O4363_5007"], color=FGDIM, marker=".", s=22, zorder=2)
    for k, s in res["static_fronts"].items():
        ax.scatter(s["CIV_OIII_phot"], s["O4363_5007"], color=WARN, marker="*", s=110, zorder=4,
                   label="static conduction front" if "1.5" in k else None)
    for i, (k, s) in enumerate(res["evaporation"]["runs"].items()):
        ax.scatter(s["CIV_OIII_phot"], s["O4363_5007"], color=DANGER, marker="X", s=50, zorder=4,
                   label="evaporation front" if i == 0 else None)
    ax.scatter([], [], color=FGDIM, marker=".", s=22, label="warm CIE layer 4e4–1.5e5 K")
    ax.add_patch(plt.Rectangle((1e-3, 0.059), 0.45 - 1e-3, 0.096 - 0.059, color=OK, alpha=0.25, lw=0))
    label(ax, 1.3e-3, 0.064, "Silk: 4363/5007 and\nall-FUV-as-C IV ceiling", color=OK, size=7.5)
    label(ax, 1.3e-3, 0.032, "cooling flows: ○ 1.6 MK→1e4 K, △ →2e4, ▽ →3e4, ◇ 3e5→1e4 K", color=FG, size=6.5)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1e-3, 50); ax.set_ylim(0.03, 0.2)
    ax.set_xlabel("C IV 1549 / [O III] 5007  (photons)"); ax.set_ylabel("[O III] 4363/5007")
    ax.legend(loc="upper left", fontsize=6.5); ax.grid(alpha=0.3, which="both")
    finish(fig, FIG)


if __name__ == "__main__":
    main()
