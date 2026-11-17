"""verify_altmodel_budget: densities, photon budgets, photoionisation(+heating) hybrids, shock-precursor
and fossil (recombining) scenarios for the Silk [O III] filament.

1. Density of the O+ gas from [O II]3729/3726 (CHIANTI v7.1.4 O II, own 35-level model), at T(O II).
2. rms density needed for the observed I(5007) on a path L = G x Delta (sheet geometry), for each
   candidate (T, O++ fraction); the implied volume filling factor (n_rms / n_[O II])^2 and the thickness of
   the emitting layer if n = n_[O II].
3. Photoionisation: the ionisation parameter required for [O II]/[O III] = 1.1-1.9 (from my one-zone PIE
   grid, verify_altmodel_photoion.json) against the photon flux of each candidate source:
   WD J0736+1444 (Q(H) = 6.6e44 s^-1, Q(>35 eV) = 1.1e44 s^-1, context.md), Monogem hot-gas EUV, generic
   field; plus the Stromgren (recombination) budget.  Also [O III]5007 photons per H recombination.
4. Photoionisation + extra heating (fixed-T PIE, T = 1.2e4-6e4 K) for 81 kK BB and 0.14 keV spectra:
   the U at which [O II]/[O III] = 1.37, the other ratios, the GALEX FUV/NUV per R and the heating needed.
5. Fossil: isochoric, unheated recombination of gas initially ionised (T0 = 2e4 or 4e4 K) at n_H = 100 and
   0.5 cm^-3, with my NEI solver (CI + RR/DR + O/H charge exchange) and the photoion cooling table.
Outputs: data/derived/verify_altmodel_budget.json, figures/verify_altmodel_budget.png
Run: ~/.venvs/astro313/bin/python -I scripts/verify_altmodel_budget.py
"""
import json
import os
import sys

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, os.path.dirname(__file__))
import verify_altmodel_atomic as A  # noqa: E402
import verify_altmodel_lines as LN  # noqa: E402
import verify_altmodel_photoion as P  # noqa: E402
from verify_altmodel_common import sheet_G  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUTJ = os.path.join(ROOT, "data/derived/verify_altmodel_budget.json")
FIG = os.path.join(ROOT, "figures/verify_altmodel_budget.png")
PC, C_L, K_B = 3.0857e18, A.C_L, A.K_B
E5007 = A.H_PL * C_L / 5.0082e-5
I_OBS = dict(fid=2.0, lo=1.4, hi=4.0)
OII_DOUBLET = dict(fid=1.30, lo=1.24, hi=1.33, subset_lo=1.17, subset_hi=1.37, sdss=1.29, sdss_lo=1.07, sdss_hi=1.59)


def doublet_density(T):
    nes = np.logspace(-1, 3.5, 200)
    r = np.array([A.emis("OII3729", T, n) / A.emis("OII3726", T, n) for n in nes])
    out = {"low_density_limit": float(r[0])}
    for k, v in OII_DOUBLET.items():
        if v >= r[0]:
            out[k] = 0.0
        else:
            out[k] = float(np.exp(np.interp(-v, -r, np.log(nes))))
    return out


def n_rms_needed(T, xOpp, I_R, G, delta_pc, ne_nH=1.2, OH=4.9e-4):
    L = G * delta_pc * PC
    q_ph = float(LN.e_line("OIII5007", T, 1.0)[0]) / E5007        # photons cm^3 s^-1 per O++ per n_e
    # I_R = 1e-6 * L * n_e * n(O++) * q_ph ;  n(O++) = xOpp OH n_H ;  n_e = ne_nH n_H
    nH = np.sqrt(I_R / (1e-6 * L * ne_nH * xOpp * OH * q_ph))
    return float(nH), float(L / PC), q_ph


def main():
    res = {"notes": __doc__}
    G, delta = sheet_G()
    res["geometry"] = dict(G=float(G), delta_pc=float(delta))
    # ---------------- 1. densities
    dens = {f"T={T:.0f}": doublet_density(T) for T in (1.0e4, 1.25e4, 2.0e4)}
    res["OII_doublet_density"] = dens
    for k, v in dens.items():
        print("[O II] doublet n_e", k, {kk: round(vv, 1) for kk, vv in v.items()})
    ne_dbl = dens["T=12500"]["fid"]
    P_k = 2.2 * ne_dbl * 1.25e4
    res["P_k_from_doublet"] = dict(fid=P_k, lo=2.2 * dens["T=12500"]["subset_hi"] * 1.25e4,
                                   hi=2.2 * dens["T=12500"]["subset_lo"] * 1.25e4)
    print("P/k from doublet:", res["P_k_from_doublet"])
    # ---------------- 2. rms densities
    cases = {"PIE_1.2e4K_xOpp0.45": (1.2e4, 0.45), "warm_6e4K_xOpp0.45": (6e4, 0.45), "shock_like_4e4K_xOpp0.6": (4e4, 0.6)}
    rms = {}
    for name, (T, x) in cases.items():
        for I in (I_OBS["lo"], I_OBS["fid"], I_OBS["hi"]):
            nH, Lpc, q = n_rms_needed(T, x, I, G, delta)
            ne_true = ne_dbl if T < 2e4 else ne_dbl * 1.25e4 / T
            f = (1.2 * nH / ne_true) ** 2
            L_em = I / (1e-6 * 1.2 * (ne_true / 1.2) ** 2 * x * 4.9e-4 * q) / PC
            rms[f"{name}|I={I}"] = dict(nH_rms=nH, ne_rms=1.2 * nH, L_pc=Lpc, filling_factor_vs_doublet=float(f),
                                       emitting_path_pc_at_doublet_density=float(L_em))
    res["rms_density"] = rms
    for k, v in rms.items():
        print("rms", k, {kk: f"{vv:.3g}" for kk, vv in v.items()})
    # ---------------- 3. photoionisation budgets
    pj = json.load(open(os.path.join(ROOT, "data/derived/verify_altmodel_photoion.json")))
    needU = {}
    for name, rows in pj["models"].items():
        lu = np.array([r["logU"] for r in rows])
        oo = np.array([r["OII_OIII"] for r in rows])
        o3hb = np.array([r["OIII_Hb"] for r in rows])
        Ts = np.array([r["T"] for r in rows])
        # OII/OIII decreases with U
        out = {}
        for tag, v in (("1.9", 1.9), ("1.37", 1.37), ("1.1", 1.1)):
            if oo.min() < v < oo.max():
                l = float(np.interp(-np.log(v), -np.log(oo), lu))
                out[tag] = dict(logU=l, OIII_Hb=float(np.interp(l, lu, o3hb)), T=float(np.interp(l, lu, Ts)))
        out["max_OIII_Hb"] = float(o3hb.max())
        needU[name] = out
    res["PIE_needed_U"] = needU
    for k, v in needU.items():
        print("PIE", k, v)
    # source fluxes
    Q_wd, Q35_wd = 6.6e44, 1.1e44
    src = {}
    for dname, dpc in (("coplanar_6.4pc", 6.4), ("Silk_at_282pc_WD_at_150pc", 132.0)):
        Phi = Q_wd / (4 * np.pi * (dpc * PC) ** 2)
        src[f"WDJ0736_{dname}"] = dict(Phi_H=float(Phi), Phi_35eV=float(Q35_wd / (4 * np.pi * (dpc * PC) ** 2)))
    # Monogem hot gas EUV at its boundary: F = n_e n_H Lambda R / 3 (uniform optically thin sphere)
    Tg, lamg = [], []
    with open(os.path.join(ROOT, "data/raw/verify_altmodel/gs07/tab13.txt")) as f:
        for l in f:
            p = l.split()
            if len(p) == 6:
                try:
                    Tg.append(float(p[0])); lamg.append(float(p[4]))
                except ValueError:
                    pass
    lam16 = float(np.exp(np.interp(np.log(1.62e6), np.log(Tg), np.log(lamg))))
    for ne in (3.1e-3, 5e-3, 7.2e-3):
        F = ne * ne / 1.2 * lam16 * 60 * PC / 3
        for Emean in (30.0, 100.0):
            src[f"Monogem_ne{ne:.1e}_Emean{Emean:.0f}eV"] = dict(F_erg=float(F), Phi_H=float(F / (Emean * A.EV)))
    src["generic_upper_1e6"] = dict(Phi_H=1e6)
    res["Lambda_CIE_1.6MK_GS07"] = lam16
    # U available at n_H = n_rms (PIE case, fid I) and at the doublet density
    nH_rms = rms["PIE_1.2e4K_xOpp0.45|I=2.0"]["nH_rms"]
    for k, v in src.items():
        v["logU_at_nrms"] = float(np.log10(v["Phi_H"] / (nH_rms * C_L)))
        v["logU_at_doublet"] = float(np.log10(v["Phi_H"] / (ne_dbl / 1.2 * C_L)))
        need = 10 ** needU["BB81kK_WDJ0736"]["1.37"]["logU"] if "WD" in k else 10 ** needU["hotgas_kT140eV"]["1.37"]["logU"]
        v["shortfall_at_nrms"] = float(need * nH_rms * C_L / v["Phi_H"])
        v["shortfall_at_doublet"] = float(need * ne_dbl / 1.2 * C_L / v["Phi_H"])
    res["sources"] = src
    for k, v in src.items():
        print("source", k, {kk: f"{vv:.3g}" for kk, vv in v.items()})
    # 5007 photons per H recombination and the face-on photon requirement
    I = I_OBS["fid"]
    ph_face = 4 * np.pi * I * 1e6 / (4 * np.pi) / G
    per_rec = {}
    for o3hb in (15.0, 47.0):
        T = 1.2e4
        hb_per_rec = float(A.j_hbeta(T) / (A.alpha_B(T) * A.H_PL * C_L / 4.8627e-5))
        n5007 = o3hb * (4861.3 / 5006.8) * hb_per_rec
        per_rec[f"OIII_Hb={o3hb}"] = dict(Hb_per_rec=hb_per_rec, n5007_per_rec=n5007,
                                         ionising_flux_needed_faceon=float(ph_face / n5007))
    res["recombination_budget"] = dict(face_on_5007_photon_flux=float(ph_face), cases=per_rec)
    print("recombination budget", res["recombination_budget"])
    # ---------------- 4. PIE + extra heating (fixed T) incl. Si
    P.ELEM = [1, 2, 6, 7, 8, 10, 14, 16]
    P.AB = {z: A.ABUND_A09[z] for z in P.ELEM}
    P.COOL = P.build_cooling_table()
    hyb = {}
    for sname, (kind, par) in {"BB81kK": ("bb", 8.13e4), "hotgas140eV": ("thermal", 140.0), "precursor28eV": ("thermal", 28.0)}.items():
        f = P.spectrum(kind, par)
        smom, hmom = P.sigma_moments(f)
        rows = []
        for Tfix in (None, 1.2e4, 1.5e4, 2e4, 3e4, 4e4, 6e4):
            lus = np.arange(-6.0, -1.49, 0.125)
            grid = []
            for lu in lus:
                T, frac, ne = P.solve(10 ** lu, smom, hmom, Tfix=Tfix)
                Ls = LN.synth(np.array([T]), {Z: frac[Z][:, None] for Z in P.ELEM}, ne_nH=ne, ne=100.0, ab=P.AB)
                Ls = {k: (v[0] if np.ndim(v) else v) for k, v in Ls.items()}
                r = LN.ratios(Ls)
                grid.append((lu, T, frac, ne, Ls, r))
            oo = np.array([g[5]["OII_OIII"] for g in grid])
            ok = np.isfinite(oo) & (oo > 0)
            if not ((oo[ok] > 1.37).any() and (oo[ok] < 1.37).any()):
                continue
            k = int(np.argmin(abs(np.log(oo[ok] / 1.37))))
            lu, T, frac, ne, Ls, r = [g for g, o in zip(grid, ok) if o][k]
            uv = LN.uv_CU_per_R5007(Ls)
            h, c = P.cool_heat(T, frac, ne, 10 ** lu, smom, hmom, P.rates_at(T)[0])
            row = dict(Tfix=Tfix, T=float(T), logU=float(lu), heat_extra_over_photo=float(max(c - h, 0) / h),
                       FUV=uv["FUV"]["total"], FUV_parts=uv["FUV"]["parts"], NUV=uv["NUV"]["total"],
                       xOpp=float(frac[8][2]), xHp=float(frac[1][1]), **{kk: float(vv) for kk, vv in r.items()})
            rows.append(row)
            print(f"hybrid {sname:13s} Tfix={Tfix} T={T:.0f} logU={lu:.2f} OIII/Hb={r['OIII_Hb']:.1f} OII/OIII={r['OII_OIII']:.2f}"
                  f" 4363={r['O4363_5007']:.3f} 7325/3727={r['OII_7325_3727']:.3f} NeIII={r['NeIII_OIII']:.2f} HeII={r['HeII_OIII']:.3f}"
                  f" NII/Ha={r['NII_Ha']:.2f} SII/Ha={r['SII_Ha']:.2f} CIV/OIII={r['CIV_OIII_phot']:.3f} FUV={uv['FUV']['total']:.0f}"
                  f" NUV={uv['NUV']['total']:.0f} heat x{row['heat_extra_over_photo']:.1f}")
        hyb[sname] = rows
    res["PIE_plus_heating"] = hyb
    # ---------------- 5. fossil recombination
    res["fossil"] = fossil_runs()
    with open(OUTJ, "w") as fo:
        json.dump(res, fo, indent=1, default=float)
    plot(res)


# ---------------------------------------------------------------------------------------------------
ELEM = [1, 2, 6, 7, 8, 10, 16]


def rate_matrices(T, nH0_ne, nHp_ne):
    out = {}
    for Z in ELEM:
        n = Z + 1
        M = np.zeros((n, n))
        for i in range(Z):
            c = float(A.ci_rate(Z, i + 1, T))
            a = float(A.rec_rate(Z, i + 2, T)[0])
            if Z == 8 and i == 0:
                c += (8.0 / 9.0) * 1e-9 * nHp_ne
                a += 1e-9 * nH0_ne
            if Z == 8 and i == 1:
                a += 1e-9 * nH0_ne
            M[i + 1, i] += c; M[i, i] -= c
            M[i, i + 1] += a; M[i + 1, i + 1] -= a
        out[Z] = M
    return out


def fossil_runs():
    P.ELEM = [1, 2, 6, 7, 8, 10, 16]
    P.AB = {z: A.ABUND_A09[z] for z in P.ELEM}
    ab = P.AB
    runs = {}
    for T0 in (2e4, 4e4):
        for nH in (100.0, 0.5):
            # initial state: H+, He+, O++ 0.7 / O+ 0.3, C++, N++, Ne++, S++ (photoionised / shock-heated gas)
            x = {}
            for Z in ELEM:
                f = np.zeros(Z + 1)
                if Z == 1:
                    f[1] = 1.0
                elif Z == 2:
                    f[1] = 1.0
                elif Z == 8:
                    f[2], f[1] = 0.7, 0.3
                else:
                    f[2] = 1.0
                x[Z] = f
            T = T0
            t = 0.0
            rows = []
            nsteps = 4000
            tmax = 3e13 / nH * 30
            dts = np.diff(np.concatenate([[0], np.logspace(np.log10(tmax) - 7, np.log10(tmax), nsteps)]))
            for dt in dts:
                nenh = sum(ab[Z] * np.dot(np.arange(Z + 1), x[Z]) for Z in ELEM)
                ntnh = sum(ab[Z] for Z in ELEM) + nenh
                ne = nenh * nH
                Ms = rate_matrices(T, x[1][0] * nH / ne, x[1][1] * nH / ne)
                for Z in ELEM:
                    x[Z] = np.clip(expm(Ms[Z] * ne * dt) @ x[Z], 0, None)
                    x[Z] /= x[Z].sum()
                frac = {Z: x[Z] for Z in ELEM}
                h, c = P.cool_heat(max(T, 3.1e3), frac, nenh, 0.0, {Z: np.zeros(Z) for Z in ELEM},
                                   {Z: np.zeros(Z) for Z in ELEM}, P.rates_at(max(T, 3.1e3))[0])
                dTdt = -(c * nH * nH) / (1.5 * ntnh * nH * K_B)
                T = max(T + dTdt * dt, 3.0e3)
                t += dt
                if len(rows) == 0 or t > rows[-1]["t_yr"] * 3.156e7 * 1.08:
                    Ls = LN.synth(np.array([T]), {Z: x[Z][:, None] for Z in ELEM}, ne_nH=nenh, ne=max(ne, 0.1), ab=ab)
                    Ls = {k: (v[0] if np.ndim(v) else v) for k, v in Ls.items()}
                    r = LN.ratios(Ls)
                    rows.append(dict(t_yr=t / 3.156e7, T=float(T), xOpp=float(x[8][2]), xOp=float(x[8][1]), xHp=float(x[1][1]),
                                     e5007_per_nH2=float(Ls["OIII5007"]), OIII_Hb=float(r["OIII_Hb"]), OII_OIII=float(r["OII_OIII"]),
                                     O4363_5007=float(r["O4363_5007"]), OII_7325_3727=float(r["OII_7325_3727"])))
            e0 = rows[0]["e5007_per_nH2"]
            tfade = next((r["t_yr"] for r in rows if r["e5007_per_nH2"] < e0 / np.e), None)
            ok = [r for r in rows if r["OIII_Hb"] >= 15 and 1.0 <= r["OII_OIII"] <= 2.0 and r["O4363_5007"] >= 0.05]
            runs[f"T0={T0:.0e}|nH={nH}"] = dict(rows=rows, t_fade_yr=tfade, n_epochs_matching=len(ok),
                                                 matching_window_yr=[ok[0]["t_yr"], ok[-1]["t_yr"]] if ok else None)
            print(f"fossil T0={T0:.0e} nH={nH}: t_fade(5007,1/e)={tfade:.3g} yr; epochs matching [O III]/Hb>=15 & 1<[O II]/[O III]<2 &"
                  f" 4363>=0.05: {len(ok)} {runs[f'T0={T0:.0e}|nH={nH}']['matching_window_yr']}")
            for r in rows[:: max(1, len(rows) // 8)]:
                print("    ", {k: (f"{v:.3g}") for k, v in r.items()})
    return runs


def plot(res):
    from silk_style import new_sky_figure, label, finish, OIII_C, WARN, OK, FGDIM, FG, ACCENT, DANGER, SII_C, NII_C
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    fig, axs = new_sky_figure(14.5, 5.6, "Photoionisation, heating and fossil tests",
                              "ionisation parameter available vs needed; photoionised + heated gas; unheated recombination", 1, 3,
                              gridspec_kw=dict(width_ratios=[1.15, 1, 1]))
    ax = axs[0]
    src = res["sources"]
    names = [k for k in src]
    y = np.arange(len(names))
    for i, k in enumerate(names):
        ax.plot([src[k]["logU_at_doublet"], src[k]["logU_at_nrms"]], [i, i], color=FGDIM, lw=1)
        ax.scatter(src[k]["logU_at_nrms"], i, color=OIII_C, s=25, zorder=3)
        ax.scatter(src[k]["logU_at_doublet"], i, color=DANGER, s=25, zorder=3)
    lu = res["PIE_needed_U"]["BB81kK_WDJ0736"]
    ax.axvspan(lu["1.9"]["logU"], lu["1.1"]["logU"], color=OK, alpha=0.35, lw=0)
    label(ax, lu["1.9"]["logU"] - 0.15, 4.0, "needed for\n[O II]/[O III]\n≈ 1.1–1.9", color=OK, size=7.5, ha="right")
    ax.set_yticks(y)
    ax.set_yticklabels([k.replace("_", " ").replace("Monogem ", "Monogem EUV ")[:34] for k in names], fontsize=6.5)
    ax.set_xlabel("log U available\n(blue: n = n_rms, red: n from [O II] 3729/3726)")
    ax.set_xlim(-12, -2); ax.grid(alpha=0.3)
    ax = axs[1]
    cols = {"BB81kK": OIII_C, "hotgas140eV": DANGER, "precursor28eV": WARN}
    for s, rows in res["PIE_plus_heating"].items():
        T = [r["T"] for r in rows]
        ax.plot(T, [r["OIII_Hb"] for r in rows], "-o", ms=3, color=cols[s], lw=1.4)
        ax.plot(T, [r["FUV"] for r in rows], "--s", ms=3, color=cols[s], lw=1.0)
        ax.plot(T, [100 * r["OII_7325_3727"] for r in rows], ":^", ms=3, color=cols[s], lw=1.0)
    ax.axhspan(15, 47, color=OIII_C, alpha=0.10, lw=0)
    ax.axhspan(68, 92, color=OK, alpha=0.22, lw=0)
    ax.axhspan(2.2, 3.3, color=SII_C, alpha=0.25, lw=0)
    ax.axhspan(59, 96, color=DANGER, alpha=0.0, lw=0)
    label(ax, 1.25e4, 100, "FUV observed (CU/R)", color=OK, size=7)
    label(ax, 1.25e4, 16.5, "[O III]/Hβ ≥ 15", color=OIII_C, size=7)
    label(ax, 2.0e4, 2.45, "7325/3727 observed (×100)", color=SII_C, size=7)
    h = [Line2D([], [], color=c, lw=2) for c in cols.values()]
    h += [Line2D([], [], color=FG, ls=ls, marker=m, ms=3) for ls, m in (("-", "o"), ("--", "s"), (":", "^"))]
    ax.legend(h, ["81 kK BB", "0.14 keV gas", "28 eV precursor", "[O III]/Hβ", "FUV CU per R", "100×7325/3727"],
              fontsize=6.5, loc="lower right", ncol=2)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(0.12, 400)
    ax.set_xlabel("T (K) of photoionised gas held at T\n(U set by [O II]/[O III] = 1.37; 4363/5007 = 0.083 needs ~4.6e4 K)")
    ax.grid(alpha=0.3, which="both")
    ax = axs[2]
    for k, run in res["fossil"].items():
        rows = run["rows"]
        t = np.array([r["t_yr"] for r in rows])
        e = np.array([r["e5007_per_nH2"] for r in rows])
        c = ACCENT if "nH=100" in k else NII_C
        ls = "-" if "2e+04" in k else "--"
        ax.plot(t, e / e[0], color=c, ls=ls, lw=1.5, label=k.replace("|", ", ").replace("e+04", "e4"))
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-4, 3)
    ax.axvline(1.1e5, color=WARN, lw=0.9, ls=":")
    label(ax, 1.2e5, 1.5, "Monogem age", color=WARN, size=7.5)
    ax.set_xlabel("time since ionisation/heating stopped (yr)"); ax.set_ylabel("[O III]5007 emissivity / initial (fixed n_H)")
    ax.legend(fontsize=7, loc="lower left"); ax.grid(alpha=0.3, which="both")
    finish(fig, FIG)


if __name__ == "__main__":
    main()
