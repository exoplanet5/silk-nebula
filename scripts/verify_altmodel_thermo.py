"""verify_altmodel_thermo: mechanism-independent 'thermometer' for the Silk [O III] gas.

For gas at a single electron temperature T in which H is fully ionised and O is in O+ / O++ only
(the most favourable case for high [O III]/Hbeta), the line ratios are
    [O III]5007/Hb = (O/H) f++ eps5007(T) / jHb(T)
    [O II]3727/[O III]5007 = (f+/f++) eps3727(T)/eps5007(T)
so for an observed r = [O II]/[O III], f++ = 1/(1 + r eps5007/eps3727) and the *maximum* possible
[O III]/Hbeta at that T follows.  Requiring it to reach the observed value gives T_min.  Any mixture
of zones has a ratio <= the largest single-zone ratio, so T_min bounds the [O III]-emitting gas.
Also: T from [O III]4363/5007, [O II]7325/3727 and [N II]5755/6583 with the same atomic data.
Outputs: data/derived/verify_altmodel_thermo.json, figures/verify_altmodel_thermometer.png
"""
import json
import os
import sys

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(__file__))
import verify_altmodel_atomic as A  # noqa: E402
from verify_altmodel_common import OBS  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "data/derived/verify_altmodel_thermo.json")
FIG = os.path.join(ROOT, "figures/verify_altmodel_thermometer.png")

NE = 10.0
TG = np.logspace(np.log10(5e3), np.log10(2e5), 120)


def tab(ne=NE):
    e5007 = np.array([A.emis("OIII5007", t, ne) for t in TG])
    e4959 = np.array([A.emis("OIII4959", t, ne) for t in TG])
    e4363 = np.array([A.emis("OIII4363", t, ne) for t in TG])
    e3727 = np.array([A.oii_3727_emis(t, ne) for t in TG])
    e7325 = np.array([A.oii_7325_emis(t, ne) for t in TG])
    n6583 = np.array([A.emis("NII6583", t, ne) for t in TG])
    n5755 = np.array([A.emis("NII5755", t, ne) for t in TG])
    hb = A.j_hbeta(TG)
    return dict(e5007=e5007, e4959=e4959, e4363=e4363, e3727=e3727, e7325=e7325, n6583=n6583, n5755=n5755, hb=hb)


def invert(ratio_curve, value):
    lr = np.log(ratio_curve)
    lt = np.log(TG)
    # ratio is monotonic increasing with T for auroral/nebular ratios
    if value <= ratio_curve.min() or value >= ratio_curve.max():
        return np.nan
    return float(np.exp(np.interp(np.log(value), lr, lt)))


def max_oiii_hb(t, r, OH, hbscale=1.0):
    fpp = 1.0 / (1.0 + r * t["e5007"] / t["e3727"])
    return OH * fpp * t["e5007"] / (t["hb"] * hbscale), fpp


def tmin(curve, target):
    if curve.max() < target:
        return np.nan
    i = np.argmax(curve >= target)
    if i == 0:
        return float(TG[0])
    return float(np.exp(np.interp(target, [curve[i - 1], curve[i]], [np.log(TG[i - 1]), np.log(TG[i])])))


def main():
    t = tab()
    res = {"ne_cm3": NE, "atomic_data": "CHIANTI v7.1.4 O II/O III/N II; Hbeta case B (Draine 2011)",
           "OH_values": {"A09_solar": 4.90e-4, "ISM_gas_phase_Meyer98": 3.19e-4, "twice_solar": 9.8e-4}}
    # --- temperature diagnostics
    r4363 = t["e4363"] / t["e5007"]
    r7325 = t["e7325"] / t["e3727"]
    r5755 = t["n5755"] / t["n6583"]
    diag = {}
    o = OBS["O4363_5007"]
    diag["T_OIII_4363"] = {k: invert(r4363, v) for k, v in (("fid", o["fid"]), ("lo", o["lo"]), ("hi", o["hi"]),
                                                             ("desi_wave1_0.056", 0.056), ("verify_desi_0.073", 0.073))}
    a, b = OBS["OII7325_OIII"], OBS["OII_OIII"]
    vals = {"fid": a["fid"] / b["fid"], "lo": a["lo"] / b["hi"], "hi": a["hi"] / b["lo"]}
    diag["OII_7325_3727_obs"] = vals
    diag["T_OII_7325"] = {k: invert(r7325, v) for k, v in vals.items()}
    v5755 = OBS["NII5755_OIII"] / OBS["NII6583_OIII"]
    diag["NII_5755_6583_obs"] = v5755
    diag["T_NII_5755"] = invert(r5755, v5755)
    # density sensitivity of the O II auroral ratio
    t100 = tab(100.0)
    diag["T_OII_7325_ne100"] = {k: invert(t100["e7325"] / t100["e3727"], v) for k, v in vals.items()}
    res["T_diagnostics"] = diag
    # --- maximum [O III]/Hbeta and T_min
    tm = {}
    curves = {}
    for ohk, OH in res["OH_values"].items():
        for r in (OBS["OII_OIII"]["lo"], OBS["OII_OIII"]["fid"], OBS["OII_OIII"]["hi"]):
            c, fpp = max_oiii_hb(t, r, OH)
            curves[f"{ohk}_r{r}"] = c.tolist()
            for target in (OBS["OIII_Hb"]["lo"], OBS["OIII_Hb"]["fid"], OBS["OIII_Hb"]["hi"]):
                tm[f"{ohk}|r={r}|OIII/Hb={target}"] = tmin(c, target)
                c06, _ = max_oiii_hb(t, r, OH, hbscale=0.6)
                tm[f"{ohk}|r={r}|OIII/Hb={target}|Hb_x0.6"] = tmin(c06, target)
    res["T_min"] = tm
    # value at 1e4 K / 1.2e4 K for reference
    for T0 in (8e3, 1e4, 1.2e4, 1.5e4):
        i = np.argmin(abs(TG - T0))
        c, fpp = max_oiii_hb(t, OBS["OII_OIII"]["fid"], 4.9e-4)
        res[f"max_OIII_Hb_at_{int(T0)}K_r1.37_A09"] = float(c[i])
        res[f"fpp_at_{int(T0)}K_r1.37"] = float(fpp[i])
        res[f"max_R23_at_{int(T0)}K_r1.37_A09"] = float(c[i] * (1 + 1 / 2.98 + OBS["OII_OIII"]["fid"]))
        res[f"OIII_Hb_allOpp_at_{int(T0)}K_A09"] = float(4.9e-4 * t["e5007"][i] / t["hb"][i])
    res["T_grid"] = TG.tolist()
    res["ratio_4363_5007"] = r4363.tolist()
    res["ratio_OII_7325_3727"] = r7325.tolist()
    res["curves_max_OIII_Hb"] = curves
    with open(OUT, "w") as f:
        json.dump(res, f, indent=1)
    plot(t, res)
    # brief printout
    for k, v in diag.items():
        print(k, v)
    for k in list(tm)[:40]:
        print(k, None if tm[k] is None else round(tm[k], 0) if tm[k] == tm[k] else "nan")
    for k, v in res.items():
        if k.startswith("max_") or k.startswith("fpp") or k.startswith("OIII_Hb_all"):
            print(k, round(v, 3))


def plot(t, res):
    from silk_style import new_sky_figure, label, finish, OIII_C, HB_C, WARN, OK, FGDIM, ACCENT, DANGER, SII_C
    import matplotlib.pyplot as plt
    fig, axs = new_sky_figure(12, 5.2, "Silk [O III] thermometer",
                              "max [O III]/Hβ of fully ionised O⁺/O⁺⁺ gas vs T, and T from auroral lines", 1, 2)
    ax = axs[0]
    cols = {"A09_solar": OIII_C, "ISM_gas_phase_Meyer98": FGDIM, "twice_solar": WARN}
    labs = {"A09_solar": "O/H 4.9e-4", "ISM_gas_phase_Meyer98": "O/H 3.2e-4", "twice_solar": "O/H 9.8e-4"}
    for ohk in cols:
        rr = (OBS["OII_OIII"]["lo"], OBS["OII_OIII"]["fid"], OBS["OII_OIII"]["hi"])
        for r, ls in zip(rr, (":", "-", "--")):
            c = np.array(res["curves_max_OIII_Hb"][f"{ohk}_r{r}"])
            ax.plot(TG, c, color=cols[ohk], ls=ls, lw=1.6 if r == rr[1] else 1.0,
                    label=f"{labs[ohk]}" if r == rr[1] else None)
    ax.axhspan(OBS["OIII_Hb"]["lo"], OBS["OIII_Hb"]["hi"], color=OK, alpha=0.18, lw=0)
    label(ax, 5.5e3, 16.5, "observed [O III]/Hβ ≥ 15 (20–47)", color=OK, size=8)
    tt = res["T_diagnostics"]["T_OIII_4363"]
    ax.axvspan(tt["lo"], tt["hi"], color=DANGER, alpha=0.12, lw=0)
    label(ax, tt["lo"] * 1.03, 1.4, "T from 4363/5007", color=DANGER, size=8)
    ax.axvspan(8e3, 1.4e4, color=ACCENT, alpha=0.10, lw=0)
    label(ax, 8.2e3, 0.55, "photoionised\nnebulae", color=ACCENT, size=7.5)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(5e3, 2e5); ax.set_ylim(0.3, 300)
    ax.set_xlabel("electron temperature T (K)"); ax.set_ylabel("max [O III]5007/Hβ  (energy)")
    ax.legend(loc="lower right", ncol=1, title="[O II]/[O III]=1.37 (solid; 1.1 dotted, 1.9 dashed)",
              title_fontsize=7)
    ax.grid(alpha=0.3, which="both")
    ax2 = axs[1]
    ax2.plot(TG, res["ratio_4363_5007"], color=OIII_C, lw=1.6, label="[O III] 4363/5007")
    ax2.plot(TG, res["ratio_OII_7325_3727"], color=SII_C, lw=1.6, label="[O II] 7325/3727")
    ax2.plot(TG, t["n5755"] / t["n6583"], color="#c792ea", lw=1.2, label="[N II] 5755/6583")
    o = OBS["O4363_5007"]
    ax2.axhspan(o["lo"], o["hi"], color=OIII_C, alpha=0.15, lw=0)
    v = res["T_diagnostics"]["OII_7325_3727_obs"]
    ax2.axhspan(v["lo"], v["hi"], color=SII_C, alpha=0.15, lw=0)
    ax2.axhline(res["T_diagnostics"]["NII_5755_6583_obs"], color="#c792ea", ls="--", lw=0.8)
    label(ax2, 5.4e3, o["hi"] * 1.1, "DESI 4363/5007 (wave 2)", color=OIII_C, size=7.5)
    label(ax2, 5.4e3, v["hi"] * 1.15, "DESI [O II] 7325/3727", color=SII_C, size=7.5)
    ax2.set_xscale("log"); ax2.set_yscale("log"); ax2.set_xlim(5e3, 2e5); ax2.set_ylim(1e-3, 1.0)
    ax2.set_xlabel("electron temperature T (K)"); ax2.set_ylabel("line ratio (energy), n_e = 10 cm⁻³")
    ax2.legend(loc="lower right"); ax2.grid(alpha=0.3, which="both")
    finish(fig, FIG)


if __name__ == "__main__":
    main()
