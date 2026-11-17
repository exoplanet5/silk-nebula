"""verify_altmodel_summary: scenario scorecard for non-shock origins of the Silk [O III] filament.

Collects the predictions of verify_altmodel_{thermo,photoion,cooling,budget}.json, adds the two
energetics numbers that need both (TML entrainment speed, X-ray brightness of a 1.6 MK cooling flow
relative to the Monogem shell), scores every scenario against the observed constraints, and writes
data/derived/verify_altmodel_scenarios.json and figures/verify_altmodel_scorecard.png.
Scores: 2 = within the observed range, 1 = within ~1.5x of it (or needs one ad hoc assumption),
0 = excluded.  Run: ~/.venvs/astro313/bin/python -I scripts/verify_altmodel_summary.py
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from verify_altmodel_common import OBS, sheet_G  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
D = lambda n: json.load(open(os.path.join(ROOT, "data/derived", n)))  # noqa: E731
OUTJ = os.path.join(ROOT, "data/derived/verify_altmodel_scenarios.json")
FIG = os.path.join(ROOT, "figures/verify_altmodel_scorecard.png")
K_B, PC = 1.380649e-16, 3.0857e18
NE_SCALE_A09 = (8.51e-5 / 4.90e-4) / (1.95e-4 / 4.57e-4)   # GS07 Ne/O -> Asplund09 Ne/O


def band(v, lo, hi, f=1.5):
    if v is None or not np.isfinite(v):
        return 0
    if lo <= v <= hi:
        return 2
    if lo / f <= v <= hi * f:
        return 1
    return 0


def upper(v, lim, f=1.5):
    return 2 if v <= lim else (1 if v <= lim * f else 0)


def lower(v, lim, f=1.5):
    return 2 if v >= lim else (1 if v >= lim / f else 0)


def main():
    cool, bud, pj = D("verify_altmodel_cooling.json"), D("verify_altmodel_budget.json"), D("verify_altmodel_photoion.json")
    G, delta = sheet_G()
    P_mono = 1.5e4 * K_B
    out = {"G": G, "delta_pc": delta, "scenarios": {}}
    # ---- energetics that combine files
    cf = cool["cooling_flows"]["isobaric_Z1"]["T1.6e+06-1e+04"]
    v_hot = cf["energy_flux"] / (2.5 * P_mono)
    out["TML_or_condensation_inflow_speed_kms_at_Monogem_P"] = v_hot / 1e5
    out["TML_I5007_R_if_v_in_10kms"] = 2.0 * 1e6 / v_hot
    # X-ray: energy radiated above 1e6 K in isobaric cooling from 1.6 MK ~ (1.6-1.0)/(1.6-0.01) of enthalpy
    f_hot = (1.6e6 - 1.0e6) / (1.6e6 - 1.0e4)
    Lam16 = bud["Lambda_CIE_1.6MK_GS07"]
    xr = {}
    for chord in (22.0, 48.0):
        for ne in (3.1e-3, 7.2e-3):
            mono = ne * ne / 1.2 * Lam16 * chord * PC
            xr[f"chord{chord:g}pc_ne{ne:.1e}"] = G * f_hot * cf["energy_flux"] / mono
    out["coolingflow_Xray_filament_over_Monogem"] = xr
    print("TML inflow speed needed (km/s):", round(v_hot / 1e5), "; I5007 at v_in=10 km/s:", round(2.0 * 1e6 / v_hot, 3), "R")
    print("1.6 MK cooling-flow X-ray filament / Monogem surface brightness:", {k: f"{v:.2g}" for k, v in xr.items()})
    # ---- scenario predictions (representative cases)
    sc = {}

    def add(name, src, p, budget, extra_fail=None, note=""):
        sc[name] = dict(source=src, pred=p, budget=budget, extra=extra_fail or {}, note=note)

    ev = cool["evaporation"]["runs"]["R60pc"]
    add("1a conduction front, evaporating (R_cl = 60 pc, best)", "own NEI steady front",
        dict(OIII_Hb=ev["OIII_Hb"], OII_OIII=ev["OII_OIII"], O4363=ev["O4363_5007"], O7325=ev["OII_7325_3727"],
             HeII=ev["HeII_OIII"], NeIII=ev["NeIII_OIII"], FUV=ev["FUV_CU_per_R"], NUV=ev["NUV_CU_per_R"], NII_Ha=ev["NII_Ha"]),
        dict(score=0, text=f"face-on {ev['I5007_faceon_R']:.1e} R: needs {ev['fronts_needed_for_2R_at_G']:.0f} fronts at G = {G:.0f}"),
        dict(density=0, xray=2), "ions under-ionised: O++ at 1-2e5 K")
    st = cool["static_fronts"]["P1.5e+04"]
    add("1b conduction front, static radiative (P/k 1.5e4)", "own CIE conductive-radiative front",
        dict(OIII_Hb=st["OIII_Hb"], OII_OIII=st["OII_OIII"], O4363=st["O4363_5007"], O7325=st["OII_7325_3727"],
             HeII=st["HeII_OIII"], NeIII=st["NeIII_OIII"] * NE_SCALE_A09, FUV=st["FUV_CU_per_R"], NUV=st["NUV_CU_per_R"], NII_Ha=st["NII_Ha"]),
        dict(score=0, text=f"face-on {st['I5007_faceon_R']:.1e} R: {st['fronts_needed_for_2R_at_G']:.0f} fronts; 86 pc thick"),
        dict(density=0, xray=1), f"N(O VI) = {st['cols']['OVI']:.1e} cm^-2 per front")
    for key, lab in (("T1.6e+06-1e+04", "2 TML / 3 NEI cooling of 1.6 MK gas, to 1e4 K"),
                     ("T1.6e+06-1e+04+tail", "3' same + H recombination tail (complete)")):
        c = cool["cooling_flows"]["isobaric_Z1"][key]
        add(lab, "GS07 isobaric Z = 1 + own emissivities",
            dict(OIII_Hb=c["OIII_Hb"], OII_OIII=c["OII_OIII"], O4363=c["O4363_5007"], O7325=c["OII_7325_3727"],
                 HeII=c["HeII_OIII"], NeIII=c["NeIII_OIII"] * NE_SCALE_A09, FUV=c["FUV_CU_per_R"], NUV=c["NUV_CU_per_R"], NII_Ha=c["NII_Ha"]),
            dict(score=0, text=f"needs {c['energy_flux']:.1e} erg/cm2/s = inflow {v_hot / 1e5:.0f} km/s at Monogem P"),
            dict(density=0, xray=0), "f(5007) = %.1e of the cooling" % c["f5007"])
    c3 = cool["cooling_flows"]["isobaric_Z1"]["T3.0e+05-1e+04"]
    add("3'' NEI cooling from 3e5 K (no hot phase)", "GS07 isobaric Z = 1, T < 3e5 K part",
        dict(OIII_Hb=c3["OIII_Hb"], OII_OIII=c3["OII_OIII"], O4363=c3["O4363_5007"], O7325=c3["OII_7325_3727"],
             HeII=c3["HeII_OIII"], NeIII=c3["NeIII_OIII"] * NE_SCALE_A09, FUV=c3["FUV_CU_per_R"], NUV=c3["NUV_CU_per_R"], NII_Ha=c3["NII_Ha"]),
        dict(score=1, text=f"needs {c3['energy_flux']:.1e} erg/cm2/s: what heated it to 3e5 K? (= a ~150 km/s shock)"),
        dict(density=1, xray=2), "this is the cooling zone of an incomplete shock")
    w = cool["warm_CIE_window_OII_OIII_1_2"][-1]
    add("3w warm CIE layer held at 6e4 K", "GS07 CIE + own emissivities",
        dict(OIII_Hb=w["OIII_Hb"], OII_OIII=w["OII_OIII"], O4363=w["O4363_5007"], O7325=w["OII_7325_3727"],
             HeII=w["HeII_OIII"], NeIII=w["NeIII_OIII"] * NE_SCALE_A09, FUV=w["FUV_CU_per_R"], NUV=w["NUV_CU_per_R"], NII_Ha=w["NII_Ha"]),
        dict(score=1, text="heating 1.7e-5 erg/cm2/s, no identified source; thermally unstable"),
        dict(density=1, xray=2), "single T cannot give T(O II) = 1.2e4 K")
    for sname, lab in (("BB81kK_WDJ0736", "4a PIE: WD J0736 (81 kK)"), ("hotgas_kT140eV", "4b PIE: Monogem EUV (0.14 keV)"),
                       ("BB35kK_DIG", "4c PIE: ISRF / DIG (35 kK)"), ("shock_kT28eV", "5 shock precursor (PIE, 28 eV)")):
        rows = pj["models"][sname]
        oo = np.array([r["OII_OIII"] for r in rows])
        k = int(np.argmin(abs(np.log(oo / 1.37))))
        r = rows[k]
        hy = {"BB81kK_WDJ0736": "BB81kK", "hotgas_kT140eV": "hotgas140eV", "shock_kT28eV": "precursor28eV"}.get(sname)
        fuv = nuv = o7325 = None
        if hy:   # equilibrium-T row of the finer (dlogU = 0.125, incl. Si) grid at [O II]/[O III] ~ 1.37
            r = [x for x in bud["PIE_plus_heating"][hy] if x["Tfix"] is None][0]
            fuv, nuv, o7325 = r["FUV"], r["NUV"], r["OII_7325_3727"]
        srcs = {"BB81kK_WDJ0736": "WDJ0736_coplanar_6.4pc", "hotgas_kT140eV": "Monogem_ne7.2e-03_Emean30eV"}
        if sname in srcs:
            s = bud["sources"][srcs[sname]]
            btxt = f"U short by {s['shortfall_at_nrms']:.0f}x (n_rms) to {s['shortfall_at_doublet']:.0e}x (n from [O II])"
        elif sname == "BB35kK_DIG":
            btxt = "field soft: max [O III]/Hb = 1.3 at any U"
        else:
            btxt = "needs the parent shock; pre-shock n ~ 1 gives the low-density [O II] doublet"
        add(lab, "own one-zone PIE + thermal balance",
            dict(OIII_Hb=r["OIII_Hb"], OII_OIII=r["OII_OIII"], O4363=r["O4363_5007"], O7325=o7325, HeII=r["HeII_OIII"],
                 NeIII=r["NeIII_OIII"], FUV=fuv, NUV=nuv, NII_Ha=r["NII_Ha"]),
            dict(score=0, text=btxt), dict(density=0, xray=2), f"T_eq = {r['T']:.0f} K at log U = {r['logU']:.2f}")
    hh = [x for x in bud["PIE_plus_heating"]["hotgas140eV"] if x["Tfix"] == 4e4][0]
    add("4d PIE + extra heating to 4e4 K", "own PIE at fixed T",
        dict(OIII_Hb=hh["OIII_Hb"], OII_OIII=hh["OII_OIII"], O4363=hh["O4363_5007"], O7325=hh["OII_7325_3727"], HeII=hh["HeII_OIII"],
             NeIII=hh["NeIII_OIII"], FUV=hh["FUV"], NUV=hh["NUV"], NII_Ha=hh["NII_Ha"]),
        dict(score=0, text=f"heating {hh['heat_extra_over_photo']:.0f}x photo-heating, plus the U deficit of 4a/4b"),
        dict(density=1, xray=2), "")
    fo = bud["fossil"]["T0=4e+04|nH=100.0"]
    rows = fo["rows"]
    mid = [r for r in rows if fo["matching_window_yr"] and fo["matching_window_yr"][0] <= r["t_yr"] <= fo["matching_window_yr"][1]]
    m = mid[len(mid) // 2] if mid else rows[0]
    add("6 fossil, recombining (hot start 4e4 K, n_H = 100)", "own NEI recombination",
        dict(OIII_Hb=m["OIII_Hb"], OII_OIII=m["OII_OIII"], O4363=m["O4363_5007"], O7325=m["OII_7325_3727"], HeII=0.0,
             NeIII=None, FUV=None, NUV=None, NII_Ha=None),
        dict(score=0, text=f"fades by e in {fo['t_fade_yr']:.0f} yr; matching window {fo['matching_window_yr'][0]:.0f}-{fo['matching_window_yr'][1]:.0f} yr"),
        dict(density=2, xray=2), "needs a heating/ionising event < ~20 yr ago")
    # ---- scoring
    cols = ["OIII_Hb", "OII_OIII", "O4363", "O7325", "HeII", "NeIII", "FUV", "NUV", "NII_Ha", "density", "xray", "budget"]
    for name, s in sc.items():
        p = s["pred"]
        sco = {}
        sco["OIII_Hb"] = lower(p["OIII_Hb"], OBS["OIII_Hb"]["lo"])
        sco["OII_OIII"] = band(p["OII_OIII"], OBS["OII_OIII"]["lo"], OBS["OII_OIII"]["hi"])
        sco["O4363"] = band(p["O4363"], OBS["O4363_5007"]["lo"], OBS["O4363_5007"]["hi"], f=1.25)
        v = OBS["OII7325_OIII"]
        sco["O7325"] = band(p["O7325"], v["lo"] / OBS["OII_OIII"]["hi"], v["hi"] / OBS["OII_OIII"]["lo"]) if p["O7325"] is not None else None
        sco["HeII"] = upper(p["HeII"], OBS["HeII_OIII_max"])
        sco["NeIII"] = band(p["NeIII"], OBS["NeIII_OIII"]["lo"], OBS["NeIII_OIII"]["hi"]) if p["NeIII"] is not None else None
        sco["FUV"] = band(p["FUV"], OBS["FUV_CU_per_R"]["lo"], OBS["FUV_CU_per_R"]["hi"]) if p["FUV"] is not None else None
        sco["NUV"] = band(p["NUV"], OBS["NUV_CU_per_R"]["lo"], OBS["NUV_CU_per_R"]["hi"]) if p["NUV"] is not None else None
        sco["NII_Ha"] = band(p["NII_Ha"], OBS["NII_Ha"]["lo"], OBS["NII_Ha"]["hi"]) if p["NII_Ha"] is not None else None
        sco["density"] = s["extra"]["density"]
        sco["xray"] = s["extra"]["xray"]
        sco["budget"] = s["budget"]["score"]
        s["scores"] = sco
        s["n_fail"] = int(sum(1 for c in cols if sco.get(c) == 0))
        print(f"{name:52s} fails={s['n_fail']:2d} " + " ".join(f"{c}={sco[c]}" for c in cols))
    out["scenarios"] = sc
    out["columns"] = cols
    with open(OUTJ, "w") as f:
        json.dump(out, f, indent=1, default=float)
    plot(out)


def plot(out):
    from silk_style import new_sky_figure, finish, label, FG, FGDIM
    from matplotlib.colors import ListedColormap
    sc = out["scenarios"]
    cols = out["columns"]
    names = list(sc)
    heads = ["[O III]\n/Hβ\n≥15", "[O II]\n/[O III]\n1.1–1.9", "4363\n/5007\n.06–.10", "7325\n/3727\n.02–.04",
             "He II\n/[O III]\n<.026", "[Ne III]\n/[O III]\n.06–.11", "FUV\nCU/R\n68–92", "NUV\nCU/R\n36–80",
             "[N II]\n/Hα\n1.0–1.4", "n_e\n[O II]\n~100", "no\nX-ray\nfil.", "energy/\nphoton\nbudget"]
    fmt = {"OIII_Hb": "{:.3g}", "OII_OIII": "{:.2g}", "O4363": "{:.3f}", "O7325": "{:.3f}", "HeII": "{:.3f}", "NeIII": "{:.2f}",
           "FUV": "{:.0f}", "NUV": "{:.0f}", "NII_Ha": "{:.2g}"}
    M = np.full((len(names), len(cols)), np.nan)
    for i, n in enumerate(names):
        for j, c in enumerate(cols):
            v = sc[n]["scores"].get(c)
            if v is not None:
                M[i, j] = v
    fig, ax = new_sky_figure(15.0, 7.8, "Non-shock origins of the Silk [O III]: scorecard",
                             "green = within observed range, amber = within ~1.5× or one ad hoc assumption, red = excluded; "
                             "cell text = model value (energy ratios)")
    cmap = ListedColormap(["#7a2b2b", "#7a5a1f", "#2f6b3a"])
    ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=-0.5, vmax=2.5, aspect="auto")
    for i, n in enumerate(names):
        for j, c in enumerate(cols):
            if c in fmt:
                v = sc[n]["pred"].get(c)
                t = "n/a" if v is None else fmt[c].format(v)
            else:
                t = {0: "✗", 1: "~", 2: "✓"}.get(sc[n]["scores"].get(c), "")
            ax.text(j, i, t, ha="center", va="center", fontsize=7.3, color=FG if np.isfinite(M[i, j]) else FGDIM)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(heads, fontsize=6.8)
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels([n if len(n) < 50 else n[:49] + "…" for n in names], fontsize=7.5)
    import textwrap
    for i, n in enumerate(names):
        t = "\n".join(textwrap.wrap(sc[n]["budget"]["text"], 46))
        label(ax, len(cols) - 0.4, i, t, color=FGDIM, size=6.3, va="center", ha="left", clip_on=False)
    ax.set_xlim(-0.5, len(cols) - 0.5)
    fig.subplots_adjust(left=0.235, right=0.80, top=0.84, bottom=0.03)
    finish(fig, FIG, tight=False)


if __name__ == "__main__":
    main()
