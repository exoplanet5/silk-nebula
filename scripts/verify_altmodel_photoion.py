"""verify_altmodel_photoion: one-zone, optically thin photoionisation-equilibrium (PIE) model with
thermal balance, written for this project (no Cloudy).

Physics
  * elements H, He, C, N, O, Ne, S (Asplund+2009 solar); all ionisation stages of each.
  * photoionisation of ground-state outer shells: Verner & Yakovlev (1995) fits (CHIANTI verner_short.txt);
    no inner-shell / Auger ionisation (photons > 280 eV carry a negligible part of the soft spectra used).
  * collisional ionisation: Dere (2007); RR + DR: CHIANTI v7.1.4 rrparams / drparams; H: case A (thin).
  * charge exchange (approximate rates, flagged): O+ + H0 <-> O0 + H+ (k = 1e-9, reverse 8/9 k),
    O++ + H0 -> O+ + H+ (k = 1e-9 cm^3 s^-1, order of Kingdon & Ferland 1996 at 1e4 K). Others neglected.
  * heating: photo-electron energy (E - IP) for every absorbing ion.
  * cooling: collisional excitation of every level in the CHIANTI models of H I, He I-II, C II-V, N I-V, O I-VI,
    Ne II-VI, S II-VI (n_e = 1 cm^-3, low-density limit); H/He recombination cooling (Draine 2011 eq 27.23 form),
    free-free (g_ff = 1.3), H collisional-ionisation cooling.
  * optional extra (non-radiative) heating Gamma_x per H, to ask how much is needed to reach T ~ 3e4 K.
Spectral shapes (photon number per unit energy, normalised above 13.6 eV):
  BB 35 kK (OB-star / DIG proxy), BB 81.3 kK (DOZ WD J0736+1444), BB 150 kK (hottest WD/PN nuclei),
  'hot gas' exp(-E/kT)/E with kT = 0.14 keV (crude stand-in for the Monogem 1.6e6 K plasma: the real
  spectrum is line-dominated at 40-250 eV), and 'shock' exp(-E/kT)/E with kT = 28, 50 eV (precursor field of
  v_s ~ 150-200 km/s shocks, likewise crude).
The ionisation parameter is U = Phi_H / (n_H c).
Outputs: data/derived/verify_altmodel_photoion.json, data/derived/verify_altmodel_coolion.npz,
         figures/verify_altmodel_photoion.png
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
OUTJ = os.path.join(ROOT, "data/derived/verify_altmodel_photoion.json")
COOLNPZ = os.path.join(ROOT, "data/derived/verify_altmodel_coolion.npz")
FIG = os.path.join(ROOT, "figures/verify_altmodel_photoion.png")

ELEM = [1, 2, 6, 7, 8, 10, 16]
AB = {z: A.ABUND_A09[z] for z in ELEM}
TGRID = np.logspace(np.log10(3e3), np.log10(3e5), 61)
K_CX = 1.0e-9
EGRID = np.logspace(np.log10(13.6), np.log10(2000.0), 1500)
C_L, EV, K_B = A.C_L, A.EV, A.K_B


# ------------------------------------------------------------------------------------------------
def build_cooling_table():
    if os.path.exists(COOLNPZ):
        d = np.load(COOLNPZ)
        return {tuple(int(v) for v in k.split("_")[1:]): d[k] for k in d.files if k.startswith("P_")}
    tab = {}
    for Z in ELEM:
        for ion in range(1, Z + 1):
            try:
                I = A.get_ion(Z, ion)
                tab[(Z, ion)] = np.array([I.all_line_power(T, 1.0) for T in TGRID])
            except (FileNotFoundError, ValueError, IndexError) as exc:
                print("skip cooling", Z, ion, exc)
                continue
            print("cooling table", Z, ion, I.n, "levels")
    np.savez(COOLNPZ, T=TGRID, **{f"P_{Z}_{i}": v for (Z, i), v in tab.items()})
    return tab


COOL = None


def rates_at(T):
    """collisional ionisation and recombination coefficients at T for all stages."""
    C, R = {}, {}
    for Z in ELEM:
        C[Z] = np.array([A.ci_rate(Z, i, T) for i in range(1, Z + 1)])        # stage i -> i+1
        R[Z] = np.array([A.rec_rate(Z, i, T)[0] for i in range(2, Z + 2)])    # stage i+1 -> i
    return C, R


def spectrum(kind, par):
    E = EGRID
    if kind == "bb":
        kT = A.K_B * par / EV
        n = E ** 2 / np.expm1(np.clip(E / kT, 1e-8, 700))
    elif kind == "thermal":
        n = np.exp(-E / par) / E
    else:
        raise ValueError(kind)
    n = n / np.trapezoid(n, E)
    return n


def sigma_moments(f):
    """<sigma> and <sigma (E - IP)> (eV cm^2) per stage for spectrum f (normalised photon pdf)."""
    s, h = {}, {}
    for Z in ELEM:
        s[Z] = np.zeros(Z)
        h[Z] = np.zeros(Z)
        for i in range(1, Z + 1):
            xs = A.phot_xs(Z, i, EGRID)
            if not xs.any():
                continue
            ip = A.ion_potential(Z, i)
            s[Z][i - 1] = np.trapezoid(f * xs, EGRID)
            h[Z][i - 1] = np.trapezoid(f * xs * np.clip(EGRID - ip, 0, None), EGRID)
    return s, h


def ion_balance(T, ne_nH, U, smom, C, R, xH_guess=None):
    """returns dict Z -> fractions; ne in units of n_H.  photoionisation rate per n_e:
    Gamma_i / n_e = U c <sigma_i> / (n_e/n_H)."""
    frac = {}
    g = U * C_L / ne_nH  # multiply by <sigma> -> rate per n_e
    # hydrogen first (needed for CX)
    Z = 1
    ion_r = g * smom[1][0] + C[1][0]
    xH1 = ion_r / (ion_r + R[1][0])
    frac[1] = np.array([1 - xH1, xH1])
    nH0, nHp = (1 - xH1) / ne_nH, xH1 / ne_nH  # n_H0/n_e, n_H+/n_e
    for Z in ELEM[1:]:
        up = g * smom[Z] + C[Z]               # stage i -> i+1, per n_e
        dn = R[Z].copy()                      # stage i+1 -> i, per n_e
        if Z == 8:
            up[0] += (8.0 / 9.0) * K_CX * nHp     # O0 + H+ -> O+ + H0
            dn[0] += K_CX * nH0                   # O+ + H0 -> O0 + H+
            dn[1] += K_CX * nH0                   # O++ + H0 -> O+ + H+
        lx = np.zeros(Z + 1)
        for i in range(Z):
            lx[i + 1] = lx[i] + np.log(max(up[i], 1e-300)) - np.log(max(dn[i], 1e-300))
        lx -= lx.max()
        x = np.exp(lx)
        frac[Z] = x / x.sum()
    ne = sum(AB[Z] * np.dot(np.arange(Z + 1), frac[Z]) for Z in ELEM)
    return frac, ne


def cool_heat(T, frac, ne, U, smom, hmom, C, gamma_x=0.0):
    """heating and cooling per n_H^2 (erg cm^3 s^-1)."""
    heat = 0.0
    for Z in ELEM:
        heat += AB[Z] * np.dot(frac[Z][:-1], hmom[Z]) * EV
    heat *= U * C_L
    heat += gamma_x  # extra heating per n_H^2
    cool = 0.0
    for (Z, ion), P in COOL.items():
        cool += ne * AB[Z] * frac[Z][ion - 1] * np.exp(np.interp(np.log(T), np.log(TGRID), np.log(np.maximum(P, 1e-300))))
    T4 = T / 1e4
    xHp = frac[1][1]
    rec = (0.684 - 0.0416 * np.log(T4)) * K_B * T * A.alpha_B(T)
    cool += ne * xHp * rec
    cool += ne * AB[2] * frac[2][1] * rec + ne * AB[2] * frac[2][2] * 2 * (0.684 - 0.0416 * np.log(T4 / 4)) * K_B * T * A.alpha_B(T / 4)
    zz = xHp + AB[2] * (frac[2][1] + 4 * frac[2][2])
    cool += 1.42e-27 * 1.3 * np.sqrt(T) * ne * zz
    cool += ne * frac[1][0] * C[1][0] * 13.6 * EV
    return heat, cool


def solve(U, smom, hmom, gamma_x=0.0, Tfix=None):
    def state(T):
        Cr, Rr = rates_at(T)
        ne = 1.0
        for _ in range(30):
            frac, ne_new = ion_balance(T, ne, U, smom, Cr, Rr)
            if abs(ne_new - ne) < 1e-5 * ne:
                break
            ne = 0.5 * (ne + ne_new)
        return frac, ne_new, Cr

    def F(lT):
        T = np.exp(lT)
        frac, ne, Cr = state(T)
        h, c = cool_heat(T, frac, ne, U, smom, hmom, Cr, gamma_x)
        return np.log(h) - np.log(c)

    if Tfix is None:
        lo, hi = np.log(3.2e3), np.log(2.9e5)
        flo, fhi = F(lo), F(hi)
        if flo < 0:
            T = 3.2e3
        elif fhi > 0:
            T = 2.9e5
        else:
            T = float(np.exp(brentq(F, lo, hi, xtol=1e-4)))
    else:
        T = Tfix
    frac, ne, Cr = state(T)
    return T, frac, ne


def lines(T, frac, ne_nH, ne_cgs=1.0):
    """energy emissivities per n_H^2 (only ratios used)."""
    def e(name):
        return A.emis(name, T, ne_cgs)
    O = AB[8]
    out = {}
    out["Hb"] = ne_nH * frac[1][1] * A.j_hbeta(T)
    out["Ha"] = ne_nH * frac[1][1] * A.j_halpha(T)
    coll = A.h_collisional_balmer(T, ne_cgs)
    out["Hb"] += ne_nH * frac[1][0] * coll["Hb"]
    out["Ha"] += ne_nH * frac[1][0] * coll["Ha"]
    out["OIII5007"] = ne_nH * O * frac[8][2] * e("OIII5007")
    out["OIII4363"] = ne_nH * O * frac[8][2] * e("OIII4363")
    out["OII3727"] = ne_nH * O * frac[8][1] * A.oii_3727_emis(T, ne_cgs)
    out["HeII4686"] = ne_nH * AB[2] * frac[2][2] * A.j_heii4686(T)
    out["NII6583"] = ne_nH * AB[7] * frac[7][1] * e("NII6583")
    out["SII"] = ne_nH * AB[16] * frac[16][1] * (e("SII6716") + e("SII6731"))
    out["NeIII3869"] = ne_nH * AB[10] * frac[10][2] * e("NeIII3869")
    out["CIV1549"] = ne_nH * AB[6] * frac[6][3] * (e("CIV1548") + e("CIV1551"))
    return out


def ratios(L):
    o3 = L["OIII5007"]
    return {"OIII_Hb": o3 / L["Hb"], "OII_OIII": L["OII3727"] / o3, "O4363_5007": L["OIII4363"] / o3,
            "HeII_OIII": L["HeII4686"] / o3, "NII_Ha": L["NII6583"] / L["Ha"], "SII_Ha": L["SII"] / L["Ha"],
            "NeIII_OIII": L["NeIII3869"] / o3, "CIV_OIII_energy": L["CIV1549"] / o3,
            "R23": (L["OII3727"] + o3 * (1 + 1 / 2.98)) / L["Hb"]}


def main():
    global COOL
    COOL = build_cooling_table()
    spectra = {"BB35kK_DIG": ("bb", 3.5e4), "BB81kK_WDJ0736": ("bb", 8.13e4), "BB150kK": ("bb", 1.5e5),
               "hotgas_kT140eV": ("thermal", 140.0), "shock_kT28eV": ("thermal", 28.0), "shock_kT50eV": ("thermal", 50.0)}
    logU = np.arange(-5.0, -0.49, 0.25)
    res = {"logU": logU.tolist(), "models": {}, "notes": __doc__}
    for name, (kind, par) in spectra.items():
        f = spectrum(kind, par)
        smom, hmom = sigma_moments(f)
        rows = []
        for lu in logU:
            T, frac, ne = solve(10 ** lu, smom, hmom)
            L = lines(T, frac, ne)
            r = ratios(L)
            r.update(T=T, xHp=float(frac[1][1]), xOpp=float(frac[8][2]), xOp=float(frac[8][1]),
                     xHepp=float(frac[2][2]), ne_nH=float(ne), logU=float(lu))
            rows.append({k: float(v) for k, v in r.items()})
        res["models"][name] = rows
        best = max(rows, key=lambda r: r["OIII_Hb"])
        print(f"{name:18s} max OIII/Hb={best['OIII_Hb']:.2f} at logU={best['logU']:.2f} T={best['T']:.0f}"
              f" OII/OIII={best['OII_OIII']:.2f} 4363/5007={best['O4363_5007']:.4f} HeII/OIII={best['HeII_OIII']:.3f}")
        # where OII/OIII ~ 1.3
        rr = np.array([r["OII_OIII"] for r in rows])
        k = np.argmin(abs(np.log(rr / 1.3)))
        r = rows[k]
        print(f"     at OII/OIII~1.3: logU={r['logU']:.2f} T={r['T']:.0f} OIII/Hb={r['OIII_Hb']:.2f} R23={r['R23']:.1f}"
              f" NII/Ha={r['NII_Ha']:.2f} SII/Ha={r['SII_Ha']:.2f} HeII/OIII={r['HeII_OIII']:.3f}")
    # sensitivity: CX rate x0 and x3 for the WD spectrum
    global K_CX
    sens = {}
    for kcx in (0.0, 3e-9):
        K_CX = kcx
        f = spectrum("bb", 8.13e4)
        smom, hmom = sigma_moments(f)
        rows = []
        for lu in (-3.0, -2.0, -1.0):
            T, frac, ne = solve(10 ** lu, smom, hmom)
            r = ratios(lines(T, frac, ne))
            rows.append(dict(logU=lu, T=T, **{k: float(v) for k, v in r.items()}))
        sens[f"K_CX={kcx}"] = rows
    K_CX = 1.0e-9
    res["sensitivity_CX_WD"] = sens
    # extra heating needed to reach T = 2e4, 3e4 K at the U where [O II]/[O III] ~ 1.3 (WD spectrum)
    f = spectrum("bb", 8.13e4)
    smom, hmom = sigma_moments(f)
    extra = {}
    for lu in (-2.5, -2.0, -1.5):
        U = 10 ** lu
        T0, frac0, ne0 = solve(U, smom, hmom)
        h0, c0 = cool_heat(T0, frac0, ne0, U, smom, hmom, rates_at(T0)[0])
        for Tt in (2e4, 3e4):
            Tq, fr, ne = solve(U, smom, hmom, Tfix=Tt)
            h, c = cool_heat(Tt, fr, ne, U, smom, hmom, rates_at(Tt)[0])
            r = ratios(lines(Tt, fr, ne))
            extra[f"logU={lu}|T={int(Tt)}"] = dict(T_eq=T0, photoheat_at_Teq=h0, cooling_needed=c, photoheat=h,
                                                    extra_over_photo=float((c - h) / h),
                                                    **{k: float(v) for k, v in r.items()})
    res["extra_heating_WD"] = extra
    with open(OUTJ, "w") as fo:
        json.dump(res, fo, indent=1)
    for k, v in extra.items():
        print(k, {kk: round(vv, 3) if isinstance(vv, float) else vv for kk, vv in v.items() if kk in ("T_eq", "extra_over_photo", "OIII_Hb", "OII_OIII", "O4363_5007")})
    plot(res)


def plot(res):
    from silk_style import new_sky_figure, label, finish, OIII_C, WARN, OK, FGDIM, ACCENT, DANGER, SII_C, HA_C
    fig, axs = new_sky_figure(12.5, 5.4, "Photoionisation equilibrium vs the Silk spectrum",
                              "one-zone optically thin PIE + thermal balance; solar abundances", 1, 2)
    cols = {"BB35kK_DIG": FGDIM, "BB81kK_WDJ0736": OIII_C, "BB150kK": "#2979ff", "hotgas_kT140eV": DANGER,
            "shock_kT28eV": WARN, "shock_kT50eV": "#c792ea"}
    labs = {"BB35kK_DIG": "35 kK BB (DIG)", "BB81kK_WDJ0736": "81 kK BB (WD J0736)", "BB150kK": "150 kK BB",
            "hotgas_kT140eV": "0.14 keV hot gas", "shock_kT28eV": "precursor kT 28 eV", "shock_kT50eV": "precursor kT 50 eV"}
    ax, ax2 = axs
    for k, rows in res["models"].items():
        x = [r["OII_OIII"] for r in rows]
        y = [r["OIII_Hb"] for r in rows]
        ax.plot(x, y, "-o", ms=2.5, color=cols[k], lw=1.3, label=labs[k])
        ax2.plot([r["logU"] for r in rows], [r["T"] for r in rows], "-", color=cols[k], lw=1.4, label=labs[k])
    ax2.legend(loc="lower right", fontsize=7)
    ax.add_patch(__import__("matplotlib").patches.Rectangle((OBS["OII_OIII"]["lo"], OBS["OIII_Hb"]["lo"]),
                 OBS["OII_OIII"]["hi"] - OBS["OII_OIII"]["lo"], OBS["OIII_Hb"]["hi"] - OBS["OIII_Hb"]["lo"],
                 color=OK, alpha=0.3, lw=0))
    label(ax, 2.0, 30, "Silk (DESI/SDSS)", color=OK, size=8.5)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.01, 300); ax.set_ylim(0.01, 100)
    ax.set_xlabel("[O II]3727 / [O III]5007"); ax.set_ylabel("[O III]5007 / Hβ")
    ax.legend(loc="lower left", fontsize=7.5); ax.grid(alpha=0.3, which="both")
    t = OBS["O4363_5007"]
    ax2.axhspan(31000, 53000, color=DANGER, alpha=0.12, lw=0)
    label(ax2, -4.9, 3.5e4, "T_e([O III]4363), DESI wave 2", color=DANGER, size=8)
    ax2.axhline(1.3e4, color=OK, ls="--", lw=0.8)
    label(ax2, -4.9, 1.33e4, "T_min for [O III]/Hβ ≥ 15", color=OK, size=7.5)
    ax2.set_yscale("log"); ax2.set_ylim(3e3, 6e4)
    ax2.set_xlabel("log U"); ax2.set_ylabel("equilibrium T (K)")
    ax2.grid(alpha=0.3, which="both")
    finish(fig, FIG)


if __name__ == "__main__":
    main()
