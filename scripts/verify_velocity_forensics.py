"""verify_velocity: forensic enumeration of blue-arm wavelength-reference mistakes for LAMOST MRS-N [O III].

What would each plausible mistake in the MRS-N blue arm (4950-5350 A; only sky line = [N I] 5198/5200 doublet)
do to a 5007 velocity?  Also: extrapolation of the red-arm (6290-6560 A) quadratic RV-calibration functions of
Ren et al. 2021 (RAA 21, 51; arXiv:2008.06236, Table 6) to 5008 A, and the LSR conversion used by the paper.

All numbers here are derivations from stated rest wavelengths (no data fitting except the doublet simulation).
Rest wavelengths:
  [N I] 2D3/2 -> 4S3/2 and 2D5/2 -> 4S3/2: vacuum from NIST level energies 19233.177 and 19224.464 cm^-1
      (=> 5199.349, 5201.705 A vac; air 5197.90, 5200.26).  Hanuschik 2003 UVES atlas (VizieR J/A+A/407/1157,
      table 6, file data/raw/vizier_hanuschik2003/table6.tsv): 5197.9282, 5200.2856 A (air), peak fluxes 4.7512,
      2.7566 (FWHM 0.1219, 0.1211 A) -> intensity ratio 1.73.
  [O III] 5006.843 air / vacuum by Morton (2000).
  OH 6553.617 air / 6555.4275 vac (Ren et al. 2021, Table 2); Hanuschik measures 6553.6260 (FWHM 0.239 A, blend).
LAMOST MRS resolution R ~ 7500 (data-production description; Ren+2021).
Heliocentric correction on 2024-02-01 (astropy, Xinglong, field centre): v_corr = -8.25..-8.38 km/s
(data/derived/verify_velocity_lamost_arms.json); LAMOST header HELIO_RV = +8.1..+8.3 = -v_corr.

Run: ~/.venvs/astro313/bin/python -I scripts/verify_velocity_forensics.py
Output: data/derived/verify_velocity_forensics.json
"""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit
from astropy.coordinates import SkyCoord, ICRS, Galactic
import astropy.units as u

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/derived/verify_velocity_forensics.json'
C = 299792.458


def n_air(lam_vac):
    s2 = (1e4 / np.asarray(lam_vac, float)) ** 2
    return 1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2)


def vac2air(lv):
    return lv / n_air(lv)


def air2vac(la):
    v = np.asarray(la, float)
    for _ in range(5):
        v = la * n_air(v)
    return v


NI_VAC = np.array([1e8 / 19233.177, 1e8 / 19224.464])          # 5198 (2D3/2), 5200 (2D5/2)
NI_AIR = vac2air(NI_VAC)
NI_HAN_AIR = np.array([5197.9282, 5200.2856])
NI_HAN_RATIO = (4.7512 * 0.1219) / (2.7566 * 0.1211)
O3_AIR = 5006.843
O3_VAC = float(air2vac(O3_AIR))
OH_AIR, OH_VAC_REN = 6553.617, 6555.4275
VCORR = -8.32          # km/s, exposure-weighted mean of the three exposures (astropy barycentric)


def doublet_centroid(ratio, lam=NI_VAC):
    return (ratio * lam[0] + lam[1]) / (1 + ratio)


def scenarios():
    """Velocity error at 5007 for combinations of (what is measured, what reference is assigned, [O III] rest).

    Model: LAMOST 1D spectra are vacuum and heliocentric (HELIO=T, VACUUM=T).  A sky line at topocentric rest
    lam0 therefore appears at lam0*(1+v_corr/c).  A secondary calibration that forces the measured sky feature to
    lam_ref rescales the wavelength axis by k = lam_ref/lam_meas (multiplicative; additive shifts differ by
    <= 3 km/s at 5007 and are listed separately).  The [O III] velocity is then c*(k*lam_O*(1+V/c)/rest - 1),
    optionally followed by a further heliocentric term.
    """
    out = []
    meas_opts = {
        'strong [N I] component (resolved)': lambda r: NI_VAC[0],
        'single-Gaussian/moment of whole doublet': doublet_centroid,
    }
    ref_opts = {
        'vac 5199.349 (correct)': NI_VAC[0],
        'Hanuschik air 5197.928': NI_HAN_AIR[0],
        'Hanuschik air->vac 5199.376': float(air2vac(NI_HAN_AIR[0])),
        'label 5198.000': 5198.0,
        'label 5199.000': 5199.0,
        'NIST air 5197.902': NI_AIR[0],
        'unweighted doublet mean vac 5200.527': NI_VAC.mean(),
        'unweighted doublet mean air 5199.08': NI_AIR.mean(),
        'Hanuschik-ratio (1.73) centroid vac': float(doublet_centroid(NI_HAN_RATIO)),
    }
    rest_opts = {'[O III] vac 5008.240': O3_VAC, '[O III] air 5006.843': O3_AIR}
    frame_opts = {
        'sky forced to 0 in heliocentric-corrected spectra (no further helio)': (True, 0.0),
        'sky line measured in topocentric frame, helio applied correctly': (False, 0.0),
        'sky forced to 0 in helio frame, then helio applied again (correct sign)': (True, +1.0),
        'sky forced to 0 in helio frame, then helio applied with wrong sign': (True, -1.0),
    }
    V0 = 0.0
    for mname, mf in meas_opts.items():
        for ratio in ([1.4, 1.73, 1.9] if 'doublet' in mname else [1.73]):
            for rname, lref in ref_opts.items():
                for oname, rest in rest_opts.items():
                    for fname, (in_helio, extra) in frame_opts.items():
                        lam_sky = mf(ratio) * ((1 + VCORR / C) if in_helio else 1.0)
                        k = lref / lam_sky
                        # heliocentric-corrected spectrum: [O III] at V0; topocentric spectrum: at V0 - v_corr
                        lam_o = O3_VAC * (1 + (V0 if in_helio else V0 - VCORR) / C)
                        # if not in_helio we work in topocentric frame then add v_corr at the end
                        v = C * (k * lam_o / rest - 1)
                        if not in_helio:
                            v += VCORR
                        v += extra * VCORR
                        # additive-shift variant
                        lam_o_add = lam_o + (lref - lam_sky)
                        v_add = C * (lam_o_add / rest - 1) + (0 if in_helio else VCORR) + extra * VCORR
                        out.append(dict(measured=mname, ratio=ratio, reference=rname, o3_rest=oname, frame=fname,
                                        dV_mult=round(float(v - V0), 2), dV_add=round(float(v_add - V0), 2)))
    return out


def ren_extrapolation():
    """Ren+2021 Table 6 quadratic RV-calibration functions f(x)=a x^2+b x+c (x in micron), fitted on red-arm sky
    lines 6287-6553 A (+ an extrapolated anchor at 6731).  Evaluate at 5008 (if wrongly applied to the blue arm)
    and at 6585 ([N II]) for comparison."""
    tab = [  # date, plate, a, b, c
        ('20180305', 'HD474690101', -738.7490, 989.4288, -334.8446),
        ('20180305', 'HD474690801', -835.4372, 1116.1746, -376.2256),
        ('20181018', 'NT054437N290059N01', -5739.5264, 7601.2190, -2518.0593),
        ('20181118', 'NT070341S003325N01', 742.0199, -991.7874, 331.5405),
        ('20181128', 'NT230608N631246N01', -135.3334, 182.7665, -61.9151),
        ('20181128', 'NT235302N592517N01', 118.2756, -154.3496, 50.2610),
        ('20181129', 'NT010552N655815N01', 1153.0985, -1519.0363, 500.4621),
        ('20181129', 'NT232020N601629N01', 571.4184, -753.3775, 248.4534),
        ('20181216', 'NT065450N031415N01', -425.1480, 547.7008, -176.2284),
        ('20181228', 'NT023239N614403N01', -663.7146, 864.6065, -281.2221),
        ('20190124', 'NT024709N603414N01', -311.3390, 408.5562, -134.5974),
        ('20190125', 'NT041422N543127N01', 69.5751, -99.9254, 35.4064),
        ('20190125', 'NT041453N482433N01', -158.2488, 197.4298, -61.5766),
        ('20190126', 'NT054421N274353N01', 58.2406, -88.3461, 32.6192),
        ('20190126', 'NT054421N274353N02', 68.5436, -103.3951, 37.9632),
    ]
    rows = []
    for d, p, a, b, c in tab:
        f = lambda x: a * x * x + b * x + c
        rows.append(dict(date=d, plate=p, f_5008=round(f(0.5008), 2), f_6565=round(f(0.6565), 2),
                         f_6585=round(f(0.6585), 2), blue_minus_red=round(f(0.5008) - f(0.6585), 2)))
    bm = np.array([r['blue_minus_red'] for r in rows])
    late = np.array([r['blue_minus_red'] for r in rows if r['date'] >= '20181118'])
    return dict(rows=rows, blue_minus_red_all=dict(median=float(np.median(bm)), min=float(bm.min()), max=float(bm.max()),
                                                   std=float(bm.std(ddof=1))),
                blue_minus_red_after_2018Nov=dict(median=float(np.median(late)), min=float(late.min()),
                                                  max=float(late.max()), std=float(late.std(ddof=1))),
                note='Sign convention of Ren: f is the systematic sky RV that is subtracted; a blue-arm application '
                     'would change V([O III]) by -(f(0.5008)) relative to applying nothing.')


def doublet_fit_sim():
    """Fit a single Gaussian + quadratic to the [N I] doublet at R=7500 in windows of different half-width,
    starting at the nominal 5199 A, to see what a careless single-line fit returns."""
    R = 7500.
    sig = NI_VAC[0] / R / 2.3548
    res = []
    rng = np.random.default_rng(3)
    for ratio in (1.4, 1.73, 1.9):
        for half in (1.5, 2.0, 3.0, 4.0, 6.0):
            for p0c in (5199.35, 5200.0, 5200.5):
                x = np.arange(5199.5 - half, 5199.5 + half, 0.1)   # ~0.1 A/pix (MRS blue 4117 pix / 400 A)
                y = ratio * np.exp(-0.5 * ((x - NI_VAC[0]) / sig) ** 2) + np.exp(-0.5 * ((x - NI_VAC[1]) / sig) ** 2)
                y = y + 0.05 + rng.normal(0, 0.02, x.size)

                def model(xx, A, mu, s, c0, c1, c2):
                    t = xx - 5199.5
                    return A * np.exp(-0.5 * ((xx - mu) / s) ** 2) + c0 + c1 * t + c2 * t * t
                try:
                    p, _ = curve_fit(model, x, y, p0=[1.5, p0c, sig * 1.2, 0.05, 0, 0], maxfev=20000)
                    mu = p[1]
                    s = abs(p[2])
                except Exception:
                    mu, s = np.nan, np.nan
                res.append(dict(ratio=ratio, half_window=half, start=p0c, mu=round(float(mu), 3),
                                sigma=round(float(s), 3), dv_vs_strong=round(float(C * (mu - NI_VAC[0]) / NI_VAC[0]), 1)))
    return res


def lsr_terms():
    c = SkyCoord(116.45 * u.deg, 16.35 * u.deg)
    # paper: standard solar motion 20 km/s toward RA 18h03m50.29s Dec +30d00m16.8s (J2000)
    apex = SkyCoord('18h03m50.29s', '+30d00m16.8s')
    std = 20.0 * np.cos(c.separation(apex).rad)
    # Schoenrich+2010 (U,V,W)=(11.1,12.24,7.25)
    l, b = c.galactic.l.rad, c.galactic.b.rad
    U, V, W = 11.1, 12.24, 7.25
    sch = U * np.cos(l) * np.cos(b) + V * np.sin(l) * np.cos(b) + W * np.sin(b)
    # IAU kinematic LSR (20 km/s to 18h, +30 B1900) as commonly used in radio (WHAM uses it too): same as std
    return dict(Vlsr_minus_Vhelio_standard=float(std), Vlsr_minus_Vhelio_Schoenrich2010=float(sch),
                l_deg=float(np.degrees(l)), b_deg=float(np.degrees(b)))


def main():
    sc = scenarios()
    out = dict(
        rest=dict(NI_vac=NI_VAC.round(3).tolist(), NI_air_from_vac=NI_AIR.round(3).tolist(),
                  NI_Hanuschik_air=NI_HAN_AIR.tolist(), NI_Hanuschik_ratio=round(NI_HAN_RATIO, 3),
                  NI_Hanuschik_minus_NIST_kms=(C * (NI_HAN_AIR - NI_AIR) / NI_AIR).round(2).tolist(),
                  O3_vac=round(O3_VAC, 3), OH6553_vac=round(float(air2vac(OH_AIR)), 4), OH6553_vac_Ren=OH_VAC_REN,
                  OH6553_Hanuschik_minus_Ren_kms=round(C * (6553.6260 - OH_AIR) / OH_AIR, 2),
                  doublet_sep_A=round(float(NI_VAC[1] - NI_VAC[0]), 3),
                  LSF_FWHM_A_at_5199=round(float(NI_VAC[0] / 7500), 3),
                  doublet_centroid_vac={str(r): round(float(doublet_centroid(r)), 3) for r in (1.4, 1.6, 1.73, 1.9)}),
        vcorr=VCORR,
        lsr=lsr_terms(),
        scenarios=sc,
        ren_extrapolation=ren_extrapolation(),
        doublet_fit_sim=doublet_fit_sim(),
    )
    # summary of scenarios that land within the observed offset window (+8..+18 km/s)
    hits = [s for s in sc if 8 <= s['dV_mult'] <= 18]
    out['scenarios_in_8_18'] = hits
    OUT.write_text(json.dumps(out, indent=1))
    print(json.dumps(out['rest'], indent=0))
    print('LSR', out['lsr'])
    print('n scenarios', len(sc), 'in +8..+18:', len(hits))
    for s in hits:
        print(f"{s['dV_mult']:7.2f} {s['dV_add']:7.2f} r={s['ratio']} | {s['measured']} | {s['reference']} | {s['o3_rest']} | {s['frame']}")
    print('small (|dV|<3):', len([s for s in sc if abs(s['dV_mult']) < 3]))
    print('Ren', out['ren_extrapolation']['blue_minus_red_all'], out['ren_extrapolation']['blue_minus_red_after_2018Nov'])
    for r in out['ren_extrapolation']['rows']:
        print(r)
    for r in out['doublet_fit_sim']:
        if r['start'] == 5199.35 or r['half_window'] >= 4:
            print(r)


if __name__ == '__main__':
    main()
