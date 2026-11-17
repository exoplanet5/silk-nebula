"""verify_altmodel_common: observational constraints, geometry and conversion helpers shared by the
verify_altmodel_* scripts.  Values are copied from the wave-1 notes / wave-2 derived files named in
SOURCES (they are *inputs* to the modelling; each is labelled measurement / derivation)."""
import numpy as np

PC = 3.0857e18
R_PH = 1e6 / (4 * np.pi)          # 1 Rayleigh in photons cm^-2 s^-1 sr^-1
H_PL, C_L, K_B, EV = 6.62607015e-27, 2.99792458e10, 1.380649e-16, 1.602176634e-12
M_H = 1.6726e-24
E5007 = H_PL * C_L / 5008.24e-8   # erg per [O III] 5007 photon (vacuum wavelength)

# ----------------------------------------------------------------------------------------------
# Observed spectrum of the filament itself (differential / ON-OFF measurements).  Energy ratios.
# ----------------------------------------------------------------------------------------------
OBS = {
    # [O III]5007 peak surface brightness (R). Wave-2: DESI south median 1.6 R (E strand 1.4-2.3, W strand to 3.9 R;
    # verify_desi.md), SDSS ON-OFF 1.39 R and per-plate peaks 1.6-2.5 R (verify_sdss.md), NSNS 2.6-3.6 R (wave 1).
    "I5007_R": dict(fid=2.0, lo=1.4, hi=4.0),
    # [O III]5007/Hbeta (energy). DESI wave-2 47 (>24 at 97.5 %, but blue null-line systematics +-0.01-0.02 of
    # [O III] and Hgamma/Hbeta = 1.5 suggest Hbeta may be underestimated); SDSS ~20 (>7 at 2 sigma); wave-1 DESI 20.
    # Adopted: robust lower bound 15, values 20 and 47 as bracketing estimates.
    "OIII_Hb": dict(fid=20.0, lo=15.0, hi=47.0),
    # [O II]3727/[O III]5007: DESI 1.37 (1.21-2.07), south 1.11; SDSS south 1.84 (1.46-2.37)
    "OII_OIII": dict(fid=1.37, lo=1.1, hi=1.9),
    # [O III]4363/5007: DESI wave-2 0.083 (0.059-0.091), south 0.092 (0.081-0.096), Hg coefficient ~0; SDSS 0.038+-0.035
    "O4363_5007": dict(fid=0.083, lo=0.059, hi=0.096),
    # He II 4686/[O III]5007 < 0.026 (97.5 %, DESI wave 2)
    "HeII_OIII_max": 0.026,
    # [Ne III]3869/[O III] = 0.099 (0.057-0.109); [S III]9069/[O III] 0.028; [O II]7320+30/[O III] 0.036 (0.033-0.042)
    "NeIII_OIII": dict(fid=0.099, lo=0.057, hi=0.109),
    "OII7325_OIII": dict(fid=0.036, lo=0.033, hi=0.042),
    # [O II]3729/3726 = 1.30 (1.24-1.33; subsets 1.17-1.37), SDSS 1.29 (1.07-1.59)
    "OII_doublet": dict(fid=1.30, lo=1.17, hi=1.37),
    # [N II]5755/[O III] = 0.014 (0.005-0.026) and [N II]6583/[O III] = 0.10 (DESI 'all'; red lines less reliable)
    "NII5755_OIII": 0.014, "NII6583_OIII": 0.10,
    # differential red-line ratios (energy): DESI south [N II]/Ha 1.08, [S II]/Ha 0.51; all 1.21, 0.75
    "NII_Ha": dict(fid=1.1, lo=1.0, hi=1.4), "SII_Ha": dict(fid=0.6, lo=0.45, hi=0.8),
    # GALEX FUV regression per R of [O III]5007 (verify_xrayfuv): 74+-8 (all), 92+-10 (south), 68 if [O III] = 1.35 x NSNS
    "FUV_CU_per_R": dict(fid=80.0, lo=68.0, hi=92.0), "NUV_CU_per_R": dict(fid=45.0, lo=36.0, hi=80.0),
    # line width: intrinsic FWHM <~33 km/s (DESI wave 1), V_LSR ~ -2..+2
    "FWHM_kms_max": 36.0,
}

GEOM = {
    "d_pc": dict(fid=280.0, lo=250.0, hi=320.0),   # Monogem-limb distance (context.md)
    "fwhm_deg": 0.28,                              # E-strand FWHM (geometry.md)
    "length_deg": 16.9,
    "R_shell_pc": 60.0,                            # Monogem radius 11.97 deg at 288 pc (Knies+2024)
}

# Monogem hot gas (Knies+2024, literature via context.md / imaging_mw.md)
HOT = {"kT_keV": 0.14, "T": 0.14 * 1.16045e7, "ne": dict(fid=5e-3, lo=3.1e-3, hi=7.2e-3), "age_yr": 1.1e5,
       "R_pc": 60.0}


def sheet_G(d_pc=280.0, fwhm_deg=0.28, R_c_pc=60.0):
    """limb-brightening factor of a sheet of thickness ~ projected FWHM, radius of curvature R_c:
    G = L_tangent / Delta with L = 2 sqrt(2 R_c Delta)."""
    delta = np.deg2rad(fwhm_deg) * d_pc
    return 2 * np.sqrt(2 * R_c_pc / delta), delta


def required_face_flux(I_R=3.0, G=10.0):
    """[O III]5007 power per unit area of the emitting layer (photons and erg cm^-2 s^-1) needed to
    give I_R Rayleigh when viewed with geometric enhancement G (edge-on path / thickness)."""
    phot = 4 * np.pi * I_R * R_PH / G
    return phot, phot * E5007


SOURCES = {
    "spectra_desi": "notes/spectra_desi.md; data/derived/spectra_desi_ratios.json",
    "verify_desi": "notes/verify_desi.md; data/derived/verify_desi_ratios.json (wave-2, final)",
    "spectra_sdss": "notes/spectra_sdss_lamost.md",
    "imaging_mw": "notes/imaging_mw.md; data/derived/imaging_mw_fuv_prediction.json",
    "verify_xrayfuv": "notes/verify_xrayfuv.md; data/derived/verify_xrayfuv_galex_measure.json (wave-2, final)",
    "verify_sdss": "notes/verify_sdss.md (wave-2, final)",
    "context": "notes/context.md (Knies+2024 Monogem; hot stars)",
}
