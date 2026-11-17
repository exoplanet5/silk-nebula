"""spectra_desi_lines.py -- line list, DESI wavelength grid and saved windows for the DESI residual analysis.

DESI DR1 healpix coadds (as served by SPARCL) are on one fixed grid: 3600.0 + 0.8*i A, i = 0..7780 (verified on
the pilot download). Wavelengths are VACUUM and barycentric: desi_proc passes --barycentric-correction to
desi_extract_spectra (desispec scripts/proc.py, scripts/extract.py on github.com/desihub/desispec, read 2026-11-16).
Barycentric minus heliocentric velocity is < 0.015 km/s, so DESI velocities compare directly with LAMOST V_helio.

Rest wavelengths: standard air values (NIST / Osterbrock & Ferland) converted to vacuum with the Morton (2000)
IAU formula below; values quoted in the task brief are reproduced to <= 0.02 A.
"""
import numpy as np

C_KMS = 299792.458
WAVE = 3600.0 + 0.8 * np.arange(7781)
NPIX = WAVE.size


def air_to_vac(lam_air):
    lam = np.atleast_1d(np.asarray(lam_air, float))
    vac = lam.copy()
    for _ in range(3):                       # iterate: formula is in terms of vacuum wavenumber
        s2 = (1e4 / vac) ** 2
        n = 1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2)
        vac = lam * n
    return vac if vac.size > 1 else float(vac[0])


# name: (air wavelength A, kind)  kind: 'neb' nebular line, 'sky' airglow check line
_AIR = {
    'OII3726': (3726.03, 'neb'), 'OII3729': (3728.82, 'neb'), 'NeIII3869': (3868.76, 'neb'),
    'Hd': (4101.74, 'neb'), 'Hg': (4340.47, 'neb'), 'OIII4363': (4363.21, 'neb'),
    'HeII4686': (4685.71, 'neb'), 'Hb': (4861.33, 'neb'), 'OIII4959': (4958.91, 'neb'),
    'OIII5007': (5006.84, 'neb'), 'NI5198': (5197.90, 'neb'), 'NI5200': (5200.26, 'neb'),
    'NII5755': (5754.64, 'neb'), 'HeI5876': (5875.62, 'neb'), 'OI6300': (6300.30, 'neb'),
    'OI6364': (6363.78, 'neb'), 'NII6548': (6548.05, 'neb'), 'Ha': (6562.80, 'neb'),
    'NII6583': (6583.45, 'neb'), 'HeI6678': (6678.15, 'neb'), 'SII6716': (6716.44, 'neb'),
    'SII6731': (6730.82, 'neb'), 'HeI7065': (7065.19, 'neb'), 'ArIII7136': (7135.80, 'neb'),
    'OII7320': (7319.99, 'neb'), 'OII7330': (7330.73, 'neb'), 'SIII9069': (9068.6, 'neb'),
    'SIII9532': (9530.6, 'neb'),
    'sky5577': (5577.34, 'sky'), 'sky6300': (6300.30, 'sky'), 'sky5199': (5199.0, 'sky'),
}
LINES = {k: air_to_vac(v[0]) for k, v in _AIR.items()}
KIND = {k: v[1] for k, v in _AIR.items()}

# Windows saved to disk: +-40 A around each line, merged when overlapping.
HALF = 40.0


def _windows():
    c = sorted(set(round(v, 2) for v in LINES.values()))
    segs = []
    for w in c:
        lo, hi = w - HALF, w + HALF
        if segs and lo <= segs[-1][1]:
            segs[-1][1] = max(segs[-1][1], hi)
        else:
            segs.append([lo, hi])
    return [(max(lo, WAVE[0]), min(hi, WAVE[-1])) for lo, hi in segs]


WINDOWS = _windows()
WIN_IDX = np.concatenate([np.where((WAVE >= lo) & (WAVE <= hi))[0] for lo, hi in WINDOWS])
WIN_WAVE = WAVE[WIN_IDX]

# Lines of extragalactic objects (vacuum rest) used to flag self-contamination of the z~0 windows
GAL_LINES = np.array([1215.67, 1240.14, 1397.6, 1549.06, 1640.4, 1908.73, 2799.12, 3346.8, 3426.9, 3728.48,
                      3869.86, 3934.78, 3969.59, 4102.89, 4305.6, 4341.69, 4862.68, 4960.30, 5008.24, 5176.7,
                      5895.6, 6302.05, 6549.86, 6564.61, 6585.27, 6718.29, 6732.67])

if __name__ == '__main__':
    for k, v in LINES.items():
        print(f'{k:10s} {v:9.2f}')
    print(len(WINDOWS), 'windows,', WIN_IDX.size, 'pixels')
    for lo, hi in WINDOWS:
        print(f'  {lo:7.1f}-{hi:7.1f}')
