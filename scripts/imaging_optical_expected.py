"""Expected broadband / photographic surface brightness of the Silk Nebula emission lines.

For a line of intensity I (Rayleigh) at lambda0 seen through a photon-counting bandpass T(lambda):
    f_nu,eff = F_line * T(lambda0) * lambda0 / (c * Int T dlambda/lambda),   m_AB = -2.5 log10 f_nu,eff - 48.6
with F_line = I * 1e6/(4 pi) * h c / lambda0  [erg s-1 cm-2 sr-1] / 4.2545e10 arcsec2 sr-1.
Filter curves: PS1 g/r capture cross-sections (Tonry et al. 2012, ApJ 750, 99, Table 3, including 1.2 airmass,
https://ipp.ifa.hawaii.edu/ps1.filters/apj425122t3_mrt.txt); SDSS g/r 2010 response (speclite,
https://github.com/desihub/speclite). Photographic bands use top-hat approximations (see notes).
Output: data/derived/imaging_optical_expected_signal.json
"""
import sys
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import numpy as np
from imaging_optical_common import RAW, DER, save_json

H, C = 6.62607e-27, 2.99792458e10
SR2AS2 = 4.25452e10
LINES = {'OIII5007': 5006.84, 'OIII4959': 4958.91, 'Halpha': 6562.8, 'NII6583': 6583.45, 'SII6716': 6716.44,
         'SII6731': 6730.82, 'Hbeta': 4861.3}


def rayleigh_to_flux(I, lam):
    """erg s-1 cm-2 arcsec-2 for I Rayleigh at lam (Angstrom)."""
    return I * 1e6 / (4 * np.pi) * H * C / (lam * 1e-8) / SR2AS2


def load_curves():
    cur = {}
    t = np.loadtxt(RAW / 'filter_curves/apj425122t3_mrt.txt', skiprows=26, usecols=(0, 2, 3))
    cur['PS1_g'] = (t[:, 0] * 10.0, t[:, 1]); cur['PS1_r'] = (t[:, 0] * 10.0, t[:, 2])
    for b in 'gr':
        rows = [l.split() for l in open(RAW / f'filter_curves/sdss2010-{b}.ecsv') if l[0].isdigit()]
        a = np.array(rows, float); cur[f'SDSS_{b}'] = (a[:, 0], a[:, 1])
    # photographic top-hats (approximate; emulsion+filter): POSS-II J (IIIaJ+GG385) 3850-5400,
    # POSS-II F (IIIaF+RG610) 5900-6900, POSS-I E (103aE+plexiglass) 6200-6700, POSS-I O (103aO, no filter) 3300-5000
    for name, (a, b) in {'POSSII_J': (3850, 5400), 'POSSII_F': (5900, 6900), 'POSSI_E': (6200, 6700),
                         'POSSI_O': (3300, 5000)}.items():
        lam = np.arange(a - 50, b + 51, 1.0); T = ((lam >= a) & (lam <= b)).astype(float)
        cur[name] = (lam, T)
    return cur


def line_mag(I, lam0, curve):
    lam, T = curve
    norm = np.trapezoid(T / lam, lam)
    T0 = np.interp(lam0, lam, T, left=0, right=0)
    if T0 <= 0:
        return np.inf, 0.0
    F = rayleigh_to_flux(I, lam0)
    fnu = F * T0 * (lam0 * 1e-8) / (C * norm * 1e-8)   # norm has units of dimensionless (A/A)
    # careful with units: lam0 in cm, c in cm/s -> lam0/c in s; norm dimensionless
    fnu = F * T0 * (lam0 * 1e-8) / C / norm
    return -2.5 * np.log10(fnu) - 48.6, T0 / T.max()


SKY = {  # typical dark-sky surface brightness, mag/arcsec^2 (literature; used only for contrast estimates)
    'PS1_g': 21.9, 'PS1_r': 20.9,        # Haleakala (Chambers et al. 2016 quote ~21.9 g, 20.9 r at zenith dark)
    'SDSS_g': 22.0, 'SDSS_r': 21.0,      # Apache Point typical
    'POSSII_J': 22.3, 'POSSII_F': 21.0,  # Palomar B_J / R_F 1980s-90s (Reid et al. 1991 era plates)
    'POSSI_E': 20.9, 'POSSI_O': 22.4,
}


def main():
    cur = load_curves()
    out = {'method': __doc__.strip().split('\n')[2:5], 'sky_assumed_mag_arcsec2': SKY, 'cases': {}}
    cases = {
        'OIII_12R_paper': [('OIII5007', 12.0), ('OIII4959', 4.0)],
        'OIII_3R_NSNS_peak': [('OIII5007', 3.0), ('OIII4959', 1.0)],
        'OIII_2R_NSNS_mean': [('OIII5007', 2.0), ('OIII4959', 0.667)],
        'Halpha_patch_1.3R_plus_SII_1R': [('Halpha', 1.3), ('SII6716', 0.5), ('SII6731', 0.5), ('NII6583', 0.8)],
        'Halpha_10R': [('Halpha', 10.0)],
    }
    for cname, comps in cases.items():
        res = {}
        for band, curve in cur.items():
            fnu_tot = 0.0; detail = {}
            for line, I in comps:
                m, trel = line_mag(I, LINES[line], curve)
                if np.isfinite(m):
                    fnu_tot += 10 ** (-0.4 * (m + 48.6))
                detail[line] = dict(I_R=I, T_rel=round(float(trel), 3))
            if fnu_tot <= 0:
                continue
            m = -2.5 * np.log10(fnu_tot) - 48.6
            res[band] = dict(mag_arcsec2=round(float(m), 2),
                             fraction_of_sky=float(10 ** (-0.4 * (m - SKY[band]))), lines=detail)
        out['cases'][cname] = res
    # handy number: 1 R of [OIII]5007 in erg/s/cm2/arcsec2
    out['flux_1R_OIII5007_cgs_arcsec2'] = rayleigh_to_flux(1, LINES['OIII5007'])
    out['flux_1R_Halpha_cgs_arcsec2'] = rayleigh_to_flux(1, LINES['Halpha'])
    # paper Table 3: [OIII]-band mag 23.6-23.9 in a 6 nm FWHM filter -> Rayleigh (AB, top-hat 60 A)
    conv = {}
    for m in (23.60, 23.87, 23.85):
        fnu = 10 ** (-0.4 * (m + 48.6)); lam = LINES['OIII5007']
        F = fnu * C / (lam * 1e-8) ** 2 * 60e-8   # erg/s/cm2/arcsec2 within 60 A (fnu*dnu)
        conv[str(m)] = F / rayleigh_to_flux(1, lam)
    out['paper_table3_mag_to_R_tophat60A'] = conv
    save_json(out, DER / 'imaging_optical_expected_signal.json')
    for c, r in out['cases'].items():
        print(c, {b: (v['mag_arcsec2'], f"{100 * v['fraction_of_sky']:.2f}%") for b, v in r.items()})
    print('paper mag->R', conv, '1R OIII cgs', out['flux_1R_OIII5007_cgs_arcsec2'])


if __name__ == '__main__':
    main()
