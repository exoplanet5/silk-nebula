"""verify_sdss stage 1: per-fibre geometry, plate-spectrograph groups, absolute sky fits per group.

Per group (release, plate, mjd, spectrograph), from the clipped mean TOTAL sky (flux+sky) of its sky fibres:
  * [O I]5577 and [O I]6300 Gaussian fits: intensity (R, as calibrated) and centroid velocity on the pipeline grid
  * continuum surface brightness at 5040-5100 A and synthetic g, r surface brightness (AB mag/arcsec^2, as calibrated)
  * absolute red-window model: CALSPEC solar model (R=5000, vacuum) broadened to the SDSS LSF x linear + additive
    constant + Gaussians for Halpha, [N II], [S II] and OH lines (positions found empirically)  -> I(Halpha)_abs
  * astropy heliocentric correction for groups with a TAI timestamp (SDSS-I/II, SDSS-V) to test the HELIO_RV sign
Writes data/derived/verify_sdss_fibres.csv, verify_sdss_groups.csv and verify_sdss_templates.npz (per-group
continuum and OH templates used by the differential fits).
"""
import sys, csv
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_sdss_lib import (ROOT, VAC, C, load_windows, load_strands, strand_distance, in_polygon, flam_to_R_per_A,
                             fit_lines, gauss_unit, linfit)

APO = EarthLocation.from_geodetic(lon=-105.8202 * u.deg, lat=32.7802 * u.deg, height=2788 * u.m)


def sun_template():
    with fits.open(ROOT / 'data/raw/calspec_sun/sun_mod_001.fits') as h:
        w = np.asarray(h[1].data['WAVELENGTH'], float); f = np.asarray(h[1].data['FLUX'], float)
    return w, f


SUNW, SUNF = sun_template()


def sun_on(lam, sig_A, v=0.0):
    """Solar model broadened by an extra Gaussian of sigma sig_A (A, array over lam) and sampled at lam."""
    out = np.empty_like(lam)
    s_model = lam / 5000.0 / 2.3548
    se = np.sqrt(np.clip(np.asarray(sig_A, float) ** 2 - s_model ** 2, 0.05 ** 2, None))
    cut = (SUNW > lam.min() - 25) & (SUNW < lam.max() + 25)
    lw = SUNW[cut] * (1 + v / C)
    fw = SUNF[cut]
    for i, (l0, s) in enumerate(zip(lam, se)):
        k = np.exp(-0.5 * ((lw - l0) / s) ** 2)
        out[i] = np.sum(k * fw) / np.sum(k)
    return out / np.median(out)


def clipped_mean(X, V, nsig=4.0, axis=0):
    """Per-pixel clipped mean over fibres; X (nf, npix) values, V variances (inf where bad)."""
    good = np.isfinite(X) & np.isfinite(V)
    Xm = np.where(good, X, np.nan)
    med = np.nanmedian(Xm, axis=0)
    mad = 1.4826 * np.nanmedian(np.abs(Xm - med), axis=0)
    ok = good & (np.abs(Xm - med) < nsig * np.maximum(mad, 1e-6) + 1e-12)
    n = ok.sum(0)
    mean = np.where(n > 0, np.nansum(np.where(ok, X, 0), 0) / np.maximum(n, 1), np.nan)
    var = np.where(n > 0, np.nansum(np.where(ok, V, 0), 0) / np.maximum(n, 1) ** 2, np.inf)
    return mean, var, n


def main():
    D = load_windows()
    lam = D['lam']
    nf = D['tot'].shape[0]
    S, patch, tr = load_strands()
    ra, dec = D['m_ra'], D['m_dec']
    dist = {}
    for k in ('E', 'W', 'W2'):
        d, end, dnear = strand_distance(ra, dec, S[k])
        dist[k] = np.where(end, np.sign(d) * np.maximum(np.abs(d), 0.0) + 0.0, d)
        dist[k + '_end'] = end
        dist[k + '_decnear'] = dnear
    absd = np.vstack([np.abs(dist['E']), np.abs(dist['W']), np.abs(dist['W2'])])
    near = np.array(['E', 'W', 'W2'])[np.argmin(absd, 0)]
    dmin = absd.min(0)
    # beyond the strand ends: treat end-point distances as valid only if > 0.3 deg (keeps them as OFF candidates)
    inpatch = in_polygon(ra, dec, patch)
    outer = np.array(tr['halpha_patch_info']['outer_polygon_30pct'], float)
    inpatch30 = in_polygon(ra, dec, outer)
    # --------------------------------------------------- fibre table
    keys = [(r, p, m, s) for r, p, m, s in zip(D['m_release'], D['m_plate'], D['m_mjd'], D['m_sp'])]
    uk = sorted(set(keys))
    gid = np.array([uk.index(k) for k in keys])
    with open(ROOT / 'data/derived/verify_sdss_fibres.csv', 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['i', 'gid', 'release', 'plate', 'mjd', 'mjdfirst', 'sp', 'fiber', 'ra', 'dec', 'instrument', 'diam',
                    'dE', 'dW', 'dW2', 'E_end', 'W_end', 'W2_end', 'nearest', 'dmin', 'in_patch', 'in_patch30'])
        for i in range(nf):
            w.writerow([i, gid[i], D['m_release'][i], D['m_plate'][i], D['m_mjd'][i], D['m_mjdfirst'][i], D['m_sp'][i],
                        D['m_fiber'][i], f'{ra[i]:.5f}', f'{dec[i]:.5f}', D['m_instrument'][i], D['m_diam'][i],
                        f"{dist['E'][i]:.4f}", f"{dist['W'][i]:.4f}", f"{dist['W2'][i]:.4f}", int(dist['E_end'][i]),
                        int(dist['W_end'][i]), int(dist['W2_end'][i]), near[i], f'{dmin[i]:.4f}', int(inpatch[i]),
                        int(inpatch30[i])])
    # --------------------------------------------------- per-group absolute fits
    wd = D['wdisp']
    TOT = D['tot'].astype(float); IV = D['ivar'].astype(float)
    sig_pix = lambda i_fib: None
    rows = []
    tmpl_cont = np.full((len(uk), lam.size), np.nan, np.float32)
    tmpl_oh = np.zeros((len(uk), lam.size), np.float32)
    # empirical OH line list near Halpha / [S II] from a high-S/N global stack (found below on first pass)
    oh_lines = None
    for g, key in enumerate(uk):
        sel = np.where(gid == g)[0]
        diam = float(D['m_diam'][sel[0]])
        X = flam_to_R_per_A(TOT[sel], lam, diam)
        V = np.where(IV[sel] > 0, flam_to_R_per_A(1.0, lam, diam) ** 2 / np.where(IV[sel] > 0, IV[sel], 1), np.inf)
        M, MV, Mn = clipped_mean(X, V)
        wdm = np.nanmedian(wd[sel], 0)
        sigA = wdm * lam * np.log(10) * 1e-4           # LSF sigma in A
        r = dict(gid=g, release=key[0], plate=key[1], mjd=key[2], sp=key[3], n_fib=len(sel),
                 mjdfirst=int(D['m_mjdfirst'][sel].min()),
                 ra=float(np.mean(ra[sel])), dec=float(np.mean(dec[sel])), instrument=D['m_instrument'][sel[0]],
                 diam=diam, helio_rv=float(np.median(D['m_helio_rv'][sel])), tai=float(np.nanmedian(D['m_tai'][sel])),
                 airmass=float(np.nanmedian(D['m_airmass'][sel])), seeing50=float(np.nanmedian(D['m_seeing50'][sel])),
                 exptime=float(np.nanmedian(D['m_exptime'][sel])))
        # synthetic g/r SB (AB mag/arcsec^2, as calibrated): f_nu(Jy) per fibre / area
        area = np.pi * (diam / 2) ** 2
        for b in 'gr':
            fn = np.nanmedian(D[f'm_{b}_tot'][sel])
            r[f'sb_{b}_cal'] = float(-2.5 * np.log10(fn / area / 3631.0)) if fn > 0 else np.nan
        # continuum 5040-5100: R/A -> f_lambda -> AB
        m = (lam > 5040) & (lam < 5100) & np.isfinite(M)
        cont = np.nanmedian(M[m])                         # R/A (photon) as calibrated
        # convert R/A to f_nu per arcsec^2: photons/cm2/s/A/arcsec2 = R/A * 1e6/(4pi) * ARCSEC2_SR
        ph = cont * 1e6 / (4 * np.pi) * (np.pi / 648000.0) ** 2
        flam = ph * 6.62607015e-27 * 2.99792458e18 / 5070.0
        fnu = flam * 5070.0 ** 2 / 2.99792458e18
        r['cont5070_RperA_cal'] = float(cont)
        r['sb5070_AB_cal'] = float(-2.5 * np.log10(fnu / 3631e-23)) if fnu > 0 else np.nan
        try:
            # [O I] 5577 and 6300: free centroid
            for nm, lo, hi in (('OI5577', 5555, 5600), ('OI6300', 6285, 6320)):
                mm = (lam > lo) & (lam < hi) & np.isfinite(M) & np.isfinite(MV)
                best = None
                for v in np.arange(-60, 60.1, 1.0):
                    f = fit_lines(lam[mm], M[mm], 1 / MV[mm], [nm], v, sigA[mm] if False else (lambda x, s=np.nanmedian(sigA[mm]): s), npoly=2)
                    if best is None or f['_chi2'] < best[1]['_chi2']:
                        best = (v, f)
                # refine with a parabola on chi2
                vs = np.arange(best[0] - 3, best[0] + 3.01, 0.25)
                ch = [fit_lines(lam[mm], M[mm], 1 / MV[mm], [nm], v, (lambda x, s=np.nanmedian(sigA[mm]): s), npoly=2)['_chi2'] for v in vs]
                pc = np.polyfit(vs, ch, 2)
                vbest = -pc[1] / (2 * pc[0]) if pc[0] > 0 else best[0]
                f = fit_lines(lam[mm], M[mm], 1 / MV[mm], [nm], vbest, (lambda x, s=np.nanmedian(sigA[mm]): s), npoly=2)
                # width-free check is skipped: SDSS LSF from wdisp
                r[f'v_{nm}'] = float(vbest)
                r[f'ev_{nm}'] = float(np.sqrt(2 / pc[0])) if pc[0] > 0 else np.nan
                r[f'I_{nm}_cal'] = f[nm][0]
                r[f'chi2_{nm}'] = f['_chi2'] / max(f['_dof'], 1)
            # [N I] 5198/5200 blend: two Gaussians at fixed separation with a common velocity scan, ratio free
            mm = (lam > 5180) & (lam < 5222) & np.isfinite(M) & np.isfinite(MV)
            sN = float(np.nanmedian(sigA[mm]))
            sunN = sun_on(lam[mm], sigA[mm])
            best = None
            for v in np.arange(-60, 60.1, 1.0):
                f = fit_lines(lam[mm], M[mm], 1 / MV[mm], ['NI5198', 'NI5200'], v, sN, npoly=1, extra={'sun': sunN})
                if best is None or f['_chi2'] < best[1]['_chi2']:
                    best = (v, f)
            r['v_NI_blend'] = float(best[0]); r['I_NI5198_cal'] = best[1]['NI5198'][0]; r['I_NI5200_cal'] = best[1]['NI5200'][0]
            # red absolute model
            mm = (lam > 6505) & (lam < 6790) & np.isfinite(M) & np.isfinite(MV)
            lr = lam[mm]; yr = M[mm]; wr = 1 / MV[mm]
            sunr = sun_on(lr, sigA[mm], v=0.0)
            if oh_lines is None:
                # topocentric vacuum OH positions identified empirically in the global mean sky (see notes)
                oh_lines = np.array([6524.0, 6535.0, 6545.9, 6555.4, 6579.1, 6598.8, 6606.2, 6636.5, 6670.0, 6701.8, 6709.5, 6746.0, 6759.0, 6772.0])
            # OH lines are topocentric: on the heliocentric grid at lam0(1 + v_corr/c); v_corr ~ -helio_rv (tested below)
            vtop = -r['helio_rv']
            sR = lambda x, s=float(np.nanmedian(sigA[mm])): s
            extra = {'sun': sunr, 'sun_x': sunr * (lr - 6650) / 100.0}
            ohcols = {f'oh{j}': gauss_unit(lr, l0 * (1 + vtop / C), sR(l0)) for j, l0 in enumerate(oh_lines)}
            extra.update(ohcols)
            f = fit_lines(lr, yr, wr, ['Ha', 'NII6548', 'NII6583', 'SII6716', 'SII6731'], vtop + 0.0, sR, npoly=1,
                          ties={'NII6548': ('NII6583', 1 / 3.06)}, extra=extra)
            # the Galactic lines (V_helio ~ +12) and geocoronal Halpha (topocentric) are unresolved; use the topocentric
            # frame for all (difference < 1 LSF sigma).
            r['I_Ha_abs_cal'] = f['Ha'][0]; r['eI_Ha_abs_cal'] = f['Ha'][1]
            r['I_NII6583_abs_cal'] = f['NII6583'][0]; r['I_SII_abs_cal'] = f['SII6716'][0] + f['SII6731'][0]
            r['I_OH6555_cal'] = f['oh2'][0]
            r['red_chi2nu'] = f['_chi2'] / max(f['_dof'], 1)
            # templates for the differential fits (on the full lam grid, red window only)
            names = f['_names']
            cvec = np.array([f[n][0] for n in names])
            A_cont = np.zeros_like(lr); A_oh = np.zeros_like(lr)
            for n, c_ in zip(names, cvec):
                if n.startswith('poly') or n.startswith('sun'):
                    col = (np.ones_like(lr) if n == 'poly0' else extra[n]) if n in ('poly0', 'sun', 'sun_x') else None
                    if col is not None:
                        A_cont += c_ * col
                if n.startswith('oh'):
                    A_oh += c_ * extra[n]
            idx = np.where(mm)[0]
            tmpl_cont[g, idx] = A_cont; tmpl_oh[g, idx] = A_oh
            # blue continuum template (solar x linear + const) for 4300-5110 and 3690-3770
            for lo, hi in ((3690, 3770), (4300, 4410), (4640, 4730), (4820, 5110)):
                mb = (lam > lo) & (lam < hi) & np.isfinite(M) & np.isfinite(MV)
                if mb.sum() < 20:
                    continue
                lb = lam[mb]
                sb = sun_on(lb, sigA[mb])
                ex = {'sun': sb, 'sun_x': sb * (lb - lb.mean()) / 100.0}
                lines = [l for l in ('OII3726', 'OII3729', 'Hgam', 'Hg4358', 'OIII4363', 'HeII4686', 'Hb', 'OIII4959', 'OIII5007')
                         if lo < VAC[l] < hi]
                fb = fit_lines(lb, M[mb], 1 / MV[mb], lines, 0.0, (lambda x, s=float(np.nanmedian(sigA[mb])): s), npoly=1, extra=ex)
                Acb = fb['poly0'][0] + fb['sun'][0] * sb + fb['sun_x'][0] * ex['sun_x']
                tmpl_cont[g, np.where(mb)[0]] = Acb
                if lo == 4820:
                    r['I_Hb_abs_cal'] = fb['Hb'][0]; r['I_OIII5007_abs_cal'] = fb['OIII5007'][0]
                    # solar-spectrum (zodiacal + starlight) component vs featureless additive component, R/A at 4965 A
                    r['sun_coef_blue'] = fb['sun'][0]; r['sun_coef_blue_err'] = fb['sun'][1]
                    r['add_coef_blue'] = fb['poly0'][0]
                if lo == 4300:
                    r['I_Hg4358_abs_cal'] = fb['Hg4358'][0]
        except Exception as ex:
            r["fit_error"] = repr(ex)[:80]
            print("fit error", key, repr(ex)[:80], flush=True)
        # astropy heliocentric correction (where TAI known)
        if np.isfinite(r['tai']) and r['tai'] > 0:
            t = Time(r['tai'] / 86400.0, format='mjd', scale='tai')
            sc = SkyCoord(r['ra'] * u.deg, r['dec'] * u.deg)
            r['vcorr_astropy'] = float(sc.radial_velocity_correction(kind='heliocentric', obstime=t, location=APO).to(u.km / u.s).value)
        else:
            r['vcorr_astropy'] = np.nan
        rows.append(r)
        if g % 25 == 0:
            print(g, key, r['n_fib'], f"5577 v={r['v_OI5577']:+.1f} I={r['I_OI5577_cal']:.0f}R  helio_rv={r['helio_rv']:+.1f} "
                  f"vc={r['vcorr_astropy']:+.1f} sb_g={r['sb_g_cal']:.2f} Ha_abs={r['I_Ha_abs_cal']:.2f}", flush=True)
    keys_out = list(rows[0].keys())
    for r in rows:
        for k in r:
            if k not in keys_out:
                keys_out.append(k)
    with open(ROOT / 'data/derived/verify_sdss_groups.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys_out)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, np.nan) for k in keys_out})
    np.savez_compressed(ROOT / 'data/derived/verify_sdss_templates.npz', cont=tmpl_cont, oh=tmpl_oh, gid_of_fibre=gid)
    print('groups', len(rows))


if __name__ == '__main__':
    main()
