"""DSS1/DSS2 photographic precovery search along the Silk Nebula.

Input : data/raw/skyview/<band>/dss_<band>_d<dec>_<E|W>.fits (SkyView, 6"/pix TAN, Clip sampler: each pixel is the SUM
        of ~36 native 1" DSS2 pixel values (density units), ~12.5 native 1.7" DSS1 pixels).
        data/raw/gaia_dr3_datalab/*.csv (G<17.5, masking); data/raw/ps1_dr2/catalog/*.csv (faint stars, calibration)
Steps per tile:
  1. mask Gaia stars with r(G) = 8" + 110" 10^(-0.2 (G-6)) (cap 300"), plus 3-sigma positive outliers vs a 1.5' median
     background (dilated by one pixel).
  2. calibrate density response: for isolated PS1 stars with 18.5 < m < 20.5 (g for blue plates, r for red),
     density 'volume' V = sum(D - D_sky) within r = 15"; fit V = k * 10^(-0.4 m). In the small-signal (linear) regime a
     diffuse density excess dD per 6" pixel then corresponds to mu = -2.5 log10(dD / (k * 36 arcsec^2)) mag/arcsec^2
     (AB, PS1 system). Faint-star peaks are not fully in the linear regime -> ~30% systematic (see notes).
  3. median-bin unmasked pixels into 1' blocks (10x10).
Per tile row (E+W tiles, ridge +-2 deg): cross-filament profile in 2' bins vs distance from the NSNS [OIII] ridge,
background = quadratic fit to 0.7 <= |x| <= 2.0 deg; template fit dD(x) = A * O3_NSNS(x) + poly2(x) gives the
[OIII]-template amplitude A (density per R of NSNS [OIII]) with a bootstrap-free error from residual scatter.
Outputs: data/derived/imaging_optical_dss_profiles.json, imaging_optical_dss_binned.npz; figures via --fig.
"""
import sys
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import glob, json
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import median_filter, binary_dilation, gaussian_filter
from imaging_optical_common import *
from imaging_optical_nsns import ridge_spine

BANDS = {'DSS2_Blue': dict(plate='POSS-II J (IIIaJ+GG385)', cal='g', native_arcsec=1.0),
         'DSS2_Red': dict(plate='POSS-II F (IIIaF+RG610)', cal='r', native_arcsec=1.0),
         'DSS1_Blue': dict(plate='POSS-I O (103aO)', cal='g', native_arcsec=1.7),
         'DSS1_Red': dict(plate='POSS-I E (103aE+plexiglass)', cal='r', native_arcsec=1.7)}
PIX = 6.0  # arcsec


def load_gaia():
    rows = []
    for f in sorted(glob.glob(str(RAW / 'gaia_dr3_datalab/gaia_g175_dec*.csv'))):
        a = np.genfromtxt(f, delimiter=',', names=True)
        rows.append(np.c_[a['ra'], a['dec'], a['phot_g_mean_mag']])
    return np.vstack(rows)


def load_ps1cat(name):
    a = np.genfromtxt(RAW / f'ps1_dr2/catalog/ps1cat_{name}.csv', delimiter=',', names=True, skip_header=0)
    return a


def star_mask(data, w, gaia):
    ny, nx = data.shape
    x, y = w.wcs_world2pix(gaia[:, 0], gaia[:, 1], 0)
    G = gaia[:, 2]
    ok = (x > -60) & (x < nx + 60) & (y > -60) & (y < ny + 60)
    m = np.zeros(data.shape, bool)
    yy, xx = np.ogrid[0:ny, 0:nx]
    for xi, yi, g in zip(x[ok], y[ok], G[ok]):
        r = min(8 + 110 * 10 ** (-0.2 * (g - 6)), 300) / PIX
        r = max(r, 1.2)
        x0, x1 = int(max(0, xi - r - 1)), int(min(nx, xi + r + 2)); y0, y1 = int(max(0, yi - r - 1)), int(min(ny, yi + r + 2))
        if x1 <= x0 or y1 <= y0:
            continue
        sub = (xx[:, x0:x1] - xi) ** 2 + (yy[y0:y1, :] - yi) ** 2 <= r * r
        m[y0:y1, x0:x1] |= sub
    # outliers
    d = np.where(m, np.nan, data)
    bg = median_filter(np.nan_to_num(d, nan=np.nanmedian(d)), size=15)
    res = data - bg
    s = 1.4826 * np.nanmedian(np.abs(res[~m] - np.nanmedian(res[~m])))
    out = (res > 3 * s) | (res < -6 * s)
    m |= binary_dilation(out, iterations=1)
    m |= ~np.isfinite(data) | (data == 0)
    return m, s


def calibrate(data, w, cat, band_letter, mask):
    """Return k = density-volume per unit flux (10^-0.4 m), its scatter, N used."""
    col = f'{band_letter}MeanPSFMag'
    mag = cat[col]; ra = cat['raMean']; de = cat['decMean']
    good = np.isfinite(mag) & (mag > 18.5) & (mag < 20.5) & (mag > 0)
    x, y = w.wcs_world2pix(ra, de, 0)
    ny, nx = data.shape
    good &= (x > 8) & (x < nx - 8) & (y > 8) & (y < ny - 8)
    # isolation: no other PS1 star (any mag in list) within 20" and pixel not in Gaia-bright mask
    from scipy.spatial import cKDTree
    t = cKDTree(np.c_[ra * np.cos(np.radians(de)), de])
    nn = t.query_ball_point(np.c_[ra * np.cos(np.radians(de)), de], r=20 / 3600)
    iso = np.array([len(n) == 1 for n in nn])
    good &= iso
    V = []; F = []
    yy, xx = np.mgrid[-7:8, -7:8]
    rr = np.hypot(xx, yy)
    ap = rr <= 15 / PIX; an = (rr >= 25 / PIX) & (rr <= 42 / PIX)
    for xi, yi, mi in zip(x[good], y[good], mag[good]):
        ix, iy = int(round(xi)), int(round(yi))
        cut = data[iy - 7:iy + 8, ix - 7:ix + 8]
        mk = mask[iy - 7:iy + 8, ix - 7:ix + 8]
        if cut.shape != (15, 15):
            continue
        # gaia bright-star mask must not touch the annulus heavily
        if mk[an].mean() > 0.5:
            continue
        skyv = np.median(cut[an & ~mk]) if (an & ~mk).sum() > 10 else np.median(cut[an])
        V.append(np.sum(cut[ap] - skyv)); F.append(10 ** (-0.4 * mi))
    V = np.array(V); F = np.array(F)
    if len(V) < 20:
        return np.nan, np.nan, len(V)
    ratio = V / F
    k = np.median(ratio)
    sc = 1.4826 * np.median(np.abs(ratio - k)) / np.sqrt(len(V))
    return float(k), float(sc), int(len(V))


def block_median(data, mask, b=10, minfrac=0.3):
    ny, nx = data.shape
    ny2, nx2 = ny // b, nx // b
    d = np.where(mask, np.nan, data)[:ny2 * b, :nx2 * b].reshape(ny2, b, nx2, b).transpose(0, 2, 1, 3).reshape(ny2, nx2, b * b)
    frac = np.isfinite(d).mean(-1)
    med = np.nanmedian(d, -1)
    med[frac < minfrac] = np.nan
    return med, frac


def process_tile(path, gaia, band):
    with fits.open(path) as h:
        data = h[0].data.astype(float); w = WCS(h[0].header)
    name = '_'.join(Path(path).stem.split('_')[-2:])
    m, s = star_mask(data, w, gaia)
    cat = load_ps1cat(name)
    k, ksc, nk = calibrate(data, w, cat, BANDS[band]['cal'], m)
    med, frac = block_median(data, m, 10)
    wb = w[::10, ::10]  # 1' WCS (pixel centres shifted; fine at 1' level)
    # shift CRPIX for block centres: block j centre = 10j+4.5 in original pixel coords
    wb = WCS(w.to_header()); wb.wcs.crpix = [(w.wcs.crpix[0] - 0.5) / 10 + 0.5 - 0.45, (w.wcs.crpix[1] - 0.5) / 10 + 0.5 - 0.45]
    wb.wcs.cdelt = w.wcs.cdelt * 10
    sky_D = float(np.nanmedian(med))
    plates = json.loads(Path(str(path).replace('.fits', '.json')).read_text()).get('plates', [])
    return dict(name=name, binned=med, frac=frac, wcs=wb, k=k, k_err=ksc, n_cal=nk, sky_D=sky_D, pix_noise=s,
                masked_frac=float(m.mean()), plates=plates)


def row_profile(tiles_row, sra, sdec, o3maps):
    """Combine the E and W binned tiles of one row; cross-filament profile with background poly + NSNS template fit."""
    xs, vals, decs, ras = [], [], [], []
    for t in tiles_row:
        b = t['binned']
        yy, xx = np.mgrid[0:b.shape[0], 0:b.shape[1]]
        ra, dec = t['wcs'].wcs_pix2world(xx, yy, 0)
        ok = np.isfinite(b)
        # convert to calibrated surface-brightness units relative to the row: dD / (k*36) = flux per arcsec^2
        vals.append((b[ok] - t['sky_D']) / (t['k'] * 36.0) if np.isfinite(t['k']) else b[ok] * np.nan)
        ras.append(ra[ok]); decs.append(dec[ok])
    ra = np.concatenate(ras); dec = np.concatenate(decs); v = np.concatenate(vals)
    across, along, beyond = along_across(ra, dec, sra, sdec)
    # NSNS [OIII] at the same positions (nearest pixel of the 1' NSNS map)
    zw = WCS(str(o3maps['wcs_header']))
    px, py = zw.wcs_world2pix(ra, dec, 0)
    ix = np.clip(np.round(px).astype(int), 0, o3maps['oiii'].shape[1] - 1); iy = np.clip(np.round(py).astype(int), 0, o3maps['oiii'].shape[0] - 1)
    o3 = gaussian_filter(np.nan_to_num(o3maps['oiii']), 1.0)[iy, ix]
    edges = np.arange(-2.0, 2.0001, 2 / 60); cen = 0.5 * (edges[1:] + edges[:-1])
    prof = np.full(len(cen), np.nan); po3 = np.full(len(cen), np.nan); nn = np.zeros(len(cen), int)
    sel = ~beyond & np.isfinite(v)
    for i in range(len(cen)):
        s = sel & (across >= edges[i]) & (across < edges[i + 1])
        nn[i] = s.sum()
        if s.sum() >= 5:
            prof[i] = np.median(v[s]); po3[i] = np.median(o3[s])
    ok = np.isfinite(prof) & np.isfinite(po3)
    bgm = ok & (np.abs(cen) >= 0.7)
    out = dict(across_deg=cen.tolist(), profile=prof.tolist(), nsns_o3_profile_R=po3.tolist(), n_per_bin=nn.tolist())
    if bgm.sum() < 10:
        return out
    c = np.polyfit(cen[bgm], prof[bgm], 2)
    res = prof - np.polyval(c, cen)
    out['profile_bgsub'] = res.tolist()
    sig_bin = 1.4826 * np.nanmedian(np.abs(res[bgm] - np.nanmedian(res[bgm])))
    # template fit: res = A*(o3 - o3_bgpoly) ; o3 background removed the same way
    c3 = np.polyfit(cen[bgm], po3[bgm], 2); t3 = po3 - np.polyval(c3, cen)
    X = np.c_[t3[ok], np.ones(ok.sum()), cen[ok], cen[ok] ** 2]
    coef, *_ = np.linalg.lstsq(X, prof[ok], rcond=None)
    resid = prof[ok] - X @ coef
    s2 = np.sum(resid ** 2) / max(ok.sum() - 4, 1)
    cov = s2 * np.linalg.inv(X.T @ X)
    # empirical error: shift the template to off-filament positions and refit (placebo amplitudes)
    placebo = []
    for sh in np.arange(-1.6, 1.61, 0.1):
        if abs(sh) < 0.8:
            continue
        t3s = np.interp(cen - sh, cen, t3, left=0, right=0)
        Xs = np.c_[t3s[ok], np.ones(ok.sum()), cen[ok], cen[ok] ** 2]
        cs, *_ = np.linalg.lstsq(Xs, prof[ok], rcond=None)
        placebo.append(cs[0])
    on = ok & (np.abs(cen) <= 0.25)
    out.update(dict(bin_scatter=float(sig_bin), A_per_R=float(coef[0]), A_err_formal=float(np.sqrt(cov[0, 0])),
                    A_err_placebo=float(np.std(placebo)) if placebo else None, n_placebo=len(placebo),
                    on_mean_bgsub=float(np.nanmean(res[on])), o3_on_mean_R=float(np.nanmean(t3[on]))))
    return out


def run(band, rows=range(11, 28, 2)):
    gaia = load_gaia()
    sra, sdec = ridge_spine()
    o3maps = np.load(DER / 'imaging_optical_nsns_maps.npz')
    results = {}; binned = {}
    for d in rows:
        tl = []
        for col in ('W', 'E'):
            p = RAW / f'skyview/{band}/dss_{band}_d{d:02d}_{col}.fits'
            if p.exists():
                tl.append(process_tile(p, gaia, band))
        if not tl:
            continue
        prof = row_profile(tl, sra, sdec, o3maps)
        info = [dict(name=t['name'], k=t['k'], k_err=t['k_err'], n_cal=t['n_cal'], sky_D=t['sky_D'], pix_noise=t['pix_noise'],
                     masked_frac=t['masked_frac'], plates=t['plates']) for t in tl]
        # surface-brightness equivalents (AB mag/arcsec^2 in the PS1 calibration band)
        if 'bin_scatter' in prof:
            sb = prof['bin_scatter']
            prof['bin_scatter_mag'] = float(-2.5 * np.log10(sb) - 48.6 - 2.5 * np.log10(3631e-23)) if sb > 0 else None
            # flux unit: 10^(-0.4 m) per arcsec^2 -> mag = -2.5 log10(value)
            prof['bin_scatter_mag'] = float(-2.5 * np.log10(sb)) if sb > 0 else None
            prof['sky_mag'] = [float(-2.5 * np.log10(t['sky_D'] / (t['k'] * 36.0))) if t['k'] > 0 else None for t in tl]
            A = prof['A_per_R']; eA = prof['A_err_placebo'] or prof['A_err_formal']
            prof['A_per_R_mag_of_1R'] = float(-2.5 * np.log10(abs(A))) if A != 0 else None
            prof['A_3sigma_upper_mag_of_1R'] = float(-2.5 * np.log10(abs(A) + 3 * eA))
        results[f'd{d:02d}'] = dict(tiles=info, **prof)
        binned[f'd{d:02d}'] = [(t['binned'], t['wcs'].to_header_string(), t['sky_D'], t['k']) for t in tl]
        print(band, d, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in prof.items()
                        if k in ('A_per_R', 'A_err_formal', 'A_err_placebo', 'bin_scatter_mag', 'sky_mag', 'A_3sigma_upper_mag_of_1R')},
              [(t['name'], t['plates'], round(t['k'], 1) if np.isfinite(t['k']) else None, t['n_cal']) for t in tl], flush=True)
    return results, binned


if __name__ == '__main__':
    import pickle  # only for our own derived objects (never for downloaded data)
    bands = sys.argv[1:] or list(BANDS)
    allres = {}
    p = DER / 'imaging_optical_dss_profiles.json'
    if p.exists():
        allres = json.loads(p.read_text())
    for band in bands:
        res, binned = run(band)
        allres[band] = dict(plate=BANDS[band]['plate'], calibration_band=BANDS[band]['cal'], rows=res)
        np.savez_compressed(DER / f'imaging_optical_dss_binned_{band}.npz',
                            **{f'{k}_{i}_img': v[i][0] for k, v in binned.items() for i in range(len(v))},
                            **{f'{k}_{i}_wcs': np.array(v[i][1]) for k, v in binned.items() for i in range(len(v))},
                            **{f'{k}_{i}_sky': np.array([v[i][2], v[i][3]]) for k, v in binned.items() for i in range(len(v))})
    save_json(allres, p)
