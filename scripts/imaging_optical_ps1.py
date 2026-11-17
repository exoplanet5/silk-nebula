"""Pan-STARRS1 DR2 stack strips across the Silk Nebula: star-masked, 1'-binned cross-filament profiles.

Input : data/raw/ps1_dr2/cutouts/strip_<S>_<k>_<b>.fits (fitscut 12' cutouts block-averaged to 2.5"/pix)
Flux  : stack pixel values are counts with zero point ZP = 25 + 2.5 log10(EXPTIME) (PS1 DR2 stack convention);
        surface brightness per arcsec^2 = value * 10^(-0.4 ZP) / 0.0625 (binned values are means of native 0.25" pixels).
Mask  : Gaia DR3 G<17.5 with r = 4" + 70" 10^(-0.2 (G-6)) (cap 200"), 3-sigma outliers vs 1' median (dilated), empty pixels.
Profile: 2' bins of distance from the NSNS [OIII] ridge; quadratic background from 0.7<=|x|<=1.7 deg;
         template fit with the NSNS [OIII] (g) or NSNS Halpha (r) profile -> amplitude per R, placebo error.
Also measures the PS1 response to extended emission around the Medusa PN (Abell 21) as a sky-subtraction test.
Outputs: data/derived/imaging_optical_ps1_profiles.json, imaging_optical_ps1_strips.npz
"""
import sys
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import glob, json, re
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
from scipy.ndimage import median_filter, binary_dilation, gaussian_filter
from imaging_optical_common import *
from imaging_optical_nsns import ridge_spine
from imaging_optical_dss import load_gaia

PIXAS = 2.5
STRIPS = {'P3': 14.72, 'P1': 16.60, 'P2': 17.74, 'N22': 22.3}


def load_cut(path, gaia):
    with fits.open(path) as h:
        d = h[0].data.astype(float); hd = h[0].header
    w = WCS(hd)
    zp = 25.0 + 2.5 * np.log10(hd['EXPTIME'])
    f = d * 10 ** (-0.4 * zp) / 0.25 ** 2           # flux (10^-0.4 m) per arcsec^2; fitscut output_size AVERAGES native 0.25" pixels
    m = ~np.isfinite(d) | (d == 0)
    ny, nx = d.shape
    x, y = w.wcs_world2pix(gaia[:, 0], gaia[:, 1], 0)
    ok = (x > -80) & (x < nx + 80) & (y > -80) & (y < ny + 80)
    yy, xx = np.ogrid[0:ny, 0:nx]
    for xi, yi, g in zip(x[ok], y[ok], gaia[ok, 2]):
        r = max(min(4 + 70 * 10 ** (-0.2 * (g - 6)), 200) / PIXAS, 1.5)
        x0, x1 = int(max(0, xi - r - 1)), int(min(nx, xi + r + 2)); y0, y1 = int(max(0, yi - r - 1)), int(min(ny, yi + r + 2))
        if x1 > x0 and y1 > y0:
            m[y0:y1, x0:x1] |= (xx[:, x0:x1] - xi) ** 2 + (yy[y0:y1, :] - yi) ** 2 <= r * r
    dd = np.where(m, np.nan, f)
    bg = median_filter(np.nan_to_num(dd, nan=np.nanmedian(dd)), size=25)
    res = f - bg
    s = 1.4826 * np.nanmedian(np.abs(res[~m] - np.nanmedian(res[~m])))
    m |= binary_dilation((res > 3 * s) | (res < -5 * s), iterations=1)
    return np.where(m, np.nan, f), w, s, hd['EXPTIME']


def strip_profile(name, band, gaia, sra, sdec, nsns):
    files = sorted(glob.glob(str(RAW / f'ps1_dr2/cutouts/strip_{name}_*_{band}.fits')))
    if not files:
        return None, None
    dec0 = STRIPS[name]
    r0 = float(np.interp(dec0, sdec, sra))
    wl, shp = tan_wcs(r0, dec0, 3.6, 0.22, PIXAS / 60)
    acc = np.zeros(shp); cnt = np.zeros(shp); noise = []; expt = []
    for f in files:
        img, w, s, et = load_cut(f, gaia)
        noise.append(s); expt.append(et)
        arr, fp = reproject_interp((np.nan_to_num(img, nan=0.0), w), wl, shape_out=shp, order='nearest-neighbor')
        msk, _ = reproject_interp((np.isfinite(img).astype(float), w), wl, shape_out=shp, order='nearest-neighbor')
        good = (fp > 0) & (msk > 0.5)
        acc[good] += arr[good]; cnt[good] += 1
    mos = np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan)
    # 1' block medians
    b = 24
    ny2, nx2 = shp[0] // b, shp[1] // b
    blk = mos[:ny2 * b, :nx2 * b].reshape(ny2, b, nx2, b).transpose(0, 2, 1, 3).reshape(ny2, nx2, b * b)
    frac = np.isfinite(blk).mean(-1)
    with np.errstate(all='ignore'):
        med = np.nanmedian(blk, -1)
    med[frac < 0.25] = np.nan
    yy, xx = np.mgrid[0:ny2, 0:nx2]
    ra, dec = wl.wcs_pix2world(xx * b + (b - 1) / 2, yy * b + (b - 1) / 2, 0)
    across, along, beyond = along_across(ra, dec, sra, sdec)
    # NSNS template (OIII for g, Halpha for r)
    zw = WCS(str(nsns['wcs_header']))
    px, py = zw.wcs_world2pix(ra, dec, 0)
    ix = np.clip(np.round(px).astype(int), 0, nsns['oiii'].shape[1] - 1); iy = np.clip(np.round(py).astype(int), 0, nsns['oiii'].shape[0] - 1)
    tmpl_map = gaussian_filter(np.nan_to_num(nsns['oiii' if band == 'g' else 'halpha']), 1.0)
    tm = tmpl_map[iy, ix]
    edges = np.arange(-1.7, 1.7001, 2 / 60); cen = 0.5 * (edges[1:] + edges[:-1])
    prof = np.full(len(cen), np.nan); pt = np.full(len(cen), np.nan)
    for i in range(len(cen)):
        s = np.isfinite(med) & (across >= edges[i]) & (across < edges[i + 1])
        if s.sum() >= 2:
            prof[i] = np.median(med[s]); pt[i] = np.median(tm[s])
    ok = np.isfinite(prof) & np.isfinite(pt)
    bgm = ok & (np.abs(cen) >= 0.7)
    out = dict(strip=name, band=band, dec=dec0, n_cutouts=len(files), exptime_median=float(np.median(expt)),
               pix_noise_2p5as_mag=float(-2.5 * np.log10(np.median(noise))),
               across_deg=cen.tolist(), profile=prof.tolist(), template_R=pt.tolist(),
               template='NSNS [OIII]' if band == 'g' else 'NSNS Halpha')
    if bgm.sum() > 8 and ok.sum() > 12:
        c = np.polyfit(cen[bgm], prof[bgm], 2); res = prof - np.polyval(c, cen)
        sc = 1.4826 * np.nanmedian(np.abs(res[bgm] - np.nanmedian(res[bgm])))
        c3 = np.polyfit(cen[bgm], pt[bgm], 2); t3 = pt - np.polyval(c3, cen)
        X = np.c_[t3[ok], np.ones(ok.sum()), cen[ok], cen[ok] ** 2]
        coef, *_ = np.linalg.lstsq(X, prof[ok], rcond=None)
        placebo = []
        for sh in np.arange(-1.4, 1.41, 0.1):
            if abs(sh) < 0.8:
                continue
            t3s = np.interp(cen - sh, cen, t3, left=0, right=0)
            Xs = np.c_[t3s[ok], np.ones(ok.sum()), cen[ok], cen[ok] ** 2]
            cs, *_ = np.linalg.lstsq(Xs, prof[ok], rcond=None); placebo.append(cs[0])
        eA = float(np.std(placebo)) if len(placebo) > 3 else np.nan
        A = float(coef[0])
        out.update(profile_bgsub=res.tolist(), bin_scatter_mag=float(-2.5 * np.log10(sc)) if sc > 0 else None,
                   A_per_R=A, A_err_placebo=eA,
                   A_3sigma_upper_mag_of_1R=float(-2.5 * np.log10(abs(A) + 3 * eA)) if np.isfinite(eA) else None)
    return out, dict(img=med, wcs=wl[::b, ::b].to_header_string() if False else None, across=across)


def medusa_test(gaia):
    """Radial profile around Abell 21 in the 25' PS1 cutouts (sky-subtraction response check)."""
    out = {}
    for b in 'gr':
        p = RAW / f'ps1_dr2/cutouts/medusa_A21_{b}.fits'
        if not p.exists():
            continue
        img, w, s, et = load_cut(p, gaia)
        yy, xx = np.mgrid[0:img.shape[0], 0:img.shape[1]]
        ra, dec = w.wcs_pix2world(xx, yy, 0)
        r = np.hypot((ra - 112.2625) * np.cos(np.radians(13.2464)), dec - 13.2464) * 60
        edges = np.arange(0, 13.01, 0.5)
        prof = [float(np.nanmedian(img[(r >= edges[i]) & (r < edges[i + 1])])) for i in range(len(edges) - 1)]
        out[b] = dict(r_arcmin=(0.5 * (edges[1:] + edges[:-1])).tolist(), profile_flux_arcsec2=prof)
    return out


if __name__ == '__main__':
    gaia = load_gaia()
    sra, sdec = ridge_spine()
    nsns = np.load(DER / 'imaging_optical_nsns_maps.npz')
    res = {}
    for nm in STRIPS:
        for b in 'gr':
            o, _ = strip_profile(nm, b, gaia, sra, sdec, nsns)
            if o is None:
                continue
            res[f'{nm}_{b}'] = o
            print(nm, b, o['n_cutouts'], {k: (round(v, 2) if isinstance(v, float) else v) for k, v in o.items()
                                          if k in ('pix_noise_2p5as_mag', 'bin_scatter_mag', 'A_3sigma_upper_mag_of_1R')},
                  'A/eA', (round(o['A_per_R'] / o['A_err_placebo'], 2) if o.get('A_err_placebo') else None), flush=True)
    res['medusa_A21'] = medusa_test(gaia)
    save_json(res, DER / 'imaging_optical_ps1_profiles.json')
