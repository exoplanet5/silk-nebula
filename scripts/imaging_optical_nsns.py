"""NSNS DR0.2 [OIII]/Halpha/[SII] (linear, Rayleighs) along the Silk Nebula.

1. Mosaics order-3 HiPS tiles (51.5" pixels) onto a 1'/pix TAN grid around the filament.
2. Traces the [OIII] ridge row by row (Dec) and fits a smooth polynomial RA(Dec).
3. Measures cross-filament profiles and box intensities at P1-P3, the bright south section, the north
   Halpha patch and P4; noise from the sn_* variance HiPS (normalized to 30" pixels) and from off-filament scatter.
Outputs: data/derived/imaging_optical_nsns_ridge.json, imaging_optical_nsns_measurements.json,
         data/derived/imaging_optical_nsns_maps.npz (for figures)
"""
import sys
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import numpy as np
from scipy.ndimage import median_filter, gaussian_filter
from imaging_optical_common import *

PIX = 1.0  # arcmin
w, shp = tan_wcs(116.3, 18.5, 12.0, 22.0, PIX)   # RA ~110-122.6, Dec 7.5-29.5


def build_maps():
    maps = {}
    for n in ['oiii', 'halpha', 'sii', 'vc', 'sn_oiii', 'sn_halpha', 'sn_sii']:
        maps[n] = hips_mosaic(n, 3, w, shp)
    return maps


def trace_ridge(o3, ra, dec):
    """Ridge of the brightest [OIII] strand. Anisotropic Gaussian matched filter (sigma 3' across, 12' along
    the roughly N-S filament), then per pixel row (1' in Dec) the peak within +-0.5 deg of the previous row's
    peak, starting from the bright section at Dec 15. Rows whose peak S/N < 4 are flagged."""
    img = np.nan_to_num(o3)
    sm = gaussian_filter(img, (12, 3))               # (rows=Dec, cols=RA) in 1' pixels
    bgs = median_filter(sm, size=(1, 121))           # 2-deg running median along rows = local background
    res = sm - bgs
    ny, nx = img.shape
    # robust noise of the filtered map away from the filament: MAD of residual in each row
    rows_out = {}
    x0 = int(np.argmin(np.abs(ra[int(np.argmin(np.abs(dec[:, nx // 2] - 15.0))), :] - 116.55)))
    ystart = int(np.argmin(np.abs(dec[:, x0] - 15.0)))
    for direction in (1, -1):
        xg = x0
        y = ystart
        while 0 <= y < ny:
            lo, hi = max(0, xg - 30), min(nx, xg + 31)
            seg = res[y, lo:hi]
            k = int(np.argmax(seg)) + lo
            noise = 1.4826 * np.median(np.abs(res[y] - np.median(res[y])))
            snr = res[y, k] / max(noise, 1e-4)
            if snr >= 4:
                xg = k
            rows_out[y] = (k, res[y, k], snr)
            y += direction
    out = []
    for y in sorted(rows_out):
        k, amp, snr = rows_out[y]
        out.append(dict(dec=float(dec[y, k]), ra=float(ra[y, k]), peak_R_filtered=float(amp), snr=float(snr)))
    return out


def main():
    import os
    if os.environ.get('NSNS_REUSE'):
        z = np.load(DER / 'imaging_optical_nsns_maps.npz'); maps = {k: z[k] for k in z.files if k not in ('ra','dec','wcs_header')}
    else:
        maps = build_maps()
    ra, dec = pixel_radec(w, shp)
    np.savez_compressed(DER / 'imaging_optical_nsns_maps.npz', ra=ra, dec=dec, **maps,
                        wcs_header=np.array(w.to_header_string()))
    rows = trace_ridge(maps['oiii'], ra, dec)
    good = [r for r in rows if r['snr'] > 5]
    d = np.array([r['dec'] for r in good]); a = np.array([r['ra'] for r in good])
    # robust cubic fit RA(Dec)
    keep = np.ones(len(d), bool)
    for it in range(4):
        p = np.polyfit(d[keep], a[keep], 4)
        res = (a - np.polyval(p, d)) * np.cos(np.radians(d))
        s = 1.4826 * np.median(np.abs(res[keep]))
        keep = np.abs(res) < max(3 * s, 0.05)
    ridge = dict(
        description='NSNS DR0.2 [OIII] ridge of the Silk Nebula, measured by imaging_optical (row-wise peak of the '
                    '4-arcmin-smoothed [OIII] HiPS order-3 map). Independent of the geometry teammate trace.',
        source='https://simg.de/nebulae3/dr0_2/oiii (HiPS order 3, 51.5 arcsec pixels), accessed 2026-11-16',
        poly_ra_of_dec=p.tolist(), dec_range=[float(d[keep].min()), float(d[keep].max())],
        rms_deg=float(s), points=rows)
    save_json(ridge, DER / 'imaging_optical_nsns_ridge.json')
    print('ridge dec range', ridge['dec_range'], 'rms', s, 'n', keep.sum())
    for r in rows[::60]:
        print(r)


def ridge_spine(ridge=None, step=0.05):
    """Dense (ra, dec) spine from the fitted ridge polynomial."""
    import json
    if ridge is None:
        ridge = json.loads((DER / 'imaging_optical_nsns_ridge.json').read_text())
    d = np.arange(ridge['dec_range'][0], ridge['dec_range'][1] + 1e-6, step)
    return np.polyval(ridge['poly_ra_of_dec'], d), d


SECTIONS = {   # name: (dec_center, half-length along filament in deg) ; P1-P3 from paper Table 1
    'P3': (14.72, 0.30), 'P1': (16.60, 0.30), 'P2': (17.74, 0.30),
    'south_bright': (13.75, 0.75), 'mid_17_19': (18.5, 0.5),
    'north_arc_21_23': (22.0, 1.0), 'north_arc_24_27': (25.5, 1.5),
}


def measure(maps, ra, dec, sra, sdec):
    across, along, beyond = along_across(ra, dec, sra, sdec)
    o = np.argsort(sdec)
    dec_on_spine = np.interp(along, np.concatenate([[0], np.cumsum(np.hypot(np.diff(sra[o]) * np.cos(np.radians(sdec[o][1:])), np.diff(sdec[o])))]), sdec[o])
    res = {}
    edges = np.arange(-1.5, 1.5001, 2 / 60)
    cen = 0.5 * (edges[1:] + edges[:-1])
    for name, (dc, hl) in SECTIONS.items():
        sel = (np.abs(dec_on_spine - dc) <= hl) & ~beyond
        r = {'dec_center': dc, 'half_length_deg': hl, 'across_bin_deg': cen.tolist()}
        for ch in ('oiii', 'halpha', 'sii', 'vc'):
            img = maps[ch]
            ok = sel & np.isfinite(img)
            # star residual clipping: drop pixels with sn variance > 99th percentile (bright-star contamination)
            var = maps.get('sn_' + ch)
            if var is not None:
                vv = var[ok]; thr = np.nanpercentile(vv, 98)
                ok &= ~(var > thr)
            prof = np.array([np.nanmedian(img[ok & (across >= edges[k]) & (across < edges[k + 1])])
                             for k in range(len(cen))])
            r[ch + '_profile_R'] = prof.tolist()
            on = ok & (np.abs(across) <= 0.25)
            off = ok & (np.abs(across) >= 0.8) & (np.abs(across) <= 1.5)
            on_v = img[on]; off_v = img[off]
            # uncertainty: scatter of 0.5 x 0.1 deg sub-box means in the off region
            nb = max(int(off.sum() / (30 * 6)), 1)
            sig_off = 1.4826 * np.median(np.abs(off_v - np.median(off_v))) / np.sqrt(max(off.sum() / 9, 1))
            r[ch] = dict(on_median_R=float(np.median(on_v)), off_median_R=float(np.median(off_v)),
                         excess_R=float(np.median(on_v) - np.median(off_v)),
                         peak_excess_R=float(np.nanmax(prof[np.abs(cen) <= 0.3]) - np.median(off_v)),
                         sigma_stat_R=float(sig_off), n_on=int(on.sum()), n_off=int(off.sum()),
                         sys_EW_R=float(0.5 * abs(np.median(img[ok & (across >= 0.8) & (across <= 1.5)]) -
                                                  np.median(img[ok & (across <= -0.8) & (across >= -1.5)]))))
        res[name] = r
    # P4 (paper box) relative to its surroundings, and the north Halpha patch (geometry v1.0 polygon)
    pat = in_patch(ra, dec)
    pc = (115.50, 20.55)
    dp = np.hypot((ra - pc[0]) * np.cos(np.radians(pc[1])), dec - pc[1])
    dra = (ra - 115.20) * np.cos(np.radians(16.18)); dde = dec - 16.18
    p4 = (np.abs(dra) <= 0.15) & (np.abs(dde) <= 0.40)
    p4ann = (np.abs(dra) <= 1.15) & (np.abs(dde) <= 1.4) & ~p4 & (np.abs(across) > 0.6)
    regions = {'P4': (p4, p4ann), 'Ha_patch_north': (pat, (dp > 0.9) & (dp < 2.0) & (np.abs(across) > 0.6))}
    for name, (inb, ann) in regions.items():
        r = {'definition': 'P4: paper box 0.3x0.8 deg at (115.20,16.18)' if name == 'P4' else
             'geometry v1.0 halpha_patch polygon; annulus 0.9-2.0 deg from (115.50,20.55), >0.6 deg from ridge'}
        for ch in ('oiii', 'halpha', 'sii', 'vc'):
            img = maps[ch]
            a1 = img[inb & np.isfinite(img)]; a2 = img[ann & np.isfinite(img)]
            sig = 1.4826 * np.median(np.abs(a2 - np.median(a2))) / np.sqrt(max(inb.sum() / 9, 1))
            r[ch] = dict(box_median_R=float(np.median(a1)), annulus_median_R=float(np.median(a2)),
                         excess_R=float(np.median(a1) - np.median(a2)), sigma_stat_R=float(sig),
                         annulus_scatter_R=float(1.4826 * np.median(np.abs(a2 - np.median(a2)))), n=int(inb.sum()))
        res[name] = r
    return res


if __name__ == '__main__':
    main()
    import json
    z = np.load(DER / 'imaging_optical_nsns_maps.npz')
    maps = {k: z[k] for k in z.files if k not in ('ra', 'dec', 'wcs_header')}
    sra, sdec = ridge_spine()
    res = measure(maps, z['ra'], z['dec'], sra, sdec)
    res['_notes'] = ('NSNS DR0.2 HiPS order 3 (51.5" pixels) mosaicked to 1 arcmin; on = |across|<=0.25 deg of the '
                     'NSNS [OIII] ridge, off = 0.8-1.5 deg both sides; medians; sigma_stat from off-region MAD '
                     'assuming 3x3-pixel (~2.6 arcmin) correlation. NSNS [OIII]/[SII] calibration is approximate '
                     '(star-colour based, see simg.de/nebulae3/dr0_2); large structures > ~3 deg are suppressed.')
    save_json(res, DER / 'imaging_optical_nsns_measurements.json')
    for k, v in res.items():
        if k.startswith('_'):
            continue
        print(k, {ch: (round(v[ch]['excess_R'], 2), round(v[ch].get('peak_excess_R', np.nan), 2), round(v[ch]['sigma_stat_R'], 3), round(v[ch].get('sys_EW_R', v[ch].get('annulus_scatter_R', np.nan)), 2)) for ch in ('oiii', 'halpha', 'sii', 'vc')})
