"""verify_precovery step 3: image-level test with SDSS DR17 frame files of one run (4458, 2004-02-16).

For every frame (g and i; data/raw/verify_precovery_sdss/frames_4458/, list + log in urls.txt / download_log.txt):
  total = HDU0 (calibrated, sky-subtracted, nmgy/pix) + HDU2 global sky (ALLSKY counts, bilinear at XINTERP/YINTERP)
          x HDU1 calibvec (nmgy/count)        -> nmgy per pixel -> / 0.396^2 -> nmgy arcsec^-2
  Blocks of 64x64 px (25"x25"): 3-sigma-clipped median of total, of HDU0 alone, and mean of the sky model.
  Stars are removed by the clipping plus a mask of pixels > 4 sigma above the block median dilated by 3 px.
Model over the run (backfitting, my own code):
  y[b] = S(tau) + O(camcol, block column) + D_c(tau) + A * T(ra, dec)
  S: common mode per 64-row block row (same time in all camcols), O: flat-field residual per block column,
  D_c: per-camcol quadratic drift in tau, T: unit Gaussian ridges (E+W, FWHM 0.28 deg) at the block centre.
Outputs: data/derived/verify_precovery_frames_blocks.npz, verify_precovery_frames.json
"""
import sys, bz2, io, glob, json
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import binary_dilation
from verify_precovery_lib import MYRAW, DER, load_strands, make_model, band_response_per_R5007

RUN = 4458
B = 64
PIXAREA = 0.396 ** 2


def clipped_median(x, nsig=3.0, it=4):
    x = x[np.isfinite(x)]
    for _ in range(it):
        if x.size < 20:
            return np.nan
        m = np.median(x); s = 1.4826 * np.median(np.abs(x - m))
        k = np.abs(x - m) < nsig * max(s, 1e-12)
        if k.all():
            break
        x = x[k]
    return float(np.median(x))


def process_frame(path):
    raw = bz2.open(path).read()
    with fits.open(io.BytesIO(raw)) as h:
        img = h[0].data.astype(float); hdr = h[0].header
        cal = h[1].data.astype(float)
        sk = h[2].data
        allsky = sk['ALLSKY'][0].astype(float); xi = sk['XINTERP'][0].astype(float); yi = sk['YINTERP'][0].astype(float)
    ny, nx = img.shape
    # bilinear interpolation of the binned sky (constant extrapolation at the ends)
    ax = np.clip(xi, 0, allsky.shape[1] - 1); ay = np.clip(yi, 0, allsky.shape[0] - 1)
    x0 = np.floor(ax).astype(int); y0 = np.floor(ay).astype(int)
    x1 = np.minimum(x0 + 1, allsky.shape[1] - 1); y1 = np.minimum(y0 + 1, allsky.shape[0] - 1)
    fx = ax - x0; fy = ay - y0
    s00 = allsky[np.ix_(y0, x0)]; s01 = allsky[np.ix_(y0, x1)]; s10 = allsky[np.ix_(y1, x0)]; s11 = allsky[np.ix_(y1, x1)]
    sky = (s00 * (1 - fx)[None] * (1 - fy)[:, None] + s01 * fx[None] * (1 - fy)[:, None]
           + s10 * (1 - fx)[None] * fy[:, None] + s11 * fx[None] * fy[:, None])
    sky_nmgy = sky * cal[None, :]
    tot = (img + sky_nmgy) / PIXAREA
    imgs = img / PIXAREA
    # star mask from the sky-subtracted image
    med = np.median(imgs); sig = 1.4826 * np.median(np.abs(imgs - med))
    mask = binary_dilation(imgs > med + 4 * sig, iterations=3)
    w = WCS(hdr)
    rows = []
    nyb, nxb = ny // B, nx // B
    for j in range(nyb):
        for i in range(nxb):
            sl = (slice(j * B, (j + 1) * B), slice(i * B, (i + 1) * B))
            good = ~mask[sl]
            if good.mean() < 0.5:
                continue
            t = clipped_median(tot[sl][good]); m0 = clipped_median(imgs[sl][good])
            rows.append((i, j, t, float(sky_nmgy[sl].mean() / PIXAREA), m0, float(good.mean())))
    rows = np.array(rows)
    xc = rows[:, 0] * B + B / 2 - 0.5; yc = rows[:, 1] * B + B / 2 - 0.5
    ra, dec = w.wcs_pix2world(xc, yc, 0)
    return dict(ix=rows[:, 0].astype(int), iy=rows[:, 1].astype(int), tot=rows[:, 2], skymodel=rows[:, 3],
                img=rows[:, 4], fgood=rows[:, 5], ra=ra, dec=dec, frame_sigma=sig, frame_median_img=med)


def build_blocks():
    out = {k: [] for k in ('band', 'camcol', 'field', 'ix', 'iy', 'tot', 'skymodel', 'img', 'fgood', 'ra', 'dec')}
    files = sorted(glob.glob(str(MYRAW / f'frames_{RUN}/frame-*-{RUN:06d}-*.fits.bz2')))
    for k, p in enumerate(files):
        name = p.split('/')[-1]
        b = name.split('-')[1]; c = int(name.split('-')[3]); f = int(name.split('-')[4][:4])
        try:
            d = process_frame(p)
        except Exception as e:      # truncated download etc.
            print('skip', name, e); continue
        n = len(d['ra'])
        out['band'] += [b] * n; out['camcol'] += [c] * n; out['field'] += [f] * n
        for kk in ('ix', 'iy', 'tot', 'skymodel', 'img', 'fgood', 'ra', 'dec'):
            out[kk] += list(d[kk])
        if k % 40 == 0:
            print(k, len(files), name, flush=True)
    arr = {k: np.array(v) for k, v in out.items()}
    np.savez_compressed(DER / 'verify_precovery_frames_blocks.npz', **arr)
    return arr


def group_median(values, gi, ng, ok):
    """Median of values per group index gi (0..ng-1) using only ok entries (0 for empty groups)."""
    out = np.zeros(ng)
    idx = np.where(ok)[0]
    o = idx[np.argsort(gi[idx], kind='stable')]
    g = gi[o]
    cuts = np.flatnonzero(np.diff(g)) + 1
    for chunk in np.split(o, cuts):
        if len(chunk):
            out[gi[chunk[0]]] = np.median(values[chunk])
    return out


def backfit(y, tau, col, cam, T, n_iter=25, drift=2, fitA=True):
    """y = S(tau) + O(col) + D_cam(tau) + A T; returns A, components, residual (with A T left in)."""
    ok = np.isfinite(y)
    taus, ti = np.unique(tau, return_inverse=True)
    cols, ci = np.unique(col, return_inverse=True)
    cams, mi = np.unique(cam, return_inverse=True)
    tn = (tau - np.median(tau)) / np.ptp(tau)
    S = np.zeros(len(taus)); O = np.zeros(len(cols)); Dc = np.zeros((len(cams), drift)); A = 0.0
    P = np.stack([tn ** (p + 1) for p in range(drift)], 1)
    m = ok.copy()
    for it in range(n_iter):
        D = (P * Dc[mi]).sum(1)
        S = group_median(y - O[ci] - D - A * T, ti, len(taus), m)
        O = group_median(y - S[ti] - D - A * T, ci, len(cols), m)
        O -= np.median(O)
        r = y - S[ti] - O[ci] - A * T
        for k in range(1, len(cams)):
            mm = (mi == k) & m
            Dc[k], *_ = np.linalg.lstsq(P[mm], r[mm], rcond=None)
        D = (P * Dc[mi]).sum(1)
        r = y - S[ti] - O[ci] - D
        if fitA:
            A = np.sum(T[m] * r[m]) / np.sum(T[m] ** 2)
        e = r - A * T
        s = 1.4826 * np.median(np.abs(e[ok] - np.median(e[ok])))
        m = ok & (np.abs(e) < 4 * s)
    resid = y - S[ti] - O[ci] - (P * Dc[mi]).sum(1)          # filament left in (A*T not subtracted)
    e = resid - A * T
    sig = float(np.std(e[m]))
    return dict(A=float(A), resid=resid, S=S, O=O, Dc=Dc, sigma=sig, mask=m)


def analyse(arr):
    strands, _ = load_strands()
    res = {}
    for band in ('g', 'i'):
        m = (arr['band'] == band) & np.isfinite(arr['tot'])
        if m.sum() == 0:
            continue
        cam = arr['camcol'][m]; fld = arr['field'][m]; ix = arr['ix'][m]; iy = arr['iy'][m]
        tau = fld * 1361 + iy * B                      # along-scan row index (same time in all camcols)
        tau = (tau // B) * B
        col = cam * 100 + ix
        ra, dec = arr['ra'][m], arr['dec'][m]
        T = make_model('gauss028', strands)(ra, dec)
        conv = band_response_per_R5007(band)[0]
        dE, _, oE, _ = strands['E'].geometry(ra, dec)
        dW, _, _, _ = strands['W'].geometry(ra, dec)
        out = {}
        for comp in ('tot', 'skymodel', 'img'):
            y = arr[comp][m]
            fitr = backfit(y, tau, col, cam, T)
            # error from shifted-template nulls (template displaced in RA, fitted on the same residual)
            r_noA = fitr['resid']
            mm = fitr['mask']
            nulls = []
            for s in (-2.0, -1.6, -1.2, 1.2, 1.6, 2.0):
                Ts = make_model('gauss028', strands, ra_shift=s)(ra, dec)
                ee = r_noA - fitr['A'] * T
                if np.sum(Ts[mm] ** 2) > 0:
                    nulls.append(np.sum(Ts[mm] * ee[mm]) / np.sum(Ts[mm] ** 2))
            # cross-ridge profile (distance to E spine), Dec bands
            edges = np.arange(-2.0, 2.0001, 0.1)
            prof = {}
            for d0, d1 in ((14.5, 16.0), (16.0, 17.5), (17.5, 19.5), (14.5, 19.5)):
                sel = mm & (dec >= d0) & (dec < d1) & (oE == 0)
                p, pe = [], []
                for j in range(len(edges) - 1):
                    k = sel & (dE >= edges[j]) & (dE < edges[j + 1])
                    v = r_noA[k]
                    p.append(float(np.median(v)) if k.sum() > 10 else np.nan)
                    pe.append(float(1.2533 * np.std(v) / np.sqrt(k.sum())) if k.sum() > 10 else np.nan)
                prof[f'{d0}-{d1}'] = dict(p=(np.array(p) / conv).tolist(), e=(np.array(pe) / conv).tolist())
            out[comp] = dict(A_R=float(fitr['A'] / conv), null_R=[float(v / conv) for v in nulls],
                             null_std_R=float(np.std(nulls) / conv), block_sigma_R=float(fitr['sigma'] / conv),
                             n_blocks=int(mm.sum()), profile=prof, profile_edges=edges.tolist())
            if comp == 'tot':
                np.savez_compressed(DER / f'verify_precovery_frames_resid_{band}.npz', ra=ra, dec=dec, resid=r_noA / conv,
                                    model=fitr['A'] * T / conv, mask=mm, camcol=cam, field=fld, ix=ix, iy=iy, dE=dE, dW=dW)
        out['conv'] = conv
        res[band] = out
        print(band, {c: (round(out[c]['A_R'], 3), round(out[c]['null_std_R'], 3)) for c in ('tot', 'skymodel', 'img')})
    return res


def compare_psp(arr):
    """Field-mean of my block 'tot' vs the Field table sky estimators (linearity / plateau check)."""
    from verify_precovery_lib import load_fields
    a = load_fields()
    m = a['run'] == RUN
    out = {}
    for band in ('g', 'i'):
        rows = []
        for c in range(1, 7):
            for f in np.unique(arr['field'][(arr['band'] == band) & (arr['camcol'] == c)]):
                k = (arr['band'] == band) & (arr['camcol'] == c) & (arr['field'] == f) & np.isfinite(arr['tot'])
                j = np.where(m & (a['camcol'] == c) & (a['field'] == f))[0]
                if k.sum() < 100 or len(j) == 0:
                    continue
                j = j[0]
                rows.append((c, f, np.median(arr['tot'][k]), np.mean(arr['tot'][k]), np.mean(arr['skymodel'][k]),
                             a[f'sky_{band}'][j], a[f'skyFrames_{band}'][j], a[f'skyFramesSub_{band}'][j]))
        rows = np.array(rows)
        res = {}
        for name, col in (('psp', 5), ('frames', 6), ('framessub', 7), ('globalsky_mean', 4)):
            # remove per-camcol offsets, then regress on my block-median field value
            x = rows[:, 2].copy(); y = rows[:, col].copy()
            for c in range(1, 7):
                k = rows[:, 0] == c
                x[k] -= np.median(x[k]); y[k] -= np.median(y[k])
            sl = np.polyfit(x, y, 1)
            res[name] = dict(slope=float(sl[0]), rms=float(np.std(y - np.polyval(sl, x))),
                             corr=float(np.corrcoef(x, y)[0, 1]))
        out[band] = res
        np.save(DER / f'verify_precovery_frames_fieldcompare_{band}.npy', rows)
    return out


if __name__ == '__main__':
    import os
    p = DER / 'verify_precovery_frames_blocks.npz'
    arr = dict(np.load(p)) if (p.exists() and '--rebuild' not in sys.argv) else build_blocks()
    res = {'run': RUN, 'block_px': B, 'analysis': analyse(arr), 'field_compare': compare_psp(arr)}
    print(json.dumps(res['field_compare'], indent=1))
    (DER / 'verify_precovery_frames.json').write_text(json.dumps(res, indent=1))
