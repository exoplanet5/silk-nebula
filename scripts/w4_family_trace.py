"""w4_family_trace.py -- wave-4 `family` step 1b: trace Monogem-periphery [O III] filaments and the G206.6+6.1 rim
in the NSNS DR0.2 [O III] order-3 HiPS (linear, Rayleigh), and measure their NSNS ridge brightness and width.

Filaments (guides read off Fesen+2026 Figs 31/32 (v2 numbering) and the NSNS mosaic; ridge positions are then
measured, not copied):
  S1  Monogem south-limb [O III] filament (Fesen+26 Fig 32 "Monogem SNR"), RA 102-113.4, Dec 3.4-5.0
  E1  short faint filament ~3 deg east of the Silk (Fesen+26 6.2.1 "a bit farther east ... much fainter"),
      RA ~119.5, Dec 13.5-16.6
  G206  rim of G206.6+6.1 ([O III]-dominated SNR, ~440 pc; control), traced radially about (105.6, +8.3)
Method: 1'/px TAN mosaic, 1-deg median background removed, Gaussian smoothing sigma 1.5'. Along the guide (0.05-deg
steps) the peak offset along the local normal is searched within +-w; S/N against the robust scatter of the outer
profile. Offsets are median-smoothed (5 samples); points with S/N >= 4 are kept. Width: Gaussian fit to the stack of
unsmoothed profiles aligned on the ridge.
Outputs: data/derived/w4_family_traces.json, figures/w4_family_traces.png
Usage: ~/.venvs/astro313/bin/python -I scripts/w4_family_trace.py
"""
import sys
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import json
import numpy as np
from scipy.ndimage import median_filter, gaussian_filter, map_coordinates, median_filter as mf1
from scipy.optimize import curve_fit
from imaging_optical_common import tan_wcs, hips_mosaic
from silk_style import *

ROOT = '/Users/mickey/astronomy/silk-nebula'
GUIDES = {
    'S1': dict(box=(107.8, 4.5, 14.0, 5.0), w=0.35, pts=[(102.1, 5.0), (103.0, 4.45), (104.0, 3.95), (105.0, 3.62),
                                                         (106.0, 3.47), (107.0, 3.50), (108.0, 3.58), (109.0, 3.68),
                                                         (110.0, 3.72), (111.0, 3.85), (112.0, 4.05), (113.4, 4.15)]),
    'E1': dict(box=(119.5, 15.0, 4.0, 5.0), w=0.40, pts=[(119.55, 13.3), (119.45, 14.5), (119.45, 15.5),
                                                         (119.5, 16.8)]),
}
G206_C = (105.6, 8.3)


def mosaic(box):
    W, S = tan_wcs(box[0], box[1], box[2], box[3], 1.0)
    m = hips_mosaic('oiii', 3, W, S)
    sn = hips_mosaic('sn_oiii', 3, W, S)
    bg = median_filter(np.nan_to_num(m, nan=np.nanmedian(m)), size=61)
    raw = m - bg
    sm = gaussian_filter(np.nan_to_num(raw), 1.5)
    sm[~np.isfinite(m)] = np.nan
    return W, raw, sm, sn


def tangent_frame(ra, dec):
    """Local east/north unit offsets (deg) -> use flat-sky with cos(dec)."""
    return np.cos(np.radians(dec))


def sample(W, img, ra, dec):
    x, y = W.wcs_world2pix(ra, dec, 0)
    return map_coordinates(np.nan_to_num(img, nan=0.0), [np.atleast_1d(y), np.atleast_1d(x)], order=1, mode='nearest')


def trace_polyline(W, raw, sm, pts, w, step=0.05):
    pts = np.array(pts, float)
    c = np.cos(np.radians(pts[:, 1]))
    xy = np.c_[pts[:, 0] * c, pts[:, 1]]           # flat approx, fine over a few deg
    seg = np.hypot(*np.diff(xy, axis=0).T)
    s = np.r_[0, np.cumsum(seg)]
    ss = np.arange(0, s[-1], step)
    gx = np.interp(ss, s, xy[:, 0]); gy = np.interp(ss, s, xy[:, 1])
    tx = np.gradient(gx); ty = np.gradient(gy)
    tn = np.hypot(tx, ty); tx /= tn; ty /= tn
    nx, ny = -ty, tx                                 # normal
    off = np.arange(-1.0, 1.0001, 1 / 60)
    out = []
    profs = []
    for k in range(len(ss)):
        cd = np.cos(np.radians(gy[k]))
        ra = (gx[k] + off * nx[k]) / cd
        dec = gy[k] + off * ny[k]
        p = sample(W, sm, ra, dec)
        pr = sample(W, raw, ra, dec)
        win = np.abs(off) <= w
        j = np.argmax(np.where(win, p, -np.inf))
        outer = np.abs(off) > w + 0.15
        noise = 1.4826 * np.median(np.abs(p[outer] - np.median(p[outer])))
        base = np.median(p[outer])
        out.append((ss[k], off[j], p[j] - base, noise, gx[k], gy[k], nx[k], ny[k]))
        profs.append(pr - np.median(pr[outer]))
    a = np.array(out)
    offs = mf1(a[:, 1], size=5, mode='nearest')
    sn = a[:, 2] / np.where(a[:, 3] > 0, a[:, 3], np.inf)
    ridge_x = a[:, 4] + offs * a[:, 6]; ridge_y = a[:, 5] + offs * a[:, 7]
    ra_r = ridge_x / np.cos(np.radians(ridge_y))
    # aligned stack for width
    stack, n = np.zeros_like(off), 0
    for k in range(len(ss)):
        if sn[k] >= 4:
            stack += np.interp(off, off - offs[k], profs[k]); n += 1
    stack /= max(n, 1)
    return dict(s=a[:, 0], ra=ra_r, dec=ridge_y, peak_R=a[:, 2], noise_R=a[:, 3], sn=sn), off, stack


def gfit(off, prof):
    f = lambda x, A, mu, s, c: A * np.exp(-0.5 * ((x - mu) / s) ** 2) + c   # noqa: E731
    m = np.abs(off) < 0.8
    p, cov = curve_fit(f, off[m], prof[m], p0=[np.max(prof), 0, 0.05, 0])
    return dict(peak_R=float(p[0]), mu=float(p[1]), sigma_deg=float(abs(p[2])), fwhm_deg=float(2.3548 * abs(p[2])),
                sigma_err=float(np.sqrt(cov[2, 2])))


def trace_ring(W, raw, sm, c0, r0=0.9, r1=2.2, step_pa=2.0):
    out, profs = [], []
    rr = np.arange(0.5, 2.6, 1 / 60)
    for pa in np.arange(0, 360, step_pa):
        t = np.radians(pa)
        dec = c0[1] + rr * np.cos(t)
        ra = c0[0] + rr * np.sin(t) / np.cos(np.radians(dec))
        p = sample(W, sm, ra, dec); pr = sample(W, raw, ra, dec)
        win = (rr >= r0) & (rr <= r1)
        j = np.argmax(np.where(win, p, -np.inf))
        outer = (rr < 0.7) | (rr > 2.4)
        noise = 1.4826 * np.median(np.abs(p - np.median(p)))
        base = np.median(p[outer]) if outer.sum() > 5 else np.median(p)
        out.append((pa, rr[j], p[j] - base, noise, ra[j], dec[j]))
        profs.append(pr - np.median(pr))
    a = np.array(out)
    sn = a[:, 2] / a[:, 3]
    return dict(pa=a[:, 0], r=a[:, 1], ra=a[:, 4], dec=a[:, 5], peak_R=a[:, 2], noise_R=a[:, 3], sn=sn), rr, profs


def main():
    res = {'method': __doc__.split('Method:')[1].split('Outputs')[0].strip(), 'filaments': {}}
    fig, axes = new_sky_figure(16, 9.6, 'Monogem-periphery [O III] filaments · NSNS DR0.2 ridge traces',
                               'S1 south limb, E1 east of the Silk, G206.6+6.1 rim (control); dots = ridge with S/N ≥ 4',
                               nrows=2, ncols=2)
    axes = axes.ravel()
    for k, (name, g) in enumerate(GUIDES.items()):
        W, raw, sm, sn = mosaic(g['box'])
        tr, off, stack = trace_polyline(W, raw, sm, g['pts'], g['w'])
        good = tr['sn'] >= 4
        wfit = gfit(off, stack)
        res['filaments'][name] = dict(
            guide=g['pts'], n_samples=int(len(tr['s'])), n_good=int(good.sum()),
            ridge=[dict(ra=round(float(a), 4), dec=round(float(d), 4), peak_R=round(float(p), 3),
                        sn=round(float(s), 1)) for a, d, p, s in zip(tr['ra'], tr['dec'], tr['peak_R'], tr['sn'])],
            ridge_peak_R_median_good=float(np.median(tr['peak_R'][good])) if good.any() else None,
            ridge_peak_R_p84_good=float(np.percentile(tr['peak_R'][good], 84)) if good.any() else None,
            stack_fit=wfit, length_good_deg=float(0.05 * good.sum()))
        ax = axes[k]
        ext = W.wcs_pix2world([[0, 0], [raw.shape[1] - 1, raw.shape[0] - 1]], 0)
        ax.imshow(sm, origin='lower', cmap='gray', vmin=-0.2, vmax=1.0, extent=[ext[0, 0], ext[1, 0], ext[0, 1], ext[1, 1]],
                  aspect=1 / np.cos(np.radians(g['box'][1])))
        ax.plot(tr['ra'][good], tr['dec'][good], '.', color=OIII_C, ms=3)
        ax.set_title(f'{name}: ridge median {res["filaments"][name]["ridge_peak_R_median_good"]:.2f} R (smoothed), '
                     f'FWHM {wfit["fwhm_deg"]:.3f}°', color=FG)
        ax.set_xlabel('RA (deg)'); ax.set_ylabel('Dec (deg)')
        ax3 = axes[3]
        ax3.plot(off, stack, label=f'{name} (n={good.sum()})')
    # ring
    W, raw, sm, sn = mosaic((105.6, 8.3, 5.5, 5.5))
    rg, rr, profs = trace_ring(W, raw, sm, G206_C)
    good = rg['sn'] >= 4
    res['filaments']['G206'] = dict(centre=G206_C, n_good=int(good.sum()),
                                    ridge=[dict(pa=float(p), r=round(float(r), 3), ra=round(float(a), 4),
                                                dec=round(float(d), 4), peak_R=round(float(q), 3), sn=round(float(s), 1))
                                           for p, r, a, d, q, s in zip(rg['pa'], rg['r'], rg['ra'], rg['dec'],
                                                                       rg['peak_R'], rg['sn'])],
                                    ridge_peak_R_median_good=float(np.median(rg['peak_R'][good])),
                                    r_median_good=float(np.median(rg['r'][good])))
    ax = axes[2]
    ext = W.wcs_pix2world([[0, 0], [raw.shape[1] - 1, raw.shape[0] - 1]], 0)
    ax.imshow(sm, origin='lower', cmap='gray', vmin=-0.2, vmax=1.0, extent=[ext[0, 0], ext[1, 0], ext[0, 1], ext[1, 1]],
              aspect=1 / np.cos(np.radians(8.3)))
    ax.plot(rg['ra'][good], rg['dec'][good], '.', color=WARN, ms=4)
    ax.set_title(f'G206.6+6.1 rim: median peak {res["filaments"]["G206"]["ridge_peak_R_median_good"]:.2f} R', color=FG)
    ax.set_xlabel('RA (deg)'); ax.set_ylabel('Dec (deg)')
    ax3 = axes[3]
    ax3.set_xlabel('offset across ridge (deg)'); ax3.set_ylabel('NSNS [O III] (R, 1° bg removed)')
    ax3.set_xlim(-0.8, 0.8); ax3.legend(); ax3.set_title('Aligned cross-profiles (unsmoothed)', color=FG)
    ax3.grid(alpha=0.3)
    finish(fig, f'{ROOT}/figures/w4_family_traces.png')
    json.dump(res, open(f'{ROOT}/data/derived/w4_family_traces.json', 'w'), indent=1)
    for n, f in res['filaments'].items():
        print(n, f.get('n_good'), f.get('ridge_peak_R_median_good'), f.get('stack_fit'))


if __name__ == '__main__':
    main()
