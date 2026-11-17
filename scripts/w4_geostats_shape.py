"""w4_geostats task 3: cross-profile shape of the NSNS [O III] ridge (south symmetric? north one-sided?) with one
estimator, placement + injection nulls, and a geometric sheet family that spans flat ribbon -> curved limb.

Estimator (primary): w4_geostats_common.stack_profile (block medians, 0.25-deg blocks, 0.025-deg bins), split normal
+ linear baseline, unweighted, |x| <= 0.9 deg -> r = sigma_W / sigma_E; non-parametric R_asym = <I-b>(-0.45..-0.15) /
<I-b>(+0.15..+0.45) with b the fitted baseline. Variants: single-block median (judge_skeptic-like), 0.02-deg bins with
bootstrap-weighted fit (judge_obs-like), fit window 0.6 deg, block means. Errors: block bootstrap (300).

Geometric model. A thin emitting layer (depth 0..w behind the front, top hat) on a sheet with curvature radius R_c,
extending a distance l along the line of sight from the tangent point. In the small-angle limit a point at LOS
coordinate z and depth u projects to h = u + z^2 / (2 R_c), so the projected column is
  N(h) ~ sqrt(min(h, q w)) - sqrt(max(h - w, 0)),  0 <= h <= w (1 + q),   q = l^2 / (2 R_c w)  (LOS sagitta / depth).
q -> 0 is a flat ribbon (top hat), q >> 1 a curved limb (sqrt rise, 1/sqrt(h) tail toward the centre of curvature).
A ONE-SIDED CHORD (only the near or far half, z in [0, l]) has exactly the same N(h) and half the column, so profile
shape cannot separate it from a full limb. The front is on the EAST; the tail runs west. The profile is convolved
with a Gaussian blur sigma_b (0.03, 0.06 or 0.10 deg: NSNS PSF + trace error + stacking of a wiggly ridge) and w is solved so the model's
FWHM equals the observed stacked FWHM of the segment. The peak limb gain relative to a face-on layer is
G = 2 sqrt(2 q R_c / w) for q <= 1 (both halves), 2 sqrt(2 R_c / w) for q >= 1.
Nulls: each model (plus a symmetric Gaussian of the observed width) is injected at pixel level, with the observed
amplitude, into 48 off-Silk NSNS placements (dRA +-1.5..5, dDec -2/0/+2; pixels within 1.0 deg of the real spine
excluded) and re-measured -> distribution of r per model -> likelihood of the observed r (Gaussian of the
injected-null mean and sd) and empirical p.
Segments: S1 Dec 12-15 (clean south), N Dec 20-28, N1 Dec 20-21.5, N3 Dec 21.5-24.25 (W2 0.6-1.3 deg W),
N2 Dec 24.25-28; north segments mask pixels within 0.9 deg of the Halpha patch centroid.
Output: data/derived/w4_geostats_shape.json, w4_geostats_shape_profiles.npz
Run: ~/.venvs/astro313/bin/python -I scripts/w4_geostats_shape.py
"""
import sys
import numpy as np
from scipy.ndimage import gaussian_filter1d
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from w4_geostats_common import (CC0, Placed, stack_profile, fit_split, shifted_cc, load_nsns, empirical_p, robust_sd,  # noqa
                                DER, json_dump, sep_pa, HA_PATCH, SEGMENTS)

RNG = np.random.default_rng(41173)
OUT = DER / 'w4_geostats_shape.json'
SEGMENTS.update({'N1': (20.0, 21.5), 'N3': (21.5, 24.25)})
SEGS = ('S1', 'N', 'N1', 'N3', 'N2')
E25 = np.arange(-1.2, 1.2001, 0.025); X25 = 0.5 * (E25[1:] + E25[:-1])
E20 = np.arange(-1.2, 1.2001, 0.02); X20 = 0.5 * (E20[1:] + E20[:-1])
NSHIFTS = [(a, b) for a in np.r_[-5:-1.4:0.5, 1.5:5.01:0.5] for b in (-2.0, 0.0, 2.0)]
QGRID = (0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 10.0)
SBS = (0.03, 0.06, 0.10)
DXF = 0.002
XF = np.arange(-2.0, 2.0001, DXF)
RC_PC = 58.0            # Monogem limb curvature used in closure W (verify_shockmodel)
PC_PER_DEG = 282.0 * np.pi / 180


def ra_asym(x, y, f):
    if f is None:
        return np.nan
    b = f['b0'] + f['b1'] * x
    w = (x >= -0.45) & (x <= -0.15); e = (x >= 0.15) & (x <= 0.45)
    yw = np.nanmean((y - b)[w]); ye = np.nanmean((y - b)[e])
    return float(yw / ye) if ye > 0 else np.nan


# ------------------------------------------------------------------------------------------------ models
def raw_profile(q, w):
    """Unconvolved projected column vs x (front at x = 0, emission at x < 0), unit peak."""
    h = -XF
    qw = max(q * w, 1e-6)
    N = np.sqrt(np.clip(np.minimum(h, qw), 0, None)) - np.sqrt(np.clip(h - w, 0, None))
    N[(h < 0) | (h > w + qw)] = 0
    N = np.clip(N, 0, None)
    return N / N.max()


def conv(prof, sb):
    p = gaussian_filter1d(prof, sb / DXF, mode='constant')
    return p / p.max()


def fwhm_of(p):
    above = np.nonzero(p >= 0.5)[0]
    return (above[-1] - above[0]) * DXF


def model(q, sb, fwhm_target):
    """Solve the layer depth w so the convolved model has the target FWHM; returns (profile peak-centred, w)."""
    lo, hi = 0.002, 1.5
    if fwhm_of(conv(raw_profile(q, lo), sb)) > fwhm_target:
        w = lo
    else:
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if fwhm_of(conv(raw_profile(q, mid), sb)) > fwhm_target:
                hi = mid
            else:
                lo = mid
        w = 0.5 * (lo + hi)
    p = conv(raw_profile(q, w), sb)
    xpk = XF[np.argmax(p)]
    return (lambda xx: np.interp(xx + xpk, XF, p, left=0, right=0)), w


def gain(q, w_deg):
    w = w_deg * PC_PER_DEG
    return float(2 * np.sqrt(2 * min(q, 1.0) * RC_PC / w)) if q > 0 else 0.0


def model_stats(M):
    """Noise-free estimator output for a model profile."""
    y = M(X25)
    f = fit_split(X25, y)
    return f['ratio_we'] if f else np.nan, ra_asym(X25, y, f)


# ------------------------------------------------------------------------------------------------ observed
def mask_patch(N):
    d, _ = sep_pa(*HA_PATCH, N['ra'], N['dec'])
    return d < 0.9


def observe(N, seg, mask):
    P = Placed(CC0, N, seg, mask=mask if seg.startswith('N') else None)
    v = N['v'][P.idx]
    st, B = stack_profile(P, v, E25, 'median')
    f = fit_split(X25, st)
    out = dict(primary=dict(fit=f, R_asym=ra_asym(X25, st, f)))
    boots = []
    for _ in range(300):
        idx = RNG.integers(0, P.nblock, P.nblock)
        with np.errstate(all='ignore'):
            y = np.nanmean(B[idx], axis=0)
        fb = fit_split(X25, y)
        boots.append((fb['ratio_we'] if fb else np.nan, ra_asym(X25, y, fb)))
    boots = np.array(boots)
    out['primary']['boot_ratio'] = dict(sd=float(np.nanstd(boots[:, 0])), p16=float(np.nanpercentile(boots[:, 0], 16)),
                                        p84=float(np.nanpercentile(boots[:, 0], 84)))
    out['primary']['boot_Rasym'] = dict(sd=float(np.nanstd(boots[:, 1])), p16=float(np.nanpercentile(boots[:, 1], 16)),
                                        p84=float(np.nanpercentile(boots[:, 1], 84)))
    # variants
    var = {}
    P1 = Placed(CC0, N, seg, mask=mask if seg.startswith('N') else None, block=100.0)
    s1, _ = stack_profile(P1, N['v'][P1.idx], E25, 'median')
    f1 = fit_split(X25, s1); var['single_block_median'] = f1['ratio_we'] if f1 else np.nan
    s2, B2 = stack_profile(P, v, E20, 'median')
    bb = []
    for _ in range(200):
        idx = RNG.integers(0, P.nblock, P.nblock)
        with np.errstate(all='ignore'):
            bb.append(np.nanmean(B2[idx], axis=0))
    e2 = np.nanstd(np.array(bb), axis=0); e2 = np.where(np.isfinite(e2) & (e2 > 0), e2, np.nanmedian(e2))
    f2 = fit_split(X20, s2, w=1 / e2); var['bins0.02_weighted'] = f2['ratio_we'] if f2 else np.nan
    f3 = fit_split(X25, st, xmax=0.6); var['xmax0.6'] = f3['ratio_we'] if f3 else np.nan
    s4, _ = stack_profile(P, v, E25, 'mean'); f4 = fit_split(X25, s4); var['block_mean'] = f4['ratio_we'] if f4 else np.nan
    out['variants'] = var
    return out, st


# ------------------------------------------------------------------------------------------------ nulls
def injection_nulls(N, seg, mask, f_obs):
    A = f_obs['A']; x0 = f_obs['x0']; fw = f_obs['fwhm']
    sg = fw / 2.3548
    models = {'gauss_sym': (lambda xx: np.exp(-0.5 * (xx / sg) ** 2)), }
    meta = {'gauss_sym': dict(q=None, sb=None, w=None, G=None)}
    for sb in SBS:
        for q in QGRID:
            M, w = model(q, sb, fw)
            key = f'q{q}_sb{sb}'
            models[key] = M
            r_nf, ra_nf = model_stats(M)
            meta[key] = dict(q=q, sb=sb, w_deg=w, w_pc=w * PC_PER_DEG, G=gain(q, w), noise_free_ratio=r_nf,
                             noise_free_Rasym=ra_nf)
    vals = {k: [] for k in models}; valsR = {k: [] for k in models}
    used = 0
    for dra, ddec in NSHIFTS:
        cc = shifted_cc(dra, ddec)
        Pn = Placed(cc, N, seg, excl_real=1.0, mask=mask if seg.startswith('N') else None)
        if Pn.idx.size < 500:
            continue
        base = N['v'][Pn.idx]
        st0, _ = stack_profile(Pn, base, E25, 'median')
        if np.isfinite(st0[np.abs(X25) <= 0.9]).mean() < 0.85:
            continue
        used += 1
        for k, M in models.items():
            st, _ = stack_profile(Pn, base + A * M(Pn.x - x0), E25, 'median')
            f = fit_split(X25, st)
            vals[k].append(f['ratio_we'] if f else np.nan)
            valsR[k].append(ra_asym(X25, st, f))
    res = {}
    r_obs = f_obs['ratio_we']
    for k in models:
        a = np.array(vals[k]); a = a[np.isfinite(a)]
        b = np.array(valsR[k]); b = b[np.isfinite(b)]
        mu, sd = float(np.mean(a)), float(np.std(a))
        res[k] = dict(meta[k], n=int(a.size), ratio_mean=mu, ratio_sd=sd, ratio_median=float(np.median(a)),
                      ratio_p16_p84=(float(np.percentile(a, 16)), float(np.percentile(a, 84))),
                      Rasym_median=float(np.median(b)) if b.size else np.nan,
                      p_ge_obs=empirical_p(a, r_obs, 'greater'), p_le_obs=empirical_p(a, r_obs, 'less'),
                      loglike=float(-0.5 * ((r_obs - mu) / sd) ** 2 - np.log(sd)), values=a.tolist())
    return res, used


def main():
    N = load_nsns()
    mask = mask_patch(N)
    res = {}; store = {}
    for seg in SEGS:
        ob, st = observe(N, seg, mask)
        store[seg] = st
        f = ob['primary']['fit']
        nul, used = injection_nulls(N, seg, mask, f)
        r_obs = f['ratio_we']
        # likelihood over q for each blur, normalised to the best
        like = {}
        for sb in SBS:
            ks = [f'q{q}_sb{sb}' for q in QGRID]
            ll = np.array([nul[k]['loglike'] for k in ks])
            like[f'sb{sb}'] = dict(q=list(QGRID), rel_like=np.exp(ll - ll.max()).tolist(),
                                   G=[nul[k]['G'] for k in ks], w_pc=[nul[k]['w_pc'] for k in ks],
                                   ratio_mean=[nul[k]['ratio_mean'] for k in ks], ratio_sd=[nul[k]['ratio_sd'] for k in ks])
        res[seg] = dict(observed=ob, injection=nul, n_placements=used, likelihood_vs_q=like,
                        sym_null_p_ge_obs=nul['gauss_sym']['p_ge_obs'],
                        sym_null_z=float((r_obs - nul['gauss_sym']['ratio_mean']) / nul['gauss_sym']['ratio_sd']))
        print(seg, 'obs r', round(r_obs, 3), ob['primary']['boot_ratio'], 'Rasym', round(ob['primary']['R_asym'], 3),
              'variants', {k: round(v, 3) for k, v in ob['variants'].items()},
              'sym null', round(nul['gauss_sym']['ratio_mean'], 3), round(nul['gauss_sym']['ratio_sd'], 3),
              'p', res[seg]['sym_null_p_ge_obs'], flush=True)
        for sb in SBS:
            print('   ', sb, [(q, round(nul[f'q{q}_sb{sb}']['ratio_mean'], 2), round(nul[f'q{q}_sb{sb}']['G'], 1),
                              round(l, 3)) for q, l in zip(QGRID, like[f'sb{sb}']['rel_like'])], flush=True)
    np.savez_compressed(DER / 'w4_geostats_shape_profiles.npz', x25=X25, **store)
    json_dump(res, OUT)
    print('saved', OUT)


if __name__ == '__main__':
    main()
