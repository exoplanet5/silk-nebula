"""paper_appendix_desi.py -- per-fiber DESI DR1 analysis around the W26 regions P1 and P3 (appendix of the comment).

What the spectra are.  NOT DESI sky fibers (those are hosted at NERSC and were not used).  The inputs are public DESI
DR1 healpix COADDS OF SCIENCE TARGETS (galaxies, QSOs; ZWARN = 0) fetched through NOIRLab SPARCL as +-40 A windows
around the nebular lines (data/derived/spectra_desi_windows, verify_desi_windows; 7,620 coadds).  Each spectrum is
the pipeline sky-subtracted, flux-calibrated flux minus the Redrock best-fit model of the target ("residual").  DESI
subtracts one sky model per petal and exposure, built from that petal's sky fibers, so a Galactic filament survives
as a POSITIVE line residual in fibers on the ridge and as a NEGATIVE one in fibers beside it on the same petal.

This script documents, for the fibers within 0.5 deg of P1 (116.65, +16.60) and P3 (116.45, +14.72):
  (a) the inventory of coadds (counts, programs, target classes, tiles/petals, exposure times, ridge offsets);
  (b) the fit of [O III] 4959+5007 in every individual residual spectrum (the paper's fitter, verify_desi_linefit.
      fit_group: simultaneous quadratic continuum + doublet tied 2.98:1, fixed LSF-based width, chi^2 velocity scan);
  (c) the selection behind Table 2 of the paper and the weighted-mean velocities, their distribution and alternative
      estimators (stacked spectra, other cuts, per tile);
  (d) the line strength: per-fiber surface brightness against ridge offset, the ridge amplitude from the paper's
      forward model (global nodes) and from local fits with the petal-mean term, bootstrap over fibers, and a
      tile-to-tile test of the fiber-acceptance calibration;
  (e) weighted stacks of on-ridge and same-petal off-ridge residual spectra and the line ratios, by two estimators
      (template regression = the paper's family; direct ON-OFF stack), with bootstrap errors and null tests
      (17 line-free windows; random OFF-OFF splits).

All existing functions are imported read-only from verify_desi_*.py; stored per-fiber products of the paper
(verify_desi_linefluxes.npz, verify_desi_coadd_meta.npz) are used for the global quality masks and re-derived here
for the region fibers as a cross-check.

Outputs (never overwrites other files in data/):
  data/derived/paper_appendix_desi_fibers.csv      one row per fiber per region, all fitted quantities
  data/derived/paper_appendix_desi_summary.json    every number the appendix may quote
  data/derived/paper_appendix_desi_stacks.npz      stacked spectra shown in the figure
  paper/figures/paper_appendix_desi.{png,pdf}      (SILK_THEME=light SILK_PAPER=1) or figures/paper_appendix_desi.png

Run: SILK_THEME=light SILK_PAPER=1 ~/.venvs/astro313/bin/python -I scripts/paper_appendix_desi.py
"""
import csv
import json
import sys
import time

import numpy as np

SCR = '/Users/mickey/astronomy/silk-nebula/scripts'
sys.path.insert(0, SCR)
from verify_desi_common import (D, C, WAVE, L, ROOT, NEB_LINES, SKY_LINES, GAL_REST, lsf_sigma,  # noqa: E402
                                load_spectra, load_meta, Exposures, Coadds, Template, petal_subtracted, good_fit,
                                unitvec, vec2radec)
from verify_desi_linefit import fit_group, groups, SIG_INT, PIX  # noqa: E402
from verify_desi_setup import one_R, OMEGA, LAMS  # noqa: E402
from verify_desi_brightness import design, node_list, tri, wls_clip  # noqa: E402
from verify_desi_ratios import iv  # noqa: E402
from silk_geometry import REGIONS  # noqa: E402

# ----------------------------------------------------------------------------- constants of the analysis
REG_NAMES = ('P1', 'P3')
RADIUS = 0.5            # deg, region radius (flat-sky separation, as in verify_desi_velocity.py)
SN_MIN = 4.0            # Table 2 cut: [O III] amplitude / error at the fixed velocity V_NEB
EV_MAX = 30.0           # km/s, maximum formal velocity error accepted
V_FLOOR = 5.0           # km/s, added in quadrature to the formal velocity errors in the weights
T_ON, T_OFF = 0.5, 0.05  # template thresholds of the on-ridge / off-ridge stacks (unit-peak ridge template)
NBOOT = 1000
V_ADOPT, EV_ADOPT = 13.3, 1.0        # adopted [O III] V_helio of the paper (km/s)
W26_OIII = {'P1': (26.8, 4.3), 'P3': (23.8, 3.3)}   # W26 Table 4 as published
RNG = np.random.default_rng(20261118)
LAM5007 = L['OIII5007']
R5008 = one_R(LAM5007)               # 1 R at 5008.24 A in erg s-1 cm-2 arcsec-2
K_FLUX2R = 1e-17 / (OMEGA * R5008)   # R per (1e-17 erg s-1 cm-2) for unit fiber acceptance
OUT_CSV = D / 'paper_appendix_desi_fibers.csv'
OUT_JSON = D / 'paper_appendix_desi_summary.json'
OUT_NPZ = D / 'paper_appendix_desi_stacks.npz'

NULLS = ['null4885', 'null4990', 'null5030', 'null6615', 'null4700', 'null3760', 'null4354', 'null4369.6',
         'null4374.6', 'null4385', 'null4395', 'null4330', 'null7100', 'null7170', 'null7305', 'null7350',
         'null5215']
# name: (fit group, column, chromatic-acceptance key)
LINES = {'OIII5007': ('o3_5007', 0, '5008'), 'OIII4959': ('o3_4959', 0, '5008'), 'Hb': ('hb', 0, '4861'),
         'HeII4686': ('he2', 0, '4686'), 'OII3726': ('o2', 0, '3727'), 'OII3729': ('o2', 1, '3727'),
         'Ha': ('han2', 0, '6563'), 'NII6583': ('han2', 1, '6583'), 'SII6716': ('s2', 0, '6724'),
         'SII6731': ('s2', 1, '6724')}
CK_ALL = ['3727', '3869', '4363', '4686', '4861', '5008', '6563', '7136', '9069']
FIT_GROUPS = ['o3', 'o3_5007', 'o3_4959', 'hb', 'he2', 'o2', 'han2', 's2'] + NULLS


def null_ck(g):
    lam = float(g[4:])
    return min(CK_ALL, key=lambda s: abs(float(s) - lam))


def null_band(lam):
    """Wavelength band of a line or null window (for the wavelength-matched null statistics)."""
    return 'red' if lam > 5800 else ('green' if lam > 4600 else 'blue')


LINE_LAM = {'OIII5007': 5008.2, 'OIII4959': 4960.3, 'Hb': 4862.7, 'HeII4686': 4687.0, 'OII3726': 3727.1,
            'OII3729': 3729.9, 'OII3727': 3728.5, 'Ha': 6564.6, 'NII6583': 6585.3, 'SII6716': 6718.3,
            'SII6731': 6732.7, 'SII6716+6731': 6725.5}


def selection_bias_sim(sig_line_A, v_neb, vscan, n=40000, seed=20261119):
    """Monte Carlo of the Table 2 selection (added by the independent check; own RNG, so no other number of this
    script changes).  A line of true velocity v_true and the fixed model width is added to white noise on the 0.8 A
    grid; the amplitude S/N is measured at the fixed velocity v_neb, the velocity from the chi^2 scan, and fibers
    with S/N > SN_MIN and sigma_V < EV_MAX are averaged with the Table 2 weights.  True S/N: exponential, scale 2.5
    (many faint fibers: the unfavorable case).  Returns the bias of the selected weighted mean and the rms of the
    velocity pulls (1 = formal errors correct for pure noise)."""
    rng = np.random.default_rng(seed)
    lam = LAM5007 + PIX * np.arange(-18, 19)
    vg = np.asarray(vscan, float)
    Gm = np.exp(-0.5 * ((lam[None, :] - LAM5007 * (1 + vg[:, None] / C)) / sig_line_A) ** 2)
    gg = (Gm ** 2).sum(1)
    i0 = int(np.argmin(np.abs(vg - v_neb)))
    dv = float(vg[1] - vg[0])
    out = {}
    for vt in (-6.0, 4.0, 9.0, 14.0, 24.0, 34.0):
        sn_true = rng.exponential(2.5, n)
        g1 = np.exp(-0.5 * ((lam - LAM5007 * (1 + vt / C)) / sig_line_A) ** 2)
        y = (sn_true / np.sqrt(gg[i0]))[:, None] * g1[None, :] + rng.standard_normal((n, lam.size))
        proj = y @ Gm.T
        chi = -proj ** 2 / gg
        sn = proj[:, i0] / np.sqrt(gg[i0])
        k = np.clip(np.argmin(chi, 1), 1, vg.size - 2)
        r_ = np.arange(n)
        a_ = 0.5 * (chi[r_, k - 1] + chi[r_, k + 1] - 2 * chi[r_, k])
        b_ = 0.5 * (chi[r_, k + 1] - chi[r_, k - 1])
        with np.errstate(invalid='ignore', divide='ignore'):
            vf = vg[k] - b_ / (2 * a_) * dv
            ev = np.sqrt(1 / a_) * dv
        sel = (sn > SN_MIN) & np.isfinite(vf) & np.isfinite(ev) & (ev < EV_MAX)
        w = 1 / (ev[sel] ** 2 + V_FLOOR ** 2)
        m_ = float(np.sum(w * vf[sel]) / w.sum())
        out['v_true_%+d' % vt] = dict(n_selected=int(sel.sum()), wmean=m_, err=float(w.sum() ** -0.5),
                                      bias_kms=m_ - vt, pull_rms=float(np.std((vf[sel] - vt) / ev[sel])))
    out['amplitude_loss_20kms_offset'] = float(1 - np.exp(-20.0 ** 2 / (4 * (sig_line_A / LAM5007 * C) ** 2)))
    out['note'] = ('fixed-velocity amplitude of a line displaced by dv: A exp(-dv^2 / 4 sigma^2); the cut pulls the '
                   'selected mean towards v_neb by <= 0.3 km/s for |v_true - v_neb| <= 10 km/s and by about 1 km/s at 20 km/s')
    return out


def jd(o):
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def rnd(x, n=3):
    if x is None:
        return None
    if isinstance(x, (list, tuple, np.ndarray)):
        return [rnd(v, n) for v in x]
    x = float(x)
    return x if not np.isfinite(x) else round(x, n)


def wmean(v, w):
    return float(np.sum(w * v) / np.sum(w))


def mad_sigma(v):
    return float(1.4826 * np.median(np.abs(v - np.median(v))))


# ----------------------------------------------------------------------------- geometry helpers
def ridge_offset(T, ra, dec, strand):
    """Signed great-circle offset (deg, + = east) from a strand of the trace, and the local HWHM (deg)."""
    for s, tree, sig, a, P in T.parts:
        if s != strand:
            continue
        dist, idx = tree.query(unitvec(ra, dec))
        pra, pdec = vec2radec(P[idx])
        sign = np.sign((ra - pra + 180.0) % 360.0 - 180.0)
        return np.degrees(2 * np.arcsin(dist / 2)) * sign, sig[idx] * 1.17741
    raise KeyError(strand)


def basis(ra, dec):
    """Node basis of the paper's forward model (unit-peak ridges x triangular 1-deg Dec nodes) at positions."""
    ra, dec = np.atleast_1d(ra).astype(float), np.atleast_1d(dec).astype(float)
    Tp = {s: Template(strands=(s,))(ra, dec) for s in ('E', 'W', 'W2')}
    return np.array([Tp[s] * tri(dec, j, nodes) for (s, j, nodes) in node_list()]).T


# ----------------------------------------------------------------------------- single-spectrum replica of fit_group
def fit_detail(R1, W1, z1, vc1, cols, v_neb, sigA, v_eval=None, half=30.0):
    """Replica of verify_desi_linefit.fit_group for ONE spectrum that also returns the fitted model (for plots).
    Identical masks, continuum (quadratic over +-30 A), profile and 5-sigma clip; cross-checked against fit_group."""
    cent = np.array([lc for c in cols for (kind, lc), _ in c])
    sel = (WAVE >= cent.min() - half) & (WAVE <= cent.max() + half)
    lam = WAVE[sel]
    r = R1[sel].astype(float)
    w = W1[sel].astype(float)
    bad = ~np.isfinite(r) | ~np.isfinite(w) | (w <= 0)
    r[bad], w[bad] = 0.0, 0.0
    own = set(lc for c in cols for (kind, lc), _ in c)
    m = np.ones(lam.size, bool)
    for l0, hw in [(q, 3.0) for q in NEB_LINES] + [(q, 2.5) for q in SKY_LINES]:
        if any(abs(l0 - o) < 0.5 for o in own):
            continue
        m &= np.abs(lam - l0) > hw
    wm = w * m
    for gl in GAL_REST:
        wm[np.abs(lam - gl * (1 + z1)) < 8.0] = 0.0
    x = (lam - lam.mean()) / half
    sig = np.sqrt(np.interp(lam, WAVE, sigA) ** 2 + (lam * SIG_INT / C) ** 2 + PIX ** 2 / 12)

    def dmat(v):
        G = []
        for c in cols:
            g = np.zeros(lam.size)
            for (kind, l0), a in c:
                mu = l0 * (1 + (vc1 if kind == 'sky' else v) / C)
                s = np.interp(l0, lam, sig)
                g += a * np.exp(-0.5 * ((lam - mu) / s) ** 2) / (np.sqrt(2 * np.pi) * s)
            G.append(g)
        return np.vstack([np.ones(lam.size), x, x * x] + G)

    def solve(v, ww):
        X = dmat(v)
        A = (X * ww) @ X.T
        p = np.linalg.solve(A, (X * ww) @ r)
        return p, X, np.linalg.inv(A)

    p, X, _ = solve(v_neb, wm)
    near = np.zeros(lam.size, bool)
    for l0 in own:
        near |= np.abs(lam - l0) < 3.0
    zz = (r - p @ X) * np.sqrt(wm)
    wm2 = wm * ((np.abs(zz) < 5) | near)
    p0, X0, Ai0 = solve(v_neb, wm2)
    out = dict(lam=lam, r=r, w=wm2, w_raw=w, p_fixed=p0, err_fixed=np.sqrt(np.diag(Ai0)), model_fixed=p0 @ X0,
               cont_fixed=p0[:3] @ X0[:3], sigma_line=float(np.interp(LAM5007, lam, sig)) if lam.min() < LAM5007 <
               lam.max() else float(sig.mean()))
    if v_eval is not None and np.isfinite(v_eval):
        p1, X1, Ai1 = solve(v_eval, wm2)
        out.update(p_free=p1, model_free=p1 @ X1, cont_free=p1[:3] @ X1[:3], err_free=np.sqrt(np.diag(Ai1)))
    return out


# ----------------------------------------------------------------------------- velocity statistics
def vel_stats(v, ev):
    """The paper's estimator (verify_desi_velocity.py) plus distribution diagnostics."""
    n = int(v.size)
    w = 1 / (ev ** 2 + V_FLOOR ** 2)
    vm = wmean(v, w)
    e_formal = float(np.sum(w) ** -0.5)
    e_scatter = float(np.std(v) / np.sqrt(n))
    med = float(np.median(v))
    mad = mad_sigma(v)
    wi = 1 / ev ** 2
    vi = wmean(v, wi)
    chi2_formal = float(np.sum(wi * (v - vi) ** 2) / max(n - 1, 1))
    chi2_floor = float(np.sum(w * (v - vm) ** 2) / max(n - 1, 1))
    # extra scatter s (km/s) that makes the reduced chi^2 about the weighted mean equal to 1
    lo, hi = 0.0, 100.0
    for _ in range(60):
        s = 0.5 * (lo + hi)
        ws = 1 / (ev ** 2 + s ** 2)
        c2 = np.sum(ws * (v - wmean(v, ws)) ** 2) / max(n - 1, 1)
        lo, hi = (s, hi) if c2 > 1 else (lo, s)
    s_extra = 0.5 * (lo + hi)
    ws = 1 / (ev ** 2 + s_extra ** 2)
    bidx = RNG.integers(0, n, (2000, n))
    bw = np.sum(w[bidx] * v[bidx], 1) / np.sum(w[bidx], 1)
    return dict(n=n, wmean=vm, err_formal=e_formal, err_scatter=e_scatter, err=max(e_formal, e_scatter),
                err_bootstrap=float(np.std(bw)), median=med, mad_sigma=mad, err_median=1.2533 * mad / np.sqrt(n),
                std=float(np.std(v)), median_formal_err=float(np.median(ev)), p16=float(np.percentile(v, 16)),
                p84=float(np.percentile(v, 84)), redchi2_formal=chi2_formal, redchi2_with_floor=chi2_floor,
                ivar_mean_no_floor=vi, ivar_err_no_floor=float(np.sum(wi) ** -0.5), extra_scatter_kms=float(s_extra),
                wmean_extra_scatter=wmean(v, ws), err_extra_scatter=float(np.sum(ws) ** -0.5))


# ----------------------------------------------------------------------------- WLS helpers
def wls(y, e, X):
    w = 1 / e ** 2
    A = X.T @ (X * w[:, None])
    Ai = np.linalg.pinv(A)
    b = Ai @ (X.T @ (w * y))
    chi2 = float(np.sum(w * (y - X @ b) ** 2))
    return b, np.sqrt(np.clip(np.diag(Ai), 0, None)), chi2


def wls_clip_boot(y, e, X, nsig=4.0, nboot=NBOOT, n_report=1):
    """WLS with iterative n-sigma clipping (as verify_desi_brightness.wls_clip), then a bootstrap over fibers
    (rows) of the clipped sample.  Returns best fit, formal errors (x sqrt(chi2/dof) if > 1) and bootstrap std."""
    keep = np.ones(len(y), bool)
    for _ in range(5):
        b, eb, _ = wls(y[keep], e[keep], X[keep])
        keep = np.abs((y - X @ b) / e) < nsig
    b, eb, chi2 = wls(y[keep], e[keep], X[keep])
    dof = int(keep.sum() - X.shape[1])
    eb = eb * np.sqrt(max(chi2 / max(dof, 1), 1.0))
    yk, ek, Xk = y[keep], e[keep], X[keep]
    n = len(yk)
    bs = np.empty((nboot, X.shape[1]))
    for k in range(nboot):
        j = RNG.integers(0, n, n)
        bs[k] = wls(yk[j], ek[j], Xk[j])[0]
    return dict(b=b[:n_report], err_formal=eb[:n_report], err_boot=np.std(bs[:, :n_report], 0), chi2=chi2, dof=dof,
                n=int(n), n_clipped=int((~keep).sum()), boot=bs[:, :n_report])


def fe_design(Tcols, glab):
    """Design [template columns | one indicator per tile-petal group] for the petal fixed-effects fit."""
    ug, inv = np.unique(glab, return_inverse=True)
    return np.c_[Tcols, np.eye(ug.size)[inv]], ug


# ----------------------------------------------------------------------------- stacks
class Stacker:
    """Inverse-variance stacks of residual spectra in surface-brightness units.

    For fiber i and pixel lambda:  e_i = r_i * kappa_i,  kappa_i(lambda) = f_i(lambda) * 1e-17 / (Omega * R_5008),
    i.e. energy surface brightness expressed in rayleighs at 5008 A (so that line ratios of stack fluxes are ENERGY
    ratios and [O III] 5007 is in true R).  Weight w_i = ivar_i / kappa_i^2, set to 0 within +-8 A of the target's
    own redshifted lines.  Stack = sum(w e) / sum(w) per pixel."""

    def __init__(self, R, I, z, fnodes, lam_nodes):
        n = R.shape[0]
        kap = np.empty_like(R, dtype=float)
        for i in range(n):
            kap[i] = np.interp(WAVE, lam_nodes, fnodes[i]) * K_FLUX2R
        w = np.where(np.isfinite(I) & (I > 0) & np.isfinite(R), I, 0.0) / kap ** 2
        own = np.zeros(R.shape, bool)
        for gl in GAL_REST:
            own |= np.abs(WAVE[None, :] - gl * (1 + z)[:, None]) < 8.0
        w[own] = 0.0
        w[~np.isfinite(w)] = 0.0
        self.E = np.where(np.isfinite(R), R, 0.0) * kap
        self.W = w
        self.WE = self.W * self.E

    def stack(self, counts):
        """counts: [B, N] resampling counts (or [N] weights 0/1).  Returns stack [B, P] and its ivar [B, P]."""
        c = np.atleast_2d(np.asarray(counts, float))
        den = c @ self.W
        num = c @ self.WE
        with np.errstate(invalid='ignore', divide='ignore'):
            S = np.where(den > 0, num / den, 0.0)
        return S, den

    def onoff(self, on_groups, off_groups, boot=0, split_null=0):
        """Per tile-petal ON minus OFF, combined over petals with per-pixel inverse-variance weights.
        on_groups/off_groups: lists of index arrays (one per tile-petal).  boot > 0: resample fibers with
        replacement within each ON and OFF set.  split_null > 0: instead, split each OFF set at random into two
        halves (pseudo-ON minus pseudo-OFF), a null that contains no ridge signal."""
        n = self.E.shape[0]
        B = max(boot, split_null, 1)
        num = np.zeros((B, WAVE.size))
        den = np.zeros((B, WAVE.size))
        for ion, ioff in zip(on_groups, off_groups):
            con = np.zeros((B, n))
            coff = np.zeros((B, n))
            if split_null:
                for b in range(B):
                    p = RNG.permutation(ioff)
                    h = len(p) // 2
                    con[b, p[:h]] = 1
                    coff[b, p[h:]] = 1
            elif boot:
                for b in range(B):
                    np.add.at(con[b], RNG.choice(ion, ion.size), 1)
                    np.add.at(coff[b], RNG.choice(ioff, ioff.size), 1)
            else:
                con[0, ion] = 1
                coff[0, ioff] = 1
            s1, d1 = self.stack(con)
            s0, d0 = self.stack(coff)
            with np.errstate(invalid='ignore', divide='ignore'):
                iv_ = np.where((d1 > 0) & (d0 > 0), 1 / (1 / d1 + 1 / d0), 0.0)
            num += iv_ * (s1 - s0)
            den += iv_
        with np.errstate(invalid='ignore', divide='ignore'):
            S = np.where(den > 0, num / den, 0.0)
        return S, den


def fit_stack(S, Iv, sigA, v_neb, names, vscan=None):
    """Fit line groups in stacked spectra [B, P] with the paper's fitter (no own-line masks: z = -1)."""
    G = groups(v_neb)
    B = S.shape[0]
    z = np.full(B, -1.0)
    vc = np.zeros(B)
    out = {}
    for g in names:
        out[g] = fit_group(S, Iv, z, vc, G[g], v_neb, sigA, vscan=vscan if g == 'o3' else None)
    return out


def stack_line_fluxes(fs):
    """[O III]5007 (tied doublet), single lines and nulls from a fit_stack result -> dict name: (amp[B], err[B])."""
    out = {'O3': (fs['o3']['amp'][:, 0], fs['o3']['err'][:, 0])}
    for nm, (g, col, ck) in LINES.items():
        out[nm] = (fs[g]['amp'][:, col], fs[g]['err'][:, col])
    for g in NULLS:
        out[g] = (fs[g]['amp'][:, 0], fs[g]['err'][:, 0])
    return out


def ratio_table_from_boot(best, boot, null_names, label):
    """Ratios to [O III]5007 from one best-fit dict and bootstrap dict of line fluxes (energy units)."""
    o3, o3b = best['O3'][0][0], boot['O3'][0]
    tab = {}

    def one(num, numb):
        r = num / o3
        rb = numb / o3b
        p = np.percentile(rb, [2.5, 16, 50, 84, 97.5])
        return dict(ratio=float(r), sig_boot=float(0.5 * (p[3] - p[1])), p2p5=float(p[0]), p16=float(p[1]),
                    p50=float(p[2]), p84=float(p[3]), p97p5=float(p[4]))
    for nm in LINES:
        tab[nm] = one(best[nm][0][0], boot[nm][0])
        tab[nm]['formal_err'] = float(best[nm][1][0] / o3)
    tab['OII3727'] = one(best['OII3726'][0][0] + best['OII3729'][0][0], boot['OII3726'][0] + boot['OII3729'][0])
    tab['SII6716+6731'] = one(best['SII6716'][0][0] + best['SII6731'][0][0],
                              boot['SII6716'][0] + boot['SII6731'][0])
    r5 = boot['OIII5007'][0] / boot['OIII4959'][0]
    tab['5007/4959'] = dict(ratio=float(best['OIII5007'][0][0] / best['OIII4959'][0][0]),
                            p16=float(np.percentile(r5, 16)), p84=float(np.percentile(r5, 84)),
                            sig_boot=float(0.5 * (np.percentile(r5, 84) - np.percentile(r5, 16))))
    r2 = boot['OII3729'][0] / boot['OII3726'][0]
    tab['OII3729/3726'] = dict(ratio=float(best['OII3729'][0][0] / best['OII3726'][0][0]),
                               p16=float(np.percentile(r2, 16)), p84=float(np.percentile(r2, 84)),
                               sig_boot=float(0.5 * (np.percentile(r2, 84) - np.percentile(r2, 16))))
    nul = {}
    for g in null_names:
        nul[g] = one(best[g][0][0], boot[g][0])
    vals = np.array([nul[g]['ratio'] for g in null_names])
    sigs = np.array([nul[g]['sig_boot'] for g in null_names])
    tab['_nulls'] = dict(lines=nul, rms_ratio=float(np.sqrt(np.mean(vals ** 2))), max_abs=float(np.max(np.abs(vals))),
                         mean_sig_boot=float(np.mean(sigs)),
                         chi_rms=float(np.sqrt(np.mean((vals / sigs) ** 2))), n=len(null_names))
    tab['_label'] = label
    return tab


def add_limits(tab, extra=None):
    """Uncertainty adopted for each ratio and 2-sigma limits.

    sig_infl  = bootstrap sigma x max(chi_rms of the 17 null windows, 1)   (bootstrap calibrated on the nulls)
    sig_adopt = max(sig_infl, rms of the null-window ratios[, OFF-OFF split-null sigma])   (no double counting:
                the null rms already contains the statistical noise)
    upper_2sig = max(ratio, 0) + 2 sig_adopt;  snr = ratio / sig_adopt."""
    infl = max(tab['_nulls']['chi_rms'], 1.0)
    floor = tab['_nulls']['rms_ratio']
    tab['_inflation'] = float(infl)
    # The 17 null windows are heterogeneous: the blue ones (3760-4395 A) are 3-6 times noisier than the red ones, so
    # the single floor above is conservative for the red lines.  Diagnostic added by the independent check: the same
    # recipe with the null windows of the line's own wavelength band (sig_matched, snr_matched).  The adopted
    # uncertainties and limits are unchanged.
    by_band = {}
    for b in ('blue', 'green', 'red'):
        nb_ = [g for g in tab['_nulls']['lines'] if null_band(float(g[4:])) == b]
        v_ = np.array([tab['_nulls']['lines'][g]['ratio'] for g in nb_])
        s_ = np.array([tab['_nulls']['lines'][g]['sig_boot'] for g in nb_])
        by_band[b] = dict(n=len(nb_), windows=[g[4:] for g in nb_], rms_ratio=float(np.sqrt(np.mean(v_ ** 2))),
                          max_abs=float(np.max(np.abs(v_))), chi_rms=float(np.sqrt(np.mean((v_ / s_) ** 2))))
    tab['_nulls']['by_band'] = by_band
    tab['_nulls']['by_band_note'] = ('blue: null windows < 4600 A; green: 4600-5800 A; red: > 5800 A.  sig_matched = '
                                     'max(sig_boot x max(chi_rms of the band, 1), rms of the band[, split-null sigma])')
    for nm, d in tab.items():
        if nm.startswith('_') or 'sig_boot' not in d:
            continue
        d['sig_infl'] = float(d['sig_boot'] * infl)
        if nm in ('5007/4959', 'OII3729/3726'):
            continue
        ex_ = extra.get(nm, 0.0) if extra else 0.0
        d['sig_adopt'] = float(max(d['sig_infl'], floor, ex_))
        d['snr'] = float(d['ratio'] / d['sig_adopt'])
        d['snr_boot_only'] = float(d['ratio'] / d['sig_infl'])
        d['upper_2sig'] = float(max(d['ratio'], 0.0) + 2 * d['sig_adopt'])
        if nm in LINE_LAM:
            bb_ = by_band[null_band(LINE_LAM[nm])]
            d['null_band'] = null_band(LINE_LAM[nm])
            d['sig_matched'] = float(max(d['sig_boot'] * max(bb_['chi_rms'], 1.0), bb_['rms_ratio'], ex_))
            d['snr_matched'] = float(d['ratio'] / d['sig_matched'])
    hb = tab['Hb']
    tab['OIII5007/Hb'] = dict(lower_2sig=float(1 / hb['upper_2sig']),
                              central=float(1 / hb['ratio']) if hb['ratio'] > 0 else None,
                              note='lower limit = 1 / (max(Hb/[O III], 0) + 2 sig_adopt)')
    return tab


# ----------------------------------------------------------------------------- main analysis
def main():
    t0 = time.time()
    lf = np.load(D / 'verify_desi_linefluxes.npz')
    cm = np.load(D / 'verify_desi_coadd_meta.npz')
    det = json.loads((D / 'verify_desi_detect.json').read_text())
    sc = det['err_scale_from_nulls']
    V_NEB = float(lf['v_neb'])
    ids_lf = lf['sparcl_id']
    ids, F, I, M, src = load_spectra()
    pos = {s: k for k, s in enumerate(ids)}
    kk = np.array([pos[s] for s in ids_lf])
    Rall = (F - M)[kk]
    Iall = I[kk]
    Fall = F[kk]
    Mall = M[kk]
    del F, M, I
    keep, meta = load_meta(ids_lf)
    assert keep.all()
    sigA, n_lsf = lsf_sigma()
    N = len(ids_lf)
    ra, dec = meta['ra'], meta['dec']
    print('loaded', N, 'coadds  %.0fs' % (time.time() - t0), flush=True)

    # stored per-fiber products of the paper (global masks need the whole sample)
    O_st, EO_st = lf['o3__amp'][:, 0], lf['o3__err'][:, 0] * sc
    kR = cm['kR_5008']
    Zf = cm['Z'] / kR
    base = good_fit(lf, 'o3', scale=sc) & np.isfinite(Zf) & cm['ok_contrib']
    snr_st = O_st / EO_st
    glab = np.array([f'{t}_{p}' for t, p in zip(cm['group_tile'], cm['group_petal'])])
    tile = cm['group_tile']

    # geometry: template parts, ridge offsets
    T = Template()
    offE, hwE = ridge_offset(T, ra, dec, 'E')
    offW, hwW = ridge_offset(T, ra, dec, 'W')
    TE, TW, Ttot = cm['T_E'], cm['T_W'], cm['T']

    # ---------------- global forward model of the paper (verify_desi_brightness.py), reproduced
    ex = Exposures()
    co = Coadds(meta, ex)
    Zre = petal_subtracted(Template(), meta, co)
    okz = np.isfinite(Zre) & np.isfinite(cm['Z'])
    X, names = design(meta, co)
    okb = np.isfinite(O_st) & np.isfinite(EO_st) & (EO_st > 0) & cm['ok_contrib'] & np.isfinite(kR) & good_fit(lf, 'o3')
    yR, eR = O_st * kR, EO_st * kR
    bnode, ebnode, chi2g, dofg, keepg, act = wls_clip(yR[okb], eR[okb], X[okb])
    bnode0 = np.where(np.isfinite(bnode), bnode, 0.0)
    good_nodes = [(n, a, s) for n, a, s, ac in zip(names, bnode, ebnode, act) if ac and s < 1.0]
    south = [a for n, a, s in good_nodes if float(n.split('@')[1]) < 19.5 and not n.startswith('W2')]
    Bfib = basis(ra, dec)
    model_resid = X @ bnode0                 # predicted petal-relative residual (R)
    model_sky = Bfib @ bnode0                # predicted sky brightness at the fiber (R)
    S_glob = model_sky - model_resid         # global-model tile-petal mean term (R)
    # bootstrap over fibers of the global fit (fixed clip mask) for the model values at the regions
    Xg, yg, eg = X[okb][keepg], yR[okb][keepg], eR[okb][keepg]
    actc = np.where(act)[0]
    nb_glob = 300
    bboot = np.zeros((nb_glob, X.shape[1]))
    for k in range(nb_glob):
        j = RNG.integers(0, len(yg), len(yg))
        bboot[k, actc] = wls(yg[j], eg[j], Xg[j][:, actc])[0]
    print('global forward model: south median %.3f R, chi2/dof %.0f/%d  %.0fs' % (np.median(south), chi2g, dofg,
                                                                              time.time() - t0), flush=True)

    # acceptance bookkeeping
    fnodes_keys = sorted(LAMS, key=lambda k_: LAMS[k_])
    lam_nodes = np.array([LAMS[k_] for k_ in fnodes_keys])
    fnodes_all = np.array([cm[f'f_chrom_{k_}'] for k_ in fnodes_keys]).T
    filled = cm['frac_ff_filled'] > 0.5
    # tile-level acceptance scatter among tiles WITH GFA measurements (for the uncertainty of filled tiles)
    tile_f = {}
    for t_ in np.unique(tile[cm['ok_contrib']]):
        m_ = (tile == t_) & cm['ok_contrib'] & (cm['n_tiles'] == 1)
        if m_.any() and not filled[m_][0]:
            tile_f[int(t_)] = (float(np.median(cm['f_chrom_5008'][m_])), str(meta['program'][m_][0]))
    facc_prog = {}
    for p_ in ('dark', 'bright'):
        v_ = np.array([f for f, pp in tile_f.values() if pp == p_])
        if v_.size:
            facc_prog[p_] = dict(n_tiles=int(v_.size), median=float(np.median(v_)), p16=float(np.percentile(v_, 16)),
                                 p84=float(np.percentile(v_, 84)), min=float(v_.min()), max=float(v_.max()))

    summary = {
        'meta': {
            'script': 'scripts/paper_appendix_desi.py', 'date': '2026-11-18',
            'data': 'public DESI DR1 healpix coadds of science targets (GALAXY/QSO, ZWARN=0) via NOIRLab SPARCL; '
                    '+-40 A line windows; residual = coadd flux - Redrock model. No DESI sky-fiber spectra used.',
            'n_coadds_total': int(N), 'n_usable_o3_global': int(base.sum()),
            'region_radius_deg': RADIUS, 'regions': {k: REGIONS[k] for k in REG_NAMES},
        },
        'constants': {
            'c_kms': C, 'lambda_vac_A': {k: rnd(L[k], 3) for k in ('OIII5007', 'OIII4959', 'Hb', 'HeII4686',
                                                                 'OII3726', 'OII3729', 'Ha', 'NII6583', 'NII6548',
                                                                 'SII6716', 'SII6731', 'sky5577', 'HgI4358')},
            'air_to_vacuum': 'Edlen 1966 / IAU (Morton 1991)',
            'pixel_A': PIX, 'wave_grid': '3600 + 0.8 i A (vacuum, barycentric)',
            'lsf_sigma_A_5008': rnd(float(np.interp(LAM5007, WAVE, sigA)), 4), 'n_lsf_spectra': int(n_lsf),
            'sigma_intrinsic_kms': SIG_INT,
            'sigma_line_A_5008': rnd(float(np.sqrt(np.interp(LAM5007, WAVE, sigA) ** 2 + (LAM5007 * SIG_INT / C) ** 2
                                                   + PIX ** 2 / 12)), 4),
            'doublet_ratio_5007_4959_tied': 2.98, 'v_neb_fixed_kms': V_NEB,
            'vscan_kms': [float(lf['vscan'][0]), float(lf['vscan'][-1]), float(lf['vscan'][1] - lf['vscan'][0])],
            'continuum': 'quadratic over [lambda_min - 30, lambda_max + 30] A, fitted simultaneously',
            'masks': 'other nebular lines +-3 A; sky lines (Hg I 4047/4358/5461/5770/5791, Na D, [O I]5577) +-2.5 A; '
                     'target lines at (1+z) +-8 A; one 5-sigma clip outside +-3 A of the fitted lines',
            'err_scale_from_nulls': sc,
            'fiber_diameter_arcsec': 1.52, 'fiber_solid_angle_arcsec2': rnd(OMEGA, 4),
            'one_R_5008_erg_s_cm2_arcsec2': float(R5008),
            'one_R_definition': '1e6/(4 pi) photons cm-2 s-1 sr-1',
            'R_per_1e-17_flux_unit_acceptance': rnd(K_FLUX2R, 4),
            'acceptance_median_all_5008': rnd(float(np.nanmedian(cm['f_chrom_5008'])), 3),
            'kR_median_all': rnd(float(np.nanmedian(kR)), 4),
            'table2_selection': f'within {RADIUS} deg; good fit; S/N > {SN_MIN} at fixed V_NEB (positive residual); '
                                f'finite v_free; ev_free < {EV_MAX} km/s',
            'table2_weights': f'1/(ev_free^2 + {V_FLOOR}^2 km2/s2); error = max(sum(w)^-1/2, std/sqrt(n))',
            'stack_T_on': T_ON, 'stack_T_off': T_OFF, 'n_boot': NBOOT,
            'v_adopted_paper': [V_ADOPT, EV_ADOPT],
        },
        'global_forward_model': {
            'south_median_R': rnd(float(np.median(south)), 3), 'chi2': rnd(chi2g, 1), 'dof': int(dofg),
            'n_fibers': int(okb.sum()),
            'nodes_near_regions': {n: [rnd(a, 3), rnd(s, 3)] for n, a, s in good_nodes
                                   if n.split('@')[0] in ('E', 'W') and 13.5 < float(n.split('@')[1]) < 18.5},
        },
        'checks': {
            'template_Z_recomputed_max_abs_diff': rnd(float(np.max(np.abs(Zre[okz] - cm['Z'][okz]))), 6),
        },
        'acceptance': {
            'tiles_with_gfa_by_program': {k: {kk_: rnd(v, 3) for kk_, v in d.items()} for k, d in facc_prog.items()},
            'frac_filled_all_usable': rnd(float(np.mean(filled[base])), 3),
            'frac_filled_usable_dec_lt_19.5': rnd(float(np.mean(filled[base & (dec < 19.5)])), 3),
            'filled_tiles': {str(int(t_)): int(c_) for t_, c_ in zip(*np.unique(tile[base & filled],
                                                                                 return_counts=True))},
            'note': 'filled = exposure table has no GFA fiber fraction (0); the program median was used',
        },
        'regions': {},
    }
    vj = json.loads((D / 'verify_desi_velocity.json').read_text())['sky_frame_test']
    rj = json.loads((D / 'spectra_desi_rvcheck.json').read_text())['three_way']
    lg = rj['LRS_minus_Gaia']['median']
    summary['zero_point_global'] = {
        'frame_test_OI5577_slope_vs_vcorr': [vj['OI5577']['slope'], vj['OI5577']['err_slope']],
        'frame_test_estimator_response_sim': vj['OI5577']['sim_slope_stacked_width'],
        'frame_test_calibrated_slope': vj['OI5577']['slope'] / vj['OI5577']['sim_slope_stacked_width'],
        'frame_test_n_coadds': vj['OI5577']['n'],
        'stars_DESI_sv1_bright_minus_Gaia': [rj['DESI_sv1/bright_minus_LRS']['median'] + lg,
                                             rj['DESI_sv1/bright_minus_LRS']['err_median'],
                                             rj['DESI_sv1/bright_minus_LRS']['n']],
        'stars_DESI_sv1_dark_minus_Gaia': [rj['DESI_sv1/dark_minus_LRS']['median'] + lg,
                                           rj['DESI_sv1/dark_minus_LRS']['err_median'],
                                           rj['DESI_sv1/dark_minus_LRS']['n']],
        'stars_DESI_main_bright_minus_Gaia': [rj['DESI_main/bright_minus_LRS']['median'] + lg,
                                              rj['DESI_main/bright_minus_LRS']['err_median'],
                                              rj['DESI_main/bright_minus_LRS']['n']],
        'stars_DESI_main_dark_minus_Gaia': [rj['DESI_main/dark_minus_LRS']['median'] + lg,
                                            rj['DESI_main/dark_minus_LRS']['err_median'],
                                            rj['DESI_main/dark_minus_LRS']['n']],
        'stars_note': 'DESI MWS stellar velocities minus LAMOST LRS, plus (LRS - Gaia DR3) = %.2f km/s; stars in the '
                      'same sky box (spectra_desi_rvcheck.json)' % lg,
        'zero_point_uncertainty_adopted_kms': 1.5,
        'note': 'values quoted from verify_desi_velocity.json and spectra_desi_rvcheck.json (not recomputed here)',
    }

    summary['selection_bias_simulation'] = selection_bias_sim(
        float(np.sqrt(np.interp(LAM5007, WAVE, sigA) ** 2 + (LAM5007 * SIG_INT / C) ** 2 + PIX ** 2 / 12)), V_NEB,
        lf['vscan'])

    # catalog-level inventory (all DR1 coadds listed by SPARCL in the box)
    cat_rows = []
    with open(D / 'spectra_desi_catalog.csv') as fh:
        for r in csv.DictReader(fh):
            d_ = float(r['dec'])
            if 13.9 < d_ < 17.3:
                cat_rows.append((float(r['ra']), d_, r['spectype'], int(r['zwarn']), r['sparcl_id']))
    cra = np.array([c[0] for c in cat_rows])
    cdec = np.array([c[1] for c in cat_rows])
    ctyp = np.array([c[2] for c in cat_rows])
    czw = np.array([c[3] for c in cat_rows])
    cid = np.array([c[4] for c in cat_rows])
    have = set(ids_lf.tolist())

    fiber_rows = []
    figdata = {}
    G = groups(V_NEB)
    vscan = lf['vscan']
    stk_out = {'wave': WAVE}

    for reg in REG_NAMES:
        x0, y0 = REGIONS[reg]
        sep = np.hypot((ra - x0) * np.cos(np.radians(y0)), dec - y0)
        inreg = sep < RADIUS
        idx = np.where(inreg)[0]
        nreg = idx.size
        res = {}
        # ---------------------------------------------------------------- (a) inventory
        csep = np.hypot((cra - x0) * np.cos(np.radians(y0)), cdec - y0) < RADIUS
        sel_cat = csep & np.isin(ctyp, ['GALAXY', 'QSO']) & (czw == 0)
        inv = {
            'dr1_coadds_in_sparcl_catalog': int(csep.sum()),
            'by_spectype_catalog': {k: int(v) for k, v in zip(*np.unique(ctyp[csep], return_counts=True))},
            'galaxy_qso_zwarn0_catalog': int(sel_cat.sum()),
            'fetched_and_analyzed': int(nreg),
            'fetched_fraction_of_selection': rnd(float(np.mean([c in have for c in cid[sel_cat]])), 3),
            'usable_good_o3_fit': int((inreg & base).sum()),
            'by_program': {k: int(v) for k, v in zip(*np.unique(meta['program'][idx], return_counts=True))},
            'by_survey': {k: int(v) for k, v in zip(*np.unique(meta['survey'][idx], return_counts=True))},
            'by_spectype': {k: int(v) for k, v in zip(*np.unique(meta['spectype'][idx], return_counts=True))},
            'redshift_range': rnd([float(np.min(meta['redshift'][idx])), float(np.median(meta['redshift'][idx])),
                                   float(np.max(meta['redshift'][idx]))], 3),
            'coadd_exptime_s_min_med_max': rnd(np.nanpercentile(meta['coadd_exptime'][idx], [0, 50, 100]), 0),
            'efftime_s_min_med_max': rnd(np.nanpercentile(cm['efftime'][idx], [0, 50, 100]), 0),
            'acceptance_5008_min_med_max': rnd(np.nanpercentile(cm['f_chrom_5008'][idx], [0, 50, 100]), 3),
            'frac_acceptance_filled': rnd(float(np.mean(filled[idx])), 3),
            'vcorr_kms_min_med_max': rnd(np.nanpercentile(cm['vcorr'][idx], [0, 50, 100]), 2),
            'sep_deg_min_med_max': rnd(np.percentile(sep[idx], [0, 50, 100]), 3),
            'offset_E_deg_min_med_max': rnd(np.percentile(offE[idx], [0, 50, 100]), 3),
            'abs_offset_nearest_strand_deg_med': rnd(float(np.median(np.minimum(np.abs(offE[idx]),
                                                                               np.abs(offW[idx])))), 3),
            'n_T_ge_0.5': int((Ttot[idx] >= T_ON).sum()), 'n_T_lt_0.05': int((Ttot[idx] < T_OFF).sum()),
            'template_at_region': {s: rnd(float(v[0]), 3) for s, v in
                                   T.eval_parts(np.array([x0]), np.array([y0])).items()},
            'E_ridge_hwhm_deg_at_region': rnd(float(ridge_offset(T, np.array([x0]), np.array([y0]), 'E')[1][0]), 3),
            'region_offset_from_E_deg': rnd(float(ridge_offset(T, np.array([x0]), np.array([y0]), 'E')[0][0]), 3),
            'tile_petals': {},
        }
        for g_ in np.unique(glab[idx]):
            q = idx[glab[idx] == g_]
            t_ = int(tile[q][0])
            js = ex.by_tile.get(t_, [])
            nights = sorted({int(float(ex.rows[j]['night'])) for j in js})
            inv['tile_petals'][g_] = dict(
                n=int(q.size), n_good=int(base[q].sum()), survey=str(meta['survey'][q][0]),
                program=str(meta['program'][q][0]), n_exp=len(js), nights=[nights[0], nights[-1]] if nights else [],
                n_nights=len(nights), exptime_s=rnd(float(np.nanmedian(meta['coadd_exptime'][q])), 0),
                efftime_s=rnd(float(np.nanmedian(cm['efftime'][q])), 0),
                f_acc_5008=rnd(float(np.nanmedian(cm['f_chrom_5008'][q])), 3), acceptance_filled=bool(filled[q][0]),
                vcorr_kms=rnd(float(np.nanmedian(cm['vcorr'][q])), 2),
                vcorr_spread_kms=rnd(float(np.nanmedian(cm['vspread'][q])), 1),
                petal_mean_template=rnd(float(np.nanmedian((cm['T'] - cm['Z'])[q])), 3),
                n_on=int((Ttot[q] >= T_ON).sum()), n_off=int((Ttot[q] < T_OFF).sum()))
        res['inventory'] = inv

        # ---------------------------------------------------------------- (b) per-fiber fits (re-run)
        Rr, Ir, zr = Rall[idx], Iall[idx], meta['redshift'][idx]
        vc = np.where(np.isfinite(cm['vcorr'][idx]), cm['vcorr'][idx], 0.0)
        fits = {}
        for g_ in FIT_GROUPS:
            fits[g_] = fit_group(Rr, Ir, zr, vc, G[g_], V_NEB, sigA, vscan=vscan if g_ == 'o3' else None)
        O = fits['o3']['amp'][:, 0]
        EO = fits['o3']['err'][:, 0] * sc
        vfree, evfree = fits['o3']['v_free'], fits['o3']['ev_free']
        dchi = fits['o3']['dchi_free']
        snr = O / EO
        chk = dict(max_abs_dflux=float(np.nanmax(np.abs(O - O_st[idx]))),
                   max_abs_dv=float(np.nanmax(np.abs(vfree - lf['o3__v_free'][idx]))))
        res['refit_vs_stored'] = {k: rnd(v, 6) for k, v in chk.items()}
        kRr = kR[idx]
        SB, eSB = O * kRr, EO * kRr
        gb = base[idx]
        gl_r = glab[idx]
        sig_line_A = float(np.sqrt(np.interp(LAM5007, WAVE, sigA) ** 2 + (LAM5007 * SIG_INT / C) ** 2
                                   + PIX ** 2 / 12))
        res['per_fiber'] = {
            'n': int(nreg), 'n_good': int(gb.sum()),
            'n_pos_sn_gt4': int((gb & (snr > SN_MIN)).sum()), 'n_neg_sn_lt_m4': int((gb & (snr < -SN_MIN)).sum()),
            'n_abs_sn_lt4': int((gb & (np.abs(snr) <= SN_MIN)).sum()),
            'flux_err_median_1e-17': rnd(float(np.median(EO[gb])), 3),
            'SB_err_median_R': rnd(float(np.median(eSB[gb])), 3),
            'SB_range_R': rnd([float(SB[gb].min()), float(SB[gb].max())], 2),
            'sn_max': rnd(float(snr[gb].max()), 1),
            'rchi2_median': rnd(float(np.median(fits['o3']['rchi2'][gb])), 3),
            'sigma_line_A': rnd(sig_line_A, 4), 'sigma_line_kms': rnd(sig_line_A / LAM5007 * C, 2),
        }

        # ---------------------------------------------------------------- (c) Table 2 velocity
        sel = gb & (snr > SN_MIN) & np.isfinite(vfree) & np.isfinite(evfree) & (evfree < EV_MAX)
        vel = {'table2': vel_stats(vfree[sel], evfree[sel])}
        vel['table2']['paper_value'] = {'P1': [11.9, 1.9, 61], 'P3': [9.2, 1.0, 128]}[reg]
        vel['table2']['n_excluded_ev_ge_30'] = int((gb & (snr > SN_MIN) & ~((evfree < EV_MAX) & np.isfinite(vfree))).sum())
        vel['table2']['sn_median'] = float(np.median(snr[sel]))
        vel['table2']['frac_T_ge_0.5'] = float(np.mean(Ttot[idx][sel] >= T_ON))
        vel['table2']['V_LSR_offset_note'] = 'V_LSR = V_helio - 12.10 (P1), - 12.58 (P3) km/s (paper Table 2 comments)'
        variants = {}

        def var(lab, m_):
            m_ = m_ & np.isfinite(vfree) & np.isfinite(evfree) & (evfree < EV_MAX)
            if m_.sum() >= 5:
                s_ = vel_stats(vfree[m_], evfree[m_])
                variants[lab] = {k: s_[k] for k in ('n', 'wmean', 'err', 'err_formal', 'err_scatter', 'median',
                                                    'mad_sigma')}
        for cut in (3, 5, 6, 8):
            var(f'pos_sn_gt{cut}', gb & (snr > cut))
        var('pos_sn_gt4_within0.3deg', gb & (snr > SN_MIN) & (sep[idx] < 0.3))
        var('pos_sn_gt4_T_ge_0.5', gb & (snr > SN_MIN) & (Ttot[idx] >= T_ON))
        var('pos_sn_gt4_E_only', gb & (snr > SN_MIN) & (TE[idx] > 0.3) & (TW[idx] < 0.1))
        var('pos_sn_gt4_W_only', gb & (snr > SN_MIN) & (TW[idx] > 0.3) & (TE[idx] < 0.1))
        var('neg_sn_lt_m4', gb & (snr < -SN_MIN))
        for t_ in np.unique(tile[idx]):
            var(f'pos_sn_gt4_tile{int(t_)}', gb & (snr > SN_MIN) & (tile[idx] == t_))
        var('pos_sn_gt4_dark', gb & (snr > SN_MIN) & (meta['program'][idx] == 'dark'))
        var('pos_sn_gt4_bright', gb & (snr > SN_MIN) & (meta['program'][idx] == 'bright'))
        var('pos_sn_gt4_QSO_or_z_gt_0.5', gb & (snr > SN_MIN) & ((meta['spectype'][idx] == 'QSO') | (zr > 0.5)))
        for lab, m_ in (('pos_sn_gt4_east_of_E_spine', sel & (offE[idx] >= 0)),
                        ('pos_sn_gt4_west_of_E_spine', sel & (offE[idx] < 0))):
            if m_.sum() >= 5:       # same estimator, computed directly (no bootstrap draw)
                w_ = 1 / (evfree[m_] ** 2 + V_FLOOR ** 2)
                variants[lab] = dict(n=int(m_.sum()), wmean=wmean(vfree[m_], w_),
                                     err=float(max(w_.sum() ** -0.5, np.std(vfree[m_]) / np.sqrt(m_.sum()))),
                                     err_formal=float(w_.sum() ** -0.5),
                                     err_scatter=float(np.std(vfree[m_]) / np.sqrt(m_.sum())),
                                     median=float(np.median(vfree[m_])), mad_sigma=mad_sigma(vfree[m_]))
        vel['variants'] = variants
        # trend of v with S/N and with ridge offset (selected fibers)
        if sel.sum() > 10:
            ww = 1 / (evfree[sel] ** 2 + V_FLOOR ** 2)
            for lab, xx in (('offset_E_deg', offE[idx][sel]), ('log10_sn', np.log10(snr[sel]))):
                Xd = np.c_[np.ones(sel.sum()), xx - np.mean(xx)]
                b_, e_, c2 = wls(vfree[sel], 1 / np.sqrt(ww), Xd)
                e_ = e_ * np.sqrt(max(c2 / (sel.sum() - 2), 1))
                vel[f'slope_vs_{lab}'] = [float(b_[1]), float(e_[1])]

        # stacks for the velocity (alternative estimator)
        stk = Stacker(Rr, Ir, zr, fnodes_all[idx], lam_nodes)
        on_T = gb & (Ttot[idx] >= T_ON)

        def stack_velocity(mask, lab, nb=400):
            c1 = mask.astype(float)[None, :]
            S, Iv = stk.stack(c1)
            f1 = fit_stack(S, Iv, sigA, V_NEB, ['o3'], vscan=vscan)['o3']
            ii = np.where(mask)[0]
            cb = np.zeros((nb, nreg))
            for b in range(nb):
                np.add.at(cb[b], RNG.choice(ii, ii.size), 1)
            Sb, Ib = stk.stack(cb)
            fb = fit_stack(Sb, Ib, sigA, V_NEB, ['o3'], vscan=vscan)['o3']
            vb = fb['v_free'][np.isfinite(fb['v_free'])]
            return dict(n=int(mask.sum()), v=float(f1['v_free'][0]), err_formal=float(f1['ev_free'][0]),
                        err_boot=float(np.std(vb)), p16=float(np.percentile(vb, 16)), p84=float(np.percentile(vb, 84)),
                        flux_R=float(f1['amp_vfree'][0]), flux_R_fixed_v=float(f1['amp'][0, 0]),
                        err_flux_R=float(f1['err'][0, 0])), S[0], Iv[0]
        vel['stack_table2_fibers'], S_sel, I_sel = stack_velocity(sel, 'sel')
        vel['stack_on_ridge_T_ge_0.5'], S_on, I_on = stack_velocity(on_T, 'on')
        res['velocity'] = vel

        # ---------------------------------------------------------------- (d) line strength
        strength = {}
        y, e = SB, eSB
        okr = gb & np.isfinite(y) & np.isfinite(e) & (e > 0)
        # model-free weighted means of the residual near the region centre (lower limits: petal-relative)
        for rad in (0.1, 0.2):
            m_ = okr & (sep[idx] < rad)
            if m_.sum() >= 3:
                w_ = 1 / e[m_] ** 2
                strength[f'resid_mean_R_r{rad}'] = [wmean(y[m_], w_), float(w_.sum() ** -0.5), int(m_.sum())]
        m_ = okr & (Ttot[idx] >= T_ON)
        w_ = 1 / e[m_] ** 2
        strength['resid_mean_R_on_ridge'] = [wmean(y[m_], w_), float(w_.sum() ** -0.5), int(m_.sum())]
        strength['resid_median_R_on_ridge'] = float(np.median(y[m_]))
        strength['resid_max_R'] = [float(np.max(y[okr])), float(e[okr][np.argmax(y[okr])])]
        m_ = okr & (Ttot[idx] < T_OFF)
        if m_.sum() >= 3:
            w_ = 1 / e[m_] ** 2
            strength['resid_mean_R_off_ridge'] = [wmean(y[m_], w_), float(w_.sum() ** -0.5), int(m_.sum())]
        # global forward model (the paper's) at the region and on the ridge at the region's Dec
        Breg = basis(np.array([x0]), np.array([y0]))[0]
        nlist = node_list()
        Erow = np.array([tri(np.array([y0]), j, nodes)[0] if s == 'E' else 0.0 for (s, j, nodes) in nlist])
        Wrow = np.array([tri(np.array([y0]), j, nodes)[0] if s == 'W' and nodes[0] <= y0 <= nodes[-1] else 0.0
                         for (s, j, nodes) in nlist])
        strength['global_model'] = dict(
            sky_SB_at_region_R=[float(Breg @ bnode0), float(np.std(bboot @ Breg))],
            E_ridge_peak_at_dec_R=[float(Erow @ bnode0), float(np.std(bboot @ Erow))],
            W_ridge_peak_at_dec_R=[float(Wrow @ bnode0), float(np.std(bboot @ Wrow))],
            petal_term_median_R=float(np.median(S_glob[idx][okr])),
            note='errors: bootstrap over fibers of the global fit (clip mask fixed); calibration +-8% not included')
        # local fits
        TEr, TWr, Tr = TE[idx], TW[idx], Ttot[idx]
        two = reg == 'P1'           # E and W strands are separated near P1 (0.27 deg); degenerate near P3
        Tcols = np.c_[TEr, TWr] if two else Tr[:, None]
        ncol = Tcols.shape[1]
        loc = {}
        # B: petal term fixed to the global model
        fB = wls_clip_boot(y[okr] + S_glob[idx][okr], e[okr], Tcols[okr], n_report=ncol)
        loc['B_global_petal_term'] = dict(amp_R=fB['b'], err_boot=fB['err_boot'], err_formal=fB['err_formal'],
                                          chi2=fB['chi2'], dof=fB['dof'], n=fB['n'], n_clipped=fB['n_clipped'])
        # C: free constant per tile-petal (fixed effects)
        gcount = {g_: int((gl_r[okr] == g_).sum()) for g_ in np.unique(gl_r[okr])}
        okc = okr & np.array([gcount.get(g_, 0) >= 5 for g_ in gl_r])
        Xc, ug = fe_design(Tcols[okc], gl_r[okc])
        fC = wls_clip_boot(y[okc], e[okc], Xc, n_report=Xc.shape[1])
        loc['C_petal_fixed_effects'] = dict(amp_R=fC['b'][:ncol], err_boot=fC['err_boot'][:ncol],
                                            err_formal=fC['err_formal'][:ncol], chi2=fC['chi2'], dof=fC['dof'],
                                            n=fC['n'], n_clipped=fC['n_clipped'],
                                            petal_constants_R={g_: [float(fC['b'][ncol + k]),
                                                                    float(fC['err_boot'][ncol + k])]
                                                               for k, g_ in enumerate(ug)})
        # D: single amplitude on the petal-subtracted template Z (uniform amplitude over the petals)
        Zr = cm['Z'][idx]
        fD = wls_clip_boot(y[okr], e[okr], Zr[okr][:, None])
        loc['D_single_amp_on_Z'] = dict(amp_R=fD['b'], err_boot=fD['err_boot'], err_formal=fD['err_formal'],
                                        chi2=fD['chi2'], dof=fD['dof'], n=fD['n'], n_clipped=fD['n_clipped'])
        # single-amplitude version of C (E+W together) for both regions
        Xc1, ug1 = fe_design(Tr[okc][:, None], gl_r[okc])
        fC1 = wls_clip_boot(y[okc], e[okc], Xc1, n_report=1)
        loc['C1_petal_fixed_effects_single_amp'] = dict(amp_R=fC1['b'], err_boot=fC1['err_boot'],
                                                        err_formal=fC1['err_formal'], chi2=fC1['chi2'],
                                                        dof=fC1['dof'], n=fC1['n'])
        strength['local_fits'] = loc
        strength['two_strands_fitted'] = bool(two)
        # per tile (calibration test): fixed-effects amplitude per tile in SB and in flux units
        per_tile = {}
        for t_ in np.unique(tile[idx][okc]):
            mt = okc & (tile[idx] == t_)
            if mt.sum() < 20 or np.sum(Tr[mt] >= T_ON) < 5:
                continue
            Xt, ugt = fe_design(Tr[mt][:, None], gl_r[mt])
            ft = wls_clip_boot(y[mt], e[mt], Xt, n_report=1, nboot=500)
            ff = wls_clip_boot(O[mt], EO[mt], Xt, n_report=1, nboot=500)
            per_tile[str(int(t_))] = dict(n=int(mt.sum()), program=str(meta['program'][idx][mt][0]),
                                          f_acc_5008=float(np.median(cm['f_chrom_5008'][idx][mt])),
                                          acceptance_filled=bool(filled[idx][mt][0]),
                                          amp_R=[float(ft['b'][0]), float(ft['err_boot'][0])],
                                          amp_flux_1e17=[float(ff['b'][0]), float(ff['err_boot'][0])])
        strength['per_tile_fixed_effects'] = per_tile
        if '80644' in per_tile and '80683' in per_tile:
            a1, a2 = per_tile['80644'], per_tile['80683']
            rf = a1['amp_flux_1e17'][0] / a2['amp_flux_1e17'][0]
            erf = rf * np.hypot(a1['amp_flux_1e17'][1] / a1['amp_flux_1e17'][0],
                                a2['amp_flux_1e17'][1] / a2['amp_flux_1e17'][0])
            strength['tile_flux_ratio_80644_over_80683'] = dict(
                measured=[float(rf), float(erf)], expected_from_adopted_acceptances=float(a2['f_acc_5008'] /
                                                                                         a1['f_acc_5008']),
                implied_f_80644_over_f_80683=[float(1 / rf), float(erf / rf ** 2)],
                note='for a uniform source F_cal = SB Omega / f, so F(80644)/F(80683) = f(80683)/f(80644)')
        # fibers of tiles WITH a measured GFA acceptance (main-survey tiles; P1 only)
        mt = okc & ~filled[idx]
        if mt.sum() >= 20:
            Xt, ugt = fe_design(Tr[mt][:, None], gl_r[mt])
            ft = wls_clip_boot(y[mt], e[mt], Xt, n_report=1, nboot=500)
            strength['measured_gfa_tiles_only'] = dict(
                tiles=sorted({int(t_) for t_ in tile[idx][mt]}), n=int(mt.sum()),
                f_acc_5008_range=[float(cm['f_chrom_5008'][idx][mt].min()), float(cm['f_chrom_5008'][idx][mt].max())],
                amp_R=[float(ft['b'][0]), float(ft['err_boot'][0])])
        # calibration variants: the two SV1 tiles (acceptance not measured) with one common acceptance f_c
        msv = okc & np.isin(tile[idx], [80644, 80683])
        if msv.sum() >= 40:
            Xt, ugt = fe_design(Tr[msv][:, None], gl_r[msv])
            fF = wls_clip_boot(O[msv], EO[msv], Xt, n_report=1, nboot=500)
            aF, eaF = float(fF['b'][0]), float(fF['err_boot'][0])
            cal = dict(sv1_amp_flux_1e17=[aF, eaF], n=int(msv.sum()),
                       note='SB amplitude (R) = flux amplitude x 0.7428 x f_c for a common acceptance f_c of '
                            'tiles 80644 and 80683')
            for lab, fc in (('f_dark_fill_0.608', 0.608), ('f_bright_fill_0.419', 0.419),
                            ('f_mean_of_fills_0.514', 0.5135)):
                cal[lab] = [aF * K_FLUX2R * fc, eaF * K_FLUX2R * fc]
            if 'measured_gfa_tiles_only' in strength:
                aM, eaM = strength['measured_gfa_tiles_only']['amp_R']
                fx = aM / (aF * K_FLUX2R)
                cal['f_common_implied_by_measured_tiles'] = [fx, fx * np.hypot(eaM / aM, eaF / aF)]
            strength['calibration_variants'] = cal
        res['strength'] = strength

        # profile across the E ridge (binned residuals and model)
        xb = np.arange(-0.6, 0.601, 0.1)
        prof = dict(x=[], n=[], mean=[], err=[], model=[])
        for lo, hi in zip(xb[:-1], xb[1:]):
            m_ = okr & (offE[idx] >= lo) & (offE[idx] < hi)
            if m_.sum() >= 3:
                w_ = 1 / e[m_] ** 2
                prof['x'].append(float(np.mean(offE[idx][m_])))
                prof['n'].append(int(m_.sum()))
                prof['mean'].append(wmean(y[m_], w_))
                prof['err'].append(float(max(w_.sum() ** -0.5, np.std(y[m_]) / np.sqrt(m_.sum()))))
                prof['model'].append(wmean(model_resid[idx][m_], w_))
        res['profile_across_E'] = {k: rnd(v, 4) for k, v in prof.items()}

        # ---------------------------------------------------------------- (e) line ratios
        ratios = {}
        # --- (e1) regression (instrumental-variable) estimator of the paper, restricted to the region
        Zf_r = Zf[idx]

        def iv_region(fe):
            defs = dict(LINES)
            defs.update({g_: (g_, 0, null_ck(g_)) for g_ in NULLS})
            cnts = np.array([np.bincount(RNG.integers(0, nreg, nreg), minlength=nreg) for _ in range(NBOOT)],
                            float)
            best, boot = {}, {}
            for nm, (g_, col, ck) in defs.items():
                a = fits[g_]['amp'][:, col]
                er = fits[g_]['err'][:, col] * sc
                cch = cm[f'f_chrom_{ck}'][idx] / cm['f_chrom_5008'][idx]
                yy, ey = a * cch, er * cch
                m_ = gb & good_fit(lf, g_, col)[idx] & np.isfinite(cch) & np.isfinite(yy) & np.isfinite(ey) & (ey > 0)
                if fe:
                    cg = {g2: int((gl_r[m_] == g2).sum()) for g2 in np.unique(gl_r[m_])}
                    m_ &= np.array([cg.get(g2, 0) >= 5 for g2 in gl_r])
                u = 1 / ey[m_] ** 2
                grp = np.unique(gl_r[m_], return_inverse=True)[1] if fe else None
                num, den = iv(yy[m_], Zf_r[m_], O[m_], u, grp)
                nb_, db_ = np.empty(NBOOT), np.empty(NBOOT)
                for b in range(NBOOT):
                    nb_[b], db_[b] = iv(yy[m_], Zf_r[m_], O[m_], u * cnts[b][m_], grp)
                best[nm] = (np.array([num / den]), np.array([np.nan]))
                boot[nm] = (nb_ / db_, None)
                best[nm + '__n'] = int(m_.sum())
            best['O3'] = (np.array([1.0]), np.array([np.nan]))
            boot['O3'] = (np.ones(NBOOT), None)
            return best, boot
        for lab, fe in (('regression_intercept', False), ('regression_petal_FE', True)):
            best, boot = iv_region(fe)
            tab = ratio_table_from_boot(best, boot, NULLS, lab)
            for nm in LINES:
                tab[nm]['n'] = best[nm + '__n']
                tab[nm].pop('formal_err', None)
            ratios[lab] = add_limits(tab)

        # --- (e2) direct ON - OFF stack per tile-petal.  OFF = same tile-petal, T < 0.05, anywhere on that petal
        all_idx = np.arange(N)
        on_local = np.where(gb & (Ttot[idx] >= T_ON))[0]          # positions within the region arrays
        on_gl = gl_r[on_local]
        pet = [g_ for g_ in np.unique(on_gl) if (on_gl == g_).sum() >= 5]
        off_global = {g_: all_idx[base & (glab == g_) & (Ttot < T_OFF)] for g_ in pet}
        pet = [g_ for g_ in pet if off_global[g_].size >= 5]
        # build a stacker over ON (region) + OFF (same petals) fibers
        on_g_idx = [idx[on_local[on_gl == g_]] for g_ in pet]
        off_g_idx = [off_global[g_] for g_ in pet]
        used = np.unique(np.concatenate(on_g_idx + off_g_idx))
        loc_of = {g_: k for k, g_ in enumerate(used)}
        st2 = Stacker(Rall[used], Iall[used], meta['redshift'][used], fnodes_all[used], lam_nodes)
        onL = [np.array([loc_of[i_] for i_ in a_]) for a_ in on_g_idx]
        offL = [np.array([loc_of[i_] for i_ in a_]) for a_ in off_g_idx]
        Sd, Id = st2.onoff(onL, offL)
        fd = fit_stack(Sd, Id, sigA, V_NEB, FIT_GROUPS, vscan=vscan)
        best_d = stack_line_fluxes(fd)
        nbd = 500
        Sb, Ib = st2.onoff(onL, offL, boot=nbd)
        fb = fit_stack(Sb, Ib, sigA, V_NEB, FIT_GROUPS, vscan=vscan)
        boot_d = stack_line_fluxes(fb)
        tabd = ratio_table_from_boot(best_d, boot_d, NULLS, 'direct_on_minus_off_stack')
        # OFF-OFF split null: scatter of line fluxes with no ridge signal, in units of the real [O III] flux
        Sn, In = st2.onoff(onL, offL, split_null=300)
        fn = fit_stack(Sn, In, sigA, V_NEB, FIT_GROUPS)
        null_s = stack_line_fluxes(fn)
        o3d = best_d['O3'][0][0]
        split = {}
        for nm in ['O3'] + list(LINES):
            a_, e_ = null_s[nm]
            okn = np.isfinite(a_) & np.isfinite(e_) & (e_ > 0)
            # scale the null scatter to the noise of the real ON-OFF configuration
            sc_ = best_d[nm][1][0] / np.median(e_[okn])
            split[nm] = dict(null_mean_over_o3=float(np.mean(a_[okn]) / o3d),
                             null_std_over_o3_scaled=float(np.std(a_[okn]) * sc_ / o3d),
                             chi_std=float(np.std(a_[okn] / e_[okn])))
        tabd['_offoff_split_null'] = split
        ex_sig = {nm: split[nm]['null_std_over_o3_scaled'] for nm in LINES}
        ex_sig['OII3727'] = float(np.hypot(ex_sig['OII3726'], ex_sig['OII3729']))
        ex_sig['SII6716+6731'] = float(np.hypot(ex_sig['SII6716'], ex_sig['SII6731']))
        tabd = add_limits(tabd, extra=ex_sig)
        for nm in LINES:      # significance against the OFF-OFF split null as well
            tabd[nm]['sig_split_null'] = split[nm]['null_std_over_o3_scaled']
        mT_on = np.mean(np.concatenate([Ttot[a_] for a_ in on_g_idx]))
        mT_off = np.mean(np.concatenate([Ttot[a_] for a_ in off_g_idx]))
        vbd = fb['o3']['v_free']
        vbd = vbd[np.isfinite(vbd)]
        tabd['_stack'] = dict(
            petals=pet, n_on=[int(a_.size) for a_ in on_g_idx], n_off=[int(a_.size) for a_ in off_g_idx],
            off_median_sep_deg=float(np.median(np.hypot((ra[np.concatenate(off_g_idx)] - x0) *
                                                        np.cos(np.radians(y0)), dec[np.concatenate(off_g_idx)] - y0))),
            OIII5007_on_minus_off_R=[float(o3d), float(np.std(boot_d['O3'][0]))],
            formal_err_R=float(best_d['O3'][1][0]),
            mean_T_on=float(mT_on), mean_T_off=float(mT_off),
            implied_ridge_peak_R=float(o3d / (mT_on - mT_off)),
            implied_ridge_peak_err_R=float(np.std(boot_d['O3'][0]) / (mT_on - mT_off)),
            v_helio=[float(fd['o3']['v_free'][0]), float(np.std(vbd))],
            v_formal_err=float(fd['o3']['ev_free'][0]))
        ratios['direct_on_minus_off_stack'] = tabd

        # ON - OFF velocities of sub-samples (the petal mean cancels exactly in ON - OFF)
        def onoff_velocity(on_mask_region, tiles=None, nb=300):
            onl = np.where(on_mask_region)[0]
            ogl = gl_r[onl]
            pp = [g_ for g_ in np.unique(ogl) if (ogl == g_).sum() >= 5 and
                  (tiles is None or int(g_.split('_')[0]) in tiles)]
            offd = {g_: all_idx[base & (glab == g_) & (Ttot < T_OFF)] for g_ in pp}
            pp = [g_ for g_ in pp if offd[g_].size >= 5]
            if not pp:
                return None
            ong = [idx[onl[ogl == g_]] for g_ in pp]
            ofg = [offd[g_] for g_ in pp]
            us = np.unique(np.concatenate(ong + ofg))
            lo_ = {g_: k for k, g_ in enumerate(us)}
            stq = Stacker(Rall[us], Iall[us], meta['redshift'][us], fnodes_all[us], lam_nodes)
            a1 = [np.array([lo_[i_] for i_ in a_]) for a_ in ong]
            a0 = [np.array([lo_[i_] for i_ in a_]) for a_ in ofg]
            Sq, Iq = stq.onoff(a1, a0)
            fq = fit_stack(Sq, Iq, sigA, V_NEB, ['o3'], vscan=vscan)['o3']
            Sqb, Iqb = stq.onoff(a1, a0, boot=nb)
            fqb = fit_stack(Sqb, Iqb, sigA, V_NEB, ['o3'], vscan=vscan)['o3']
            vb_ = fqb['v_free'][np.isfinite(fqb['v_free'])]
            return dict(n_on=int(sum(a_.size for a_ in ong)), n_off=int(sum(a_.size for a_ in ofg)), petals=pp,
                        v=float(fq['v_free'][0]), err_boot=float(np.std(vb_)), err_formal=float(fq['ev_free'][0]),
                        flux_R=float(fq['amp_vfree'][0]), err_flux_boot_R=float(np.nanstd(fqb['amp_vfree'])))
        on_m = gb & (Ttot[idx] >= T_ON)
        sub = {'all_on_ridge': onoff_velocity(on_m),
               'tile_80644': onoff_velocity(on_m, tiles=[80644]),
               'tile_80683': onoff_velocity(on_m, tiles=[80683]),
               'tiles_with_measured_acceptance': onoff_velocity(on_m & ~filled[idx]),
               'table2_fibers_sn_gt4': onoff_velocity(sel),
               'on_ridge_sn_le4': onoff_velocity(on_m & ~(snr > SN_MIN)),
               'east_of_E_spine': onoff_velocity(on_m & (offE[idx] >= 0)),
               'west_of_E_spine': onoff_velocity(on_m & (offE[idx] < 0)),
               'south_half': onoff_velocity(on_m & (dec[idx] < y0)),
               'north_half': onoff_velocity(on_m & (dec[idx] >= y0)),
               'E_only_TW_lt_0.1': onoff_velocity(on_m & (TW[idx] < 0.1)),
               'W_only_TE_lt_0.1': onoff_velocity(on_m & (TE[idx] < 0.1))}
        res['velocity']['on_minus_off_stacks'] = {k: v_ for k, v_ in sub.items() if v_ is not None}
        tsel = Ttot[idx][sel]
        psel = (cm['T'] - cm['Z'])[idx][sel]
        res['velocity']['petal_mean_amplification'] = dict(
            median_T=float(np.median(tsel)), median_petal_mean_T=float(np.median(psel)),
            median_amplification=float(np.median(tsel / (tsel - psel))),
            note='a residual a_i G(v_i) - abar G(vbar) has centroid (a_i v_i - abar vbar)/(a_i - abar): '
                 'velocity differences from the petal mean are amplified by T/(T - <T>_petal)')
        # wavelength zero point: tile-to-tile ON - OFF difference (independent nights) and the sky-line test
        zp = {}
        if sub.get('tile_80644') and sub.get('tile_80683'):
            a_, b_ = sub['tile_80644'], sub['tile_80683']
            zp['on_off_tile80644_minus_tile80683'] = [a_['v'] - b_['v'], float(np.hypot(a_['err_boot'],
                                                                                        b_['err_boot']))]
        if sub.get('tiles_with_measured_acceptance') and sub.get('tile_80683'):
            a_ = sub['tiles_with_measured_acceptance']
            w1 = 1 / sub['tile_80644']['err_boot'] ** 2
            w2 = 1 / sub['tile_80683']['err_boot'] ** 2
            vsv = (w1 * sub['tile_80644']['v'] + w2 * sub['tile_80683']['v']) / (w1 + w2)
            zp['on_off_main_survey_tiles_minus_sv1'] = [a_['v'] - vsv, float(np.hypot(a_['err_boot'],
                                                                                    (w1 + w2) ** -0.5))]
        sky = {}
        l5 = L['sky5577']
        s5 = (WAVE > l5 - 12) & (WAVE < l5 + 12)
        lam5 = WAVE[s5]
        side5 = (np.abs(lam5 - l5) > 5) & (np.abs(lam5 - l5) < 11)
        for t_ in np.unique(tile[idx]):
            mt_ = (tile == t_) & (cm['n_tiles'] == 1)
            with np.errstate(divide='ignore', invalid='ignore'):
                Vv = np.where(Iall[mt_][:, s5] > 0, 1 / Iall[mt_][:, s5], np.nan)
                prof5 = np.nanmedian(Vv / np.nanmedian(Vv[:, side5], 1)[:, None] - 1, 0)
            try:
                from scipy.optimize import curve_fit as _cf
                pp_, pc_ = _cf(lambda x_, a, c0, sg, b0: a * np.exp(-0.5 * ((x_ - c0) / sg) ** 2) + b0, lam5, prof5,
                               p0=[np.nanmax(prof5), l5, 0.9, 0.0])
                js = ex.by_tile[int(t_)]
                sky[str(int(t_))] = dict(n_coadds=int(mt_.sum()), n_exp=len(js),
                                         centroid_kms=float((pp_[1] / l5 - 1) * C),
                                         err_kms=float(np.sqrt(pc_[1, 1]) / l5 * C),
                                         vcorr_efftime_mean=float(np.sum(ex.eff[js] * ex.vb[js]) / np.sum(ex.eff[js])),
                                         vcorr_range=[float(ex.vb[js].min()), float(ex.vb[js].max())])
            except Exception:
                pass
        zp['OI5577_variance_peak_by_tile'] = sky
        zp['OI5577_note'] = ('centroid of the [O I]5577 excess in 1/ivar (median profile of all coadds of the tile); '
                             'expected at the barycentric correction; a saturated, non-Gaussian profile on 0.8 A '
                             'pixels: scatter ~10 km/s, not a zero-point constraint; multi-night coadds smear it')
        res['zero_point'] = zp
        # ON-only stack (relative to the petal sky mean), same fibers, for comparison
        c_on = np.zeros((1, used.size))
        c_on[0, np.concatenate(onL)] = 1
        S1, I1 = st2.stack(c_on)
        f1 = fit_stack(S1, I1, sigA, V_NEB, FIT_GROUPS, vscan=vscan)
        b1 = stack_line_fluxes(f1)
        c_off = np.zeros((1, used.size))
        c_off[0, np.concatenate(offL)] = 1
        S0, I0 = st2.stack(c_off)
        f0 = fit_stack(S0, I0, sigA, V_NEB, FIT_GROUPS, vscan=vscan)
        b0 = stack_line_fluxes(f0)
        ratios['stack_fluxes_R'] = {
            'on': {nm: [float(b1[nm][0][0]), float(b1[nm][1][0])] for nm in ['O3'] + list(LINES)},
            'off': {nm: [float(b0[nm][0][0]), float(b0[nm][1][0])] for nm in ['O3'] + list(LINES)},
            'on_minus_off': {nm: [float(best_d[nm][0][0]), float(best_d[nm][1][0]),
                                  float(np.std(boot_d[nm][0]))] for nm in ['O3'] + list(LINES)},
            'units': 'energy surface brightness in R-equivalent at 5008 A (true R for [O III] 5007); '
                     '[value, formal error(, bootstrap std)]',
            'v_on': float(f1['o3']['v_free'][0]), 'v_off': float(f0['o3']['v_free'][0]),
            'ev_on': float(f1['o3']['ev_free'][0]), 'ev_off': float(f0['o3']['ev_free'][0]),
        }
        res['ratios'] = ratios
        # free-width fit of [O III] 5007 in the ON-OFF stack (diagnostic of the line width)
        try:
            from scipy.optimize import curve_fit
            mm = (WAVE > LAM5007 - 14) & (WAVE < LAM5007 + 14) & (Id[0] > 0)
            s_inst = float(np.hypot(np.interp(LAM5007, WAVE, sigA), PIX / np.sqrt(12)))

            def gmod(x, a, mu, s, c0, c1):
                return a * np.exp(-0.5 * ((x - mu) / s) ** 2) / (np.sqrt(2 * np.pi) * s) + c0 + c1 * (x - LAM5007)
            pw, cw = curve_fit(gmod, WAVE[mm], Sd[0][mm], p0=[o3d, LAM5007 * (1 + 12 / C), 0.76, 0, 0],
                               sigma=1 / np.sqrt(Id[0][mm]))
            sb_ = []
            for b in range(0, nbd, 2):
                mb = mm & (Ib[b] > 0)
                try:
                    pb_, _ = curve_fit(gmod, WAVE[mb], Sb[b][mb], p0=pw, sigma=1 / np.sqrt(Ib[b][mb]))
                    sb_.append(abs(pb_[2]))
                except Exception:
                    pass
            sb_ = np.array(sb_)
            s_obs = abs(float(pw[2]))

            def fwhm_int(s):
                return 2.3548 * np.sqrt(max(s ** 2 - s_inst ** 2, 0.0)) / LAM5007 * C
            res['width'] = dict(sigma_obs_A=s_obs, sigma_obs_boot_A=float(np.std(sb_)), sigma_inst_A=s_inst,
                                fwhm_intrinsic_kms=fwhm_int(s_obs),
                                fwhm_intrinsic_p16_p84=[fwhm_int(np.percentile(sb_, 16)),
                                                        fwhm_int(np.percentile(sb_, 84))],
                                fwhm_intrinsic_upper_97p5=fwhm_int(np.percentile(sb_, 97.5)),
                                v_free_sigma_fit=float((pw[1] / LAM5007 - 1) * C),
                                note='instrumental sigma = LSF (+) pixel; LSF known to 2-3 %, so only an upper limit')
        except Exception as exc:        # pragma: no cover
            res['width'] = dict(error=str(exc))

        # compact headline numbers (conservative over the three ratio estimators)
        est3 = ('regression_intercept', 'regression_petal_FE', 'direct_on_minus_off_stack')
        vt = res['velocity']['table2']
        vo_ = res['velocity']['on_minus_off_stacks']['all_on_ridge']
        alt = [res['velocity']['stack_table2_fibers']['v'], vo_['v'], vt['median'],
               res['velocity']['variants']['pos_sn_gt5']['wmean'], res['velocity']['variants']['pos_sn_gt3']['wmean'],
               res['velocity']['on_minus_off_stacks']['table2_fibers_sn_gt4']['v']]
        res['headline'] = {
            'n_coadds': int(nreg), 'n_good': int(gb.sum()),
            'V_helio_table2': [vt['wmean'], vt['err'], vt['n']],
            'V_helio_on_minus_off_stack': [vo_['v'], vo_['err_boot']],
            'V_helio_alternative_estimators_min_max': [float(min(alt)), float(max(alt))],
            'W26_table4_minus_table2_kms': [W26_OIII[reg][0] - vt['wmean'], float(np.hypot(W26_OIII[reg][1],
                                                                                           vt['err']))],
            'global_model_sky_SB_at_region_R': strength['global_model']['sky_SB_at_region_R'],
            'E_ridge_peak_R_global_B_C': [strength['global_model']['E_ridge_peak_at_dec_R'][0],
                                          float(loc['B_global_petal_term']['amp_R'][0]),
                                          float(loc['C_petal_fixed_effects']['amp_R'][0])],
            'single_amp_fixed_effects_R': [float(fC1['b'][0]), float(fC1['err_boot'][0])],
            'on_minus_off_R': tabd['_stack']['OIII5007_on_minus_off_R'],
            'OIII5007_over_4959': {k_: [ratios[k_]['5007/4959']['ratio'], ratios[k_]['5007/4959']['sig_boot']]
                                   for k_ in est3},
            'OII3727_over_OIII': {k_: [ratios[k_]['OII3727']['ratio'], ratios[k_]['OII3727']['sig_adopt']]
                                  for k_ in est3},
            'OII3729_over_3726': {k_: [ratios[k_]['OII3729/3726']['ratio'], ratios[k_]['OII3729/3726']['sig_boot']]
                                  for k_ in est3},
            'OIII_over_Hb_lower_2sig_min_max': [float(min(ratios[k_]['OIII5007/Hb']['lower_2sig'] for k_ in est3)),
                                                float(max(ratios[k_]['OIII5007/Hb']['lower_2sig'] for k_ in est3))],
            'upper_2sig_over_OIII_max': {nm: float(max(ratios[k_][nm]['upper_2sig'] for k_ in est3))
                                         for nm in ('Hb', 'HeII4686', 'Ha', 'NII6583', 'SII6716+6731')},
            'max_abs_snr': {nm: float(max(abs(ratios[k_][nm]['snr']) for k_ in est3))
                            for nm in ('Hb', 'HeII4686', 'Ha', 'NII6583', 'SII6716+6731')},
            'snr_matched_nulls_by_estimator': {nm: [float(ratios[k_][nm]['snr_matched']) for k_ in est3]
                                               for nm in ('Hb', 'HeII4686', 'Ha', 'NII6583', 'SII6716+6731')},
            'sig_matched_nulls_by_estimator': {nm: [float(ratios[k_][nm]['sig_matched']) for k_ in est3]
                                               for nm in ('Hb', 'HeII4686', 'Ha', 'NII6583', 'SII6716+6731')},
            'OIII_over_Ha_photon_lower_2sig': float(LAM5007 / L['Ha'] /
                                                    max(ratios[k_]['Ha']['upper_2sig'] for k_ in est3)),
        }
        summary['regions'][reg] = res

        # ---------------------------------------------------------------- CSV rows
        det_cache = {}
        for k, i in enumerate(idx):
            row = dict(
                region=reg, sparcl_id=str(ids_lf[i]), targetid=int(meta['targetid'][i]), ra=round(float(ra[i]), 6),
                dec=round(float(dec[i]), 6), survey=str(meta['survey'][i]), program=str(meta['program'][i]),
                spectype=str(meta['spectype'][i]), redshift=round(float(meta['redshift'][i]), 5),
                tile=int(tile[i]), petal=int(cm['group_petal'][i]),
                coadd_exptime_s=round(float(meta['coadd_exptime'][i]), 1), efftime_s=round(float(cm['efftime'][i]), 1),
                vcorr_kms=round(float(cm['vcorr'][i]), 3), f_acc_5008=round(float(cm['f_chrom_5008'][i]), 4),
                acceptance_filled=int(filled[i]), kR_5008=round(float(kR[i]), 4),
                sep_deg=round(float(sep[i]), 4), offset_E_deg=round(float(offE[i]), 4),
                offset_W_deg=round(float(offW[i]), 4), T_E=round(float(TE[i]), 4), T_W=round(float(TW[i]), 4),
                T=round(float(Ttot[i]), 4), Z=round(float(cm['Z'][i]), 4),
                petal_mean_T=round(float(cm['T'][i] - cm['Z'][i]), 4),
                model_resid_R=round(float(model_resid[i]), 4), model_sky_R=round(float(model_sky[i]), 4),
                good_fit=int(gb[k]), on_ridge=int(Ttot[i] >= T_ON), off_ridge=int(Ttot[i] < T_OFF),
                in_table2=int(sel[k]),
                F5007_1e17=round(float(O[k]), 4), eF5007_1e17=round(float(EO[k]), 4), sn=round(float(snr[k]), 3),
                SB_R=round(float(SB[k]), 4), eSB_R=round(float(eSB[k]), 4),
                v_helio_kms=round(float(vfree[k]), 3) if np.isfinite(vfree[k]) else '',
                ev_kms=round(float(evfree[k]), 3) if np.isfinite(evfree[k]) else '',
                lam_centroid_A=round(float(LAM5007 * (1 + vfree[k] / C)), 4) if np.isfinite(vfree[k]) else '',
                dchi2_free_v=round(float(dchi[k]), 3) if np.isfinite(dchi[k]) else '',
                F5007_free_v_1e17=round(float(fits['o3']['amp_vfree'][k]), 4),
                sigma_line_A=round(sig_line_A, 4), rchi2=round(float(fits['o3']['rchi2'][k]), 3),
                npix=int(fits['o3']['npix'][k]))
            for nm, (g_, col, ck) in LINES.items():
                row[f'F_{nm}_1e17'] = round(float(fits[g_]['amp'][k, col]), 4)
                row[f'eF_{nm}_1e17'] = round(float(fits[g_]['err'][k, col] * sc), 4)
            fiber_rows.append(row)

        # ---------------------------------------------------------------- data for the figure
        # example spectra: brightest, median and near-threshold fibers of the Table 2 selection
        si = np.where(sel)[0]
        order = si[np.argsort(-snr[si])]
        picks = [order[0], order[len(order) // 2], order[-1]]
        exs = []
        for k in picks:
            d_ = fit_detail(Rr[k], Ir[k], zr[k], vc[k], G['o3'], V_NEB, sigA, v_eval=vfree[k])
            assert abs(d_['p_fixed'][3] - O[k]) < 1e-6 * max(1.0, abs(O[k])), 'fit_detail differs from fit_group'
            exs.append(dict(lam=d_['lam'], r=d_['r'] * kRr[k], w=d_['w'], w_raw=d_['w_raw'],
                            model=d_['model_free'] * kRr[k], cont=d_['cont_free'] * kRr[k], sn=float(snr[k]),
                            v=float(vfree[k]), ev=float(evfree[k]), SB=float(SB[k]), offE=float(offE[idx][k]),
                            program=str(meta['program'][idx][k]), spectype=str(meta['spectype'][idx][k]),
                            z=float(zr[k]), id=str(ids_lf[idx][k])))
        figdata[reg] = dict(ids=ids_lf[idx], ra=ra[idx], dec=dec[idx], SB=SB, eSB=eSB, good=gb, sel=sel, snr=snr, v=vfree, ev=evfree,
                            offE=offE[idx], T=Tr, examples=exs, Sd=Sd[0], Id=Id[0], S1=S1[0], I1=I1[0], S0=S0[0],
                            I0=I0[0], fd=fd, model_resid=model_resid[idx], okr=okr)
        stk_out.update({f'{reg}_on_minus_off': Sd[0], f'{reg}_on_minus_off_ivar': Id[0], f'{reg}_on': S1[0],
                        f'{reg}_on_ivar': I1[0], f'{reg}_off': S0[0], f'{reg}_off_ivar': I0[0],
                        f'{reg}_table2_fibers': S_sel, f'{reg}_table2_fibers_ivar': I_sel})
        print(reg, 'done  %.0fs' % (time.time() - t0), flush=True)

    # ---------------------------------------------------------------- write products
    with open(OUT_CSV, 'w', newline='') as fh:
        wri = csv.DictWriter(fh, fieldnames=list(fiber_rows[0].keys()))
        wri.writeheader()
        wri.writerows(fiber_rows)
    OUT_JSON.write_text(json.dumps(summary, indent=1, default=jd))
    np.savez_compressed(OUT_NPZ, **stk_out)
    print('wrote', OUT_CSV.name, len(fiber_rows), 'rows;', OUT_JSON.name, ';', OUT_NPZ.name)
    brief(summary)
    make_figure(figdata, summary, T, sigA, V_NEB)


def brief(s):
    for reg in REG_NAMES:
        r = s['regions'][reg]
        v = r['velocity']['table2']
        print(f"\n=== {reg} ===  n={r['inventory']['fetched_and_analyzed']} good={r['inventory']['usable_good_o3_fit']}")
        print('  Table 2: %+.2f +- %.2f (n=%d; formal %.2f, scatter %.2f, boot %.2f) median %+.2f MAD %.1f '
              'med.err %.1f chi2 %.2f/%.2f extra %.1f' % (v['wmean'], v['err'], v['n'], v['err_formal'],
                                                         v['err_scatter'], v['err_bootstrap'], v['median'],
                                                         v['mad_sigma'], v['median_formal_err'], v['redchi2_formal'],
                                                         v['redchi2_with_floor'], v['extra_scatter_kms']))
        for k in ('stack_table2_fibers', 'stack_on_ridge_T_ge_0.5'):
            q = r['velocity'][k]
            print('  %s: %+.2f +- %.2f (formal %.2f) n=%d flux %.3f R' % (k, q['v'], q['err_boot'], q['err_formal'],
                                                                        q['n'], q['flux_R']))
        for k, q in r['velocity']['variants'].items():
            print('     %-28s n=%3d %+6.2f +- %.2f  median %+6.2f' % (k, q['n'], q['wmean'], q['err'], q['median']))
        st = r['strength']
        print('  strength:', json.dumps({k: rnd(v, 3) for k, v in st.items() if k.startswith('resid')}, default=jd))
        print('   global:', json.dumps({k: rnd(v, 3) for k, v in st['global_model'].items() if k != 'note'},
                                       default=jd))
        for k, q in st['local_fits'].items():
            print('   %-36s amp %s +- %s (formal %s) chi2/dof %.0f/%d' % (k, rnd(q['amp_R'], 3), rnd(q['err_boot'], 3),
                                                                        rnd(q['err_formal'], 3), q['chi2'], q['dof']))
        print('   per tile:', json.dumps({k: {kk: rnd(vv, 3) if not isinstance(vv, (str, bool)) else vv
                                              for kk, vv in q.items()} for k, q in
                                          st['per_tile_fixed_effects'].items()}, default=jd))
        for k in ('tile_flux_ratio_80644_over_80683', 'measured_gfa_tiles_only', 'calibration_variants'):
            if k in st:
                print('   ', k, json.dumps({kk: rnd(vv, 3) if not isinstance(vv, str) else vv
                                            for kk, vv in st[k].items()}, default=jd))
        for k, q in r['velocity']['on_minus_off_stacks'].items():
            print('     ON-OFF %-32s n_on %3d n_off %3d  v %+6.2f +- %.2f  flux %.2f R' % (
                k, q['n_on'], q['n_off'], q['v'], q['err_boot'], q['flux_R']))
        print('   amplification', json.dumps(rnd_dict(r['velocity']['petal_mean_amplification'])))
        print('   zero point', json.dumps(r['zero_point'], default=lambda o: rnd(o, 2))[:900])
        for lab in ('regression_intercept', 'regression_petal_FE', 'direct_on_minus_off_stack'):
            t = r['ratios'][lab]
            print('  ratios [%s]  infl %.2f  null rms %.3f max %.3f chi %.2f' % (lab, t['_inflation'],
                                                                              t['_nulls']['rms_ratio'],
                                                                              t['_nulls']['max_abs'],
                                                                              t['_nulls']['chi_rms']))
            for nm in ('OIII4959', 'Hb', 'HeII4686', 'OII3727', 'OII3726', 'OII3729', 'Ha', 'NII6583',
                       'SII6716+6731'):
                d = t[nm]
                print('     %-13s %+.4f +- %.4f (boot) %.4f (adopt)  snr %+.1f  up2s %.3f' % (
                    nm, d['ratio'], d['sig_boot'], d['sig_adopt'], d['snr'], d['upper_2sig']))
            print('     5007/4959 %.3f (%.3f-%.3f)   3729/3726 %.3f (%.3f-%.3f)   [O III]/Hb > %.1f' % (
                t['5007/4959']['ratio'], t['5007/4959']['p16'], t['5007/4959']['p84'], t['OII3729/3726']['ratio'],
                t['OII3729/3726']['p16'], t['OII3729/3726']['p84'], t['OIII5007/Hb']['lower_2sig']))
        print('  direct stack:', json.dumps({k: rnd(v, 3) if not isinstance(v, list) or not isinstance(v[0], str)
                                             else v for k, v in r['ratios']['direct_on_minus_off_stack']['_stack'].
                                             items()}, default=jd))
        print('  width:', json.dumps(rnd_dict(r['width']), default=jd))


def rnd_dict(d):
    return {k: (rnd(v, 3) if not isinstance(v, str) else v) for k, v in d.items()}


def make_figure(figdata, summary, T, sigA, V_NEB):
    """Appendix figure: rows 1-2 = P1, P3 (map; individual residual spectra; per-fiber velocities; ON-OFF stack);
    row 3 = cross-ridge brightness profile, [O II] and Halpha-[N II] stacks, line ratios.  7.1 in wide."""
    import silk_style as ss
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    from matplotlib.lines import Line2D
    FS = dict(tick=6.5, label=7.0, tag=7.2, note=6.5)
    CREG = {'P1': ss.OIII_C, 'P3': ss.SII_C}
    G = groups(V_NEB)
    lw_ax = 0.6
    fig = plt.figure(figsize=(7.1, 7.35))
    gs = fig.add_gridspec(3, 4, left=0.068, right=0.992, top=0.972, bottom=0.058, hspace=0.36, wspace=0.42,
                          height_ratios=[1, 1, 0.92])
    letters = iter('abcdefghijkl')

    def style(ax):
        ax.tick_params(labelsize=FS['tick'], length=2.5, width=lw_ax, pad=2)
        for sp in ax.spines.values():
            sp.set_linewidth(lw_ax)

    def tag(ax, text, x=0.0, y=1.035):
        ax.text(x, y, text, transform=ax.transAxes, ha='left', va='bottom', fontsize=FS['tag'], fontweight='bold',
                color=ss.TITLE)

    def note(ax, x, y, text, **kw):
        kw.setdefault('fontsize', FS['note'])
        kw.setdefault('color', ss.FG)
        return ax.text(x, y, text, transform=ax.transAxes, path_effects=[ss.HALO], **kw)

    def cont_sub(S, Iv, gname, v=None):
        d = fit_detail(S, Iv, -1.0, 0.0, G[gname], V_NEB, sigA, v_eval=v)
        if v is not None:
            return d['lam'], d['r'] - d['cont_free'], d['model_free'] - d['cont_free']
        return d['lam'], d['r'] - d['cont_fixed'], d['model_fixed'] - d['cont_fixed']

    strands = {s: vec2radec(P) for s, tree, sig, a, P in T.parts}
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-1.5, vmax=3.0)
    cmap = plt.get_cmap('RdBu_r')
    for row, reg in enumerate(REG_NAMES):
        fd = figdata[reg]
        R = summary['regions'][reg]
        x0, y0 = REGIONS[reg]
        cd = np.cos(np.radians(y0))
        # ---------------- map
        ax = fig.add_subplot(gs[row, 0])
        style(ax)
        hw_map = 0.56
        g = fd['good']
        order = np.argsort(np.abs(fd['SB'][g]))
        sc_ = ax.scatter(fd['ra'][g][order], fd['dec'][g][order], c=fd['SB'][g][order], s=5.5, cmap=cmap, norm=norm,
                         linewidths=0.25, edgecolors=ss.FGDIM, zorder=3)
        for sname, (sra, sdec) in strands.items():
            if sname == 'W2':
                continue
            ks = sdec > y0 - hw_map - 0.02
            ax.plot(sra[ks], sdec[ks], ':', color=ss.FG, lw=0.9, zorder=2)
        th = np.linspace(0, 2 * np.pi, 181)
        ax.plot(x0 + RADIUS * np.cos(th) / cd, y0 + RADIUS * np.sin(th), '--', color=ss.FG, lw=0.7, zorder=4)
        ax.plot([x0], [y0], 's', mfc='none', mec=ss.FG, mew=0.9, ms=5, zorder=5)
        for e_ in fd['examples']:
            k = np.where(np.array([str(i_) for i_ in fd['ids']]) == e_['id'])[0]
            ax.plot(fd['ra'][k], fd['dec'][k], 'o', mfc='none', mec=ss.FG, mew=0.8, ms=6.5, zorder=6)
        hw = hw_map
        bb = ax.get_position()
        hwy = hw * (bb.height * fig.get_figheight()) / (bb.width * fig.get_figwidth())
        ax.set_xlim(x0 + hw / cd, x0 - hw / cd)
        ax.set_ylim(y0 - hwy, y0 + hwy)          # true aspect; margins hold the color bar and strand labels
        ax.set_xlabel('RA (deg)', fontsize=FS['label'], labelpad=1.5)
        ax.set_ylabel('Dec (deg)', fontsize=FS['label'], labelpad=2)
        ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
        ax.yaxis.set_major_locator(plt.MultipleLocator(0.5))
        tag(ax, '(%s) %s fibers' % (next(letters), reg))
        cax = ax.inset_axes([0.14, 0.115, 0.66, 0.035])
        cb = fig.colorbar(sc_, cax=cax, orientation='horizontal', ticks=[-1, 0, 1, 2, 3])
        cb.ax.tick_params(labelsize=FS['tick'], length=1.5, width=0.5, pad=1)
        cb.outline.set_linewidth(0.5)
        cax.xaxis.set_ticks_position('bottom')
        cax.text(1.06, 0.5, 'R', transform=cax.transAxes, fontsize=FS['note'], va='center', ha='left', color=ss.FG,
                 path_effects=[ss.HALO])
        for lab, (sra, sdec) in strands.items():
            if lab == 'W2':
                continue
            yy = ax.get_ylim()[1] - 0.07
            k = np.argmin(np.abs(sdec - yy))
            if abs(sdec[k] - yy) < 0.05 and abs(sra[k] - x0) * cd < hw - 0.05:
                ax.text(sra[k] - 0.035 / cd, yy, lab, fontsize=FS['note'], ha='left', va='center', color=ss.FG,
                        path_effects=[ss.HALO], zorder=7)

        # ---------------- individual spectra
        ax = fig.add_subplot(gs[row, 1])
        style(ax)
        base_y = 0.0
        ymax = 0
        for e_ in fd['examples'][::-1]:
            m = (e_['lam'] > 4949.5) & (e_['lam'] < 5020.5)
            rr = (e_['r'] - e_['cont'])[m]
            mm = (e_['model'] - e_['cont'])[m]
            used = e_['w'][m] > 0
            rr_plot = np.where(e_['w_raw'][m] > 0, rr, np.nan)
            ax.step(e_['lam'][m], rr_plot + base_y, where='mid', color=ss.FGDIM, lw=0.7)
            ax.plot(e_['lam'][m], mm + base_y, color=ss.OIII_C, lw=1.0)
            ax.axhline(base_y, color=ss.GRID, lw=0.5, zorder=0)
            top = max(np.nanmax(mm), np.nanpercentile(rr_plot[np.isfinite(rr_plot)], 99))
            ax.text(4984.0, base_y + 0.27, 'S/N %.1f\n%+.2f°\n%+.0f km/s' % (e_['sn'], e_['offE'], e_['v']),
                    fontsize=FS['note'], va='bottom', ha='center', color=ss.FG, path_effects=[ss.HALO],
                    linespacing=1.05)
            ymax = base_y + max(top, 0.95)
            base_y += max(top, 0.95) + 0.30
        for l0 in (L['OIII4959'], L['OIII5007']):
            ax.axvline(l0, color=ss.FGDIM, lw=0.5, ls=':', zorder=0)
        ax.set_xlim(4949.5, 5020.5)
        ax.set_ylim(-0.30, ymax * 1.10 + 0.05)
        ax.xaxis.set_major_locator(plt.MultipleLocator(20))
        ax.xaxis.set_minor_locator(plt.MultipleLocator(5))
        ax.set_xlabel('vacuum λ (Å)', fontsize=FS['label'], labelpad=1.5)
        ax.set_ylabel('residual (R Å⁻¹) + offset', fontsize=FS['label'], labelpad=2)
        tag(ax, '(%s) single fibers' % next(letters))

        # ---------------- per-fiber velocities
        ax = fig.add_subplot(gs[row, 2])
        style(ax)
        ok = fd['good'] & np.isfinite(fd['v']) & np.isfinite(fd['ev']) & (fd['ev'] < EV_MAX) & (fd['snr'] > 2.5)
        lo = ok & ~fd['sel']
        ax.errorbar(fd['snr'][lo], fd['v'][lo], yerr=fd['ev'][lo], fmt='o', ms=1.8, color=ss.FGDIM, alpha=0.45,
                    lw=0.4, mew=0, zorder=2)
        ax.errorbar(fd['snr'][fd['sel']], fd['v'][fd['sel']], yerr=fd['ev'][fd['sel']], fmt='o', ms=2.2,
                    color=ss.OIII_C, lw=0.45, mew=0, zorder=3)
        v2 = R['velocity']['table2']
        ax.axhspan(V_ADOPT - EV_ADOPT, V_ADOPT + EV_ADOPT, color=ss.FGDIM, alpha=0.30, lw=0, zorder=1)
        ax.axhline(v2['wmean'], color=ss.FG, lw=1.1, zorder=4)
        ax.axhline(W26_OIII[reg][0], color=ss.DANGER, lw=0.9, ls='--', zorder=4)
        ax.axvline(SN_MIN, color=ss.FGDIM, lw=0.6, ls=':', zorder=1)
        ax.set_xscale('log')
        ax.set_xlim(2.5, 26)
        ax.set_ylim(-52, 78)
        ax.set_xticks([3, 4, 6, 10, 20])
        ax.set_xticklabels(['3', '4', '6', '10', '20'])
        ax.xaxis.set_minor_locator(plt.NullLocator())
        ax.yaxis.set_major_locator(plt.MultipleLocator(25))
        ax.set_xlabel('[O III] S/N', fontsize=FS['label'], labelpad=1.5)
        ax.set_ylabel('V$_{helio}$ (km s⁻¹)', fontsize=FS['label'], labelpad=0.5)
        note(ax, 0.975, 0.965, 'W26 %+.1f' % W26_OIII[reg][0], ha='right', va='top', color=ss.DANGER)
        vo = R['velocity']['on_minus_off_stacks']['all_on_ridge']
        note(ax, 0.975, 0.035, 'mean %+.1f ± %.1f (%d)\nON−OFF %+.1f ± %.1f' % (v2['wmean'], v2['err'], v2['n'],
                                                                           vo['v'], vo['err_boot']),
             ha='right', va='bottom', linespacing=1.1)
        tag(ax, '(%s) fiber velocities' % next(letters))

        # ---------------- ON-OFF stack around Hbeta, 4959, 5007 (three sub-windows)
        sub = gs[row, 3].subgridspec(1, 3, wspace=0.10)
        vstk = R['ratios']['direct_on_minus_off_stack']['_stack']['v_helio'][0]
        segs = [('hb', L['Hb'], 'Hβ'), ('o3', L['OIII4959'], '4959'), ('o3', L['OIII5007'], '5007')]
        axes3 = []
        for k, (gname, l0, lab) in enumerate(segs):
            a3 = fig.add_subplot(sub[0, k])
            style(a3)
            for S_, I_, col, lw_, z_ in ((fd['S1'], fd['I1'], ss.OIII_C, 0.7, 2), (fd['S0'], fd['I0'], ss.SII_C, 0.7, 2),
                                         (fd['Sd'], fd['Id'], ss.FG, 0.9, 3)):
                lam, r_, m_ = cont_sub(S_, I_, gname, v=vstk if gname == 'o3' else None)
                q = np.abs(lam - l0) < 8.5
                a3.step(lam[q], r_[q], where='mid', color=col, lw=lw_, zorder=z_)
                if col == ss.FG:
                    a3.plot(lam[q], m_[q], color=ss.DANGER, lw=0.9, zorder=4)
            a3.axhline(0, color=ss.GRID, lw=0.5, zorder=0)
            if k:
                a3.patch.set_alpha(0.0)
            a3.set_xlim(l0 - 7.5, l0 + 7.5)
            a3.set_xticks([np.round(l0)])
            a3.set_xticklabels(['%d' % np.round(l0)])
            a3.xaxis.set_minor_locator(plt.MultipleLocator(5))
            if k:
                a3.tick_params(labelleft=False, left=False)
                a3.spines['left'].set_visible(False)
            if k < 2:
                a3.spines['right'].set_visible(False)
            a3.text(0.5, 0.97, lab, transform=a3.transAxes, ha='center', va='top', fontsize=FS['note'], color=ss.FG,
                    path_effects=[ss.HALO])
            axes3.append(a3)
        ylo = min(a.get_ylim()[0] for a in axes3)
        yhi = max(a.get_ylim()[1] for a in axes3)
        for a in axes3:
            a.set_ylim(ylo, yhi * 1.12)
        axes3[0].set_ylabel('stack (R Å⁻¹)', fontsize=FS['label'], labelpad=1)
        axes3[1].set_xlabel('vacuum λ (Å)', fontsize=FS['label'], labelpad=1.5)
        tag(axes3[0], '(%s) %s stacks' % (next(letters), reg))
        if row == 0:
            hl = [Line2D([], [], color=ss.FG, lw=0.9), Line2D([], [], color=ss.OIII_C, lw=0.7),
                  Line2D([], [], color=ss.SII_C, lw=0.7), Line2D([], [], color=ss.DANGER, lw=0.9)]
            axes3[0].legend(hl, ['ON−OFF', 'ON', 'OFF', 'fit'], loc='center left', bbox_to_anchor=(-0.02, 0.62),
                            fontsize=FS['note'], frameon=False, handlelength=1.0, handletextpad=0.4,
                            labelspacing=0.2, borderpad=0)

    # ---------------- row 3
    # (i) brightness against offset from the E spine
    ax = fig.add_subplot(gs[2, 0])
    style(ax)
    for reg, mk in zip(REG_NAMES, ('o', 's')):
        pr = summary['regions'][reg]['profile_across_E']
        ax.plot(pr['x'], pr['model'], '-', color=CREG[reg], lw=0.9, alpha=0.9, zorder=2)
        ax.errorbar(pr['x'], pr['mean'], yerr=pr['err'], fmt=mk, ms=3.0, color=CREG[reg], lw=0.7, mew=0, zorder=3,
                    label=reg)
    ax.axhline(0, color=ss.FGDIM, lw=0.5, zorder=1)
    ax.set_xlim(-0.62, 0.62)
    ax.set_xlabel('offset from E spine (deg)', fontsize=FS['label'], labelpad=1.5)
    ax.set_ylabel('[O III] residual (R)', fontsize=FS['label'], labelpad=2)
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
    ax.legend(loc='upper right', fontsize=FS['note'], frameon=False, handlelength=1.0, handletextpad=0.3,
              borderpad=0.1, labelspacing=0.2)
    note(ax, 0.975, 0.74, 'lines:\nmodel', ha='right', va='top', color=ss.FGDIM, linespacing=1.05)
    tag(ax, '(%s) cross-ridge profile' % next(letters))

    # (j), (k) ON-OFF stacks at [O II] and Halpha-[N II]
    for col, (gname, lo_, hi_, marks, title, loc_) in enumerate((
            ('o2', 3717.5, 3739.5, [(L['OII3726'], '3726 '), (L['OII3729'], ' 3729')], '[O II] stacks', 5),
            ('han2', 6541.0, 6594.0, [(L['NII6548'], ''), (L['Ha'], 'Hα'), (L['NII6583'], '[N II]')],
             'Hα, [N II] stacks', 20))):
        ax = fig.add_subplot(gs[2, 1 + col])
        style(ax)
        off = 0.0
        tops = []
        for reg in REG_NAMES[::-1]:
            fd = figdata[reg]
            lam, r_, m_ = cont_sub(fd['Sd'], fd['Id'], gname)
            q = (lam >= lo_) & (lam <= hi_) & (fd['Id'][np.searchsorted(WAVE, lam)] > 0)
            rr = np.where(q, r_, np.nan)
            ax.step(lam, rr + off, where='mid', color=CREG[reg], lw=0.8, zorder=3)
            ax.plot(lam[q], m_[q] + off, color=ss.FG, lw=0.8, zorder=4)
            ax.axhline(off, color=ss.GRID, lw=0.5, zorder=0)
            ax.text(lo_ + 0.03 * (hi_ - lo_), off + (0.10 if gname == 'o2' else 0.035), reg, fontsize=FS['note'],
                    color=CREG[reg], ha='left', va='bottom', path_effects=[ss.HALO])
            tops.append(off + np.nanmax(np.r_[rr[np.isfinite(rr)], m_[q]]))
            off += (0.52 if gname == 'o2' else 0.16)
        for l0, lab in marks:
            ax.axvline(l0, color=ss.FGDIM, lw=0.5, ls=':', zorder=0)
        ax.set_xlim(lo_, hi_)
        ax.set_ylim(-0.17 if gname == 'o2' else -0.075, max(tops) * 1.08 + (0.05 if gname == 'o2' else 0.03))
        ax.xaxis.set_major_locator(plt.MultipleLocator(loc_))
        ax.set_xlabel('vacuum λ (Å)', fontsize=FS['label'], labelpad=1.5)
        ax.set_ylabel('ON−OFF (R Å⁻¹) + offset', fontsize=FS['label'], labelpad=2)
        for l0, lab in marks:
            if lab:
                ax.text(l0, 0.975, lab, transform=ax.get_xaxis_transform(), va='top',
                        ha='right' if lab.endswith(' ') else 'left' if lab.startswith(' ') else 'center',
                        fontsize=FS['note'], color=ss.FG, path_effects=[ss.HALO])
        tag(ax, '(%s) %s' % (next(letters), title))

    # (l) ratios to [O III] 5007: regression (filled) and direct stack (open), P1 and P3
    ax = fig.add_subplot(gs[2, 3])
    style(ax)
    rows = [('OII3727', '[O II]'), ('OIII4959', '4959'), ('Hb', 'Hβ'), ('Ha', 'Hα'), ('NII6583', '[N II]'),
            ('SII6716+6731', '[S II]'), ('HeII4686', 'He II')]
    for iy, (nm, lab) in enumerate(rows):
        yb = len(rows) - 1 - iy
        k = 0
        for reg in REG_NAMES:
            for est, filled_ in (('regression_intercept', True), ('direct_on_minus_off_stack', False)):
                d = summary['regions'][reg]['ratios'][est][nm]
                yy = yb + 0.27 - 0.18 * k
                ax.errorbar([d['ratio']], [yy], xerr=[d['sig_adopt']], fmt='o' if reg == 'P1' else 's', ms=3.0,
                            mfc=CREG[reg] if filled_ else ss.AXBG, mec=CREG[reg], mew=0.7, color=CREG[reg], lw=0.7,
                            zorder=3)
                k += 1
        if iy:
            ax.axhline(yb + 0.5, color=ss.GRID, lw=0.5, zorder=0)
    ax.axvline(0, color=ss.FGDIM, lw=0.5, zorder=1)
    ax.axvline(1 / 2.98, color=ss.FGDIM, lw=0.5, ls=':', zorder=1)
    ax.set_xscale('symlog', linthresh=0.2, linscale=1.0)
    ax.set_xlim(-0.16, 1.6)
    ax.set_xticks([-0.1, 0, 0.1, 0.3, 1.0])
    ax.set_xticklabels(['−0.1', '0', '0.1', '0.3', '1'])
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([lab for nm, lab in rows][::-1], fontsize=FS['tick'])
    ax.set_ylim(-0.55, len(rows) - 0.45)
    ax.tick_params(axis='y', length=0, pad=1.5)
    ax.set_xlabel('ratio to [O III] λ5007', fontsize=FS['label'], labelpad=1.5)
    hl = [Line2D([], [], marker='o', ls='', ms=3, mfc=CREG['P1'], mec=CREG['P1']),
          Line2D([], [], marker='s', ls='', ms=3, mfc=CREG['P3'], mec=CREG['P3']),
          Line2D([], [], marker='o', ls='', ms=3, mfc=ss.FG, mec=ss.FG),
          Line2D([], [], marker='o', ls='', ms=3, mfc=ss.AXBG, mec=ss.FG, mew=0.7)]
    ax.legend(hl, ['P1', 'P3', 'regression', 'ON−OFF'], loc='upper left', fontsize=FS['note'], frameon=True,
              framealpha=0.9, handlelength=0.8, handletextpad=0.3, borderpad=0.3, labelspacing=0.2,
              edgecolor=ss.EDGE).get_frame().set_linewidth(0.5)
    tag(ax, '(%s) line ratios' % next(letters))

    out = ss.save_paper(fig, ROOT / 'figures/paper_appendix_desi.png', dpi=300)
    plt.close(fig)
    print('figure ->', out)


if __name__ == '__main__':
    main()
