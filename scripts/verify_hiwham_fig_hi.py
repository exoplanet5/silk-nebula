"""verify_hiwham: figure of the GALFA-HI cross-spine tests at the revised [O III] velocity.

Panels
 (a) whole-spine cross-profile 'PV' image: N(HI) per 2-km/s LSR slice vs offset d, each row minus its mean,
     divided by the coherent-null robust sigma of that row's binned values (i.e. in null-sigma units)
 (b) same for the south (spine Dec <= 17.6)
 (c) coherent-null z of the filament amplitude A_fil and the local step S_loc vs velocity (whole, south, north)
 (d) profiles in the revised-velocity window (V_helio +10.5..+16.5) and in the wave-1 window (V_LSR +11..+16),
     with the 16-84 % envelope of the coherent nulls
 (e) wave-1 'wall' finder: per-segment offsets (+11..+16) vs the single-segment null distribution, and the
     11-segment mean vs its null distribution
Output: figures/verify_hiwham_hi.png
"""
import sys

import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import new_sky_figure, finish, label, OIII_C, HA_C, FGDIM, WARN, OK, ACCENT, FG
from verify_hiwham_common import DER, FIG, strand, Trace, galfa_cube, nhi_slice, vlsr_minus_vhelio

Z = np.load(DER / 'verify_hiwham_hi_profiles.npz')
names = list(Z['names']); groups = list(Z['groups']); STATS = list(Z['stats'])
CEN = Z['cen']; P_real = Z['P_real']; STr = Z['ST_real']; STc = Z['ST_coh']; STs = Z['ST_seg']; seg_of = Z['seg_of']
slices = Z['slices']
nsl = len(slices)
vc = slices.mean(axis=1)

# ---------------------------------------------------------- recompute coherent null profiles (group 'all', 'south')
cube, v, gra, gdec = galfa_cube()
RA, DEC = np.meshgrid(gra, gdec); pra, pdec = RA.ravel(), DEC.ravel()
import importlib.util
spec = importlib.util.spec_from_file_location('hi', '/Users/mickey/astronomy/silk-nebula/scripts/verify_hiwham_hi.py')
# (we do not re-run the analysis module; profiles are rebuilt here with the same binning)
from scipy import sparse
maps = [nhi_slice(cube, v, lo, hi) for lo, hi in slices]
ny, nx = RA.shape
corr = vlsr_minus_vhelio(pra, pdec).reshape(ny, nx)


def helio_window(lo_h, hi_h):
    dv = abs(v[1] - v[0]); e0, e1 = v - dv / 2, v + dv / 2
    lo = lo_h + corr; hi = hi_h + corr
    o = np.zeros((ny, nx)); bad = np.zeros((ny, nx), bool)
    for k in range(len(v)):
        w = np.clip(np.minimum(e1[k], hi) - np.maximum(e0[k], lo), 0, None)
        if np.any(w > 0):
            o += w * np.nan_to_num(cube[k]); bad |= (w > 0) & ~np.isfinite(cube[k])
    o *= 1.823e18; o[bad] = np.nan
    return o


maps += [helio_window(10.5, 16.5), nhi_slice(cube, v, 11.0, 16.0)]
M = np.array([m.ravel() for m in maps]).T
F = np.isfinite(M).astype(float); M0 = np.nan_to_num(M)
ra_E, dec_E = strand('E')
REAL = Trace(ra_E, dec_E)
s_split = float(Z['s_split'])
EDGES = np.arange(-3.0, 3.0 + 1e-9, 0.05)


def prof2(tr):
    d, s, ins = tr.coords(pra, pdec, dmax=3.0)
    m = ins & np.isfinite(d) & (np.abs(d) < 3.0)
    b = ((d[m] + 3) / 0.05).astype(int); idx = np.where(m)[0]; sm = s[m]
    rows = np.concatenate([b, len(CEN) + b[sm <= s_split]])
    cols = np.concatenate([idx, idx[sm <= s_split]])
    A = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(2 * len(CEN), len(pra)))
    S = A @ M0; N = A @ F
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(N >= 4, S / N, np.nan).reshape(2, len(CEN), -1)


Pr = prof2(REAL)
pl = Z['coh_placements']
Pn = np.array([prof2(Trace(ra_E, dec_E, rot=p[0], mirror=bool(p[1]), centre=(p[2], p[3]))) for p in pl])
print('null profiles', Pn.shape)


def norm_rows(P, Pn, g):
    """row-wise: (profile - its mean over |d|<3) / robust sigma of the same quantity over nulls (per bin)."""
    x = P[g] - np.nanmean(P[g], axis=0, keepdims=True)
    xn = Pn[:, g] - np.nanmean(Pn[:, g], axis=1, keepdims=True)
    sig = 0.5 * (np.nanpercentile(xn, 84, axis=0) - np.nanpercentile(xn, 16, axis=0))
    med = np.nanmedian(xn, axis=0)
    return ((x - med) / sig).T, x, med, sig     # (nmap, nbin)


fig, axes = new_sky_figure(15, 11.5, 'Silk Nebula · GALFA-HI across the spine at the revised [O III] velocity',
                           'N(HI) cross-profiles vs 120 coherent and 1500 segment trace-placement nulls; d > 0 = east (outside)',
                           nrows=2, ncols=3)
ax = axes.ravel()
v_true = (13.5 + vlsr_minus_vhelio(116.6, 15.0)[0], 13.5 + vlsr_minus_vhelio(115.5, 25.0)[0])
for k, (g, ttl) in enumerate(((0, 'whole spine'), (1, 'south (Dec ≤ 17.6)'))):
    zimg, x, med, sig = norm_rows(Pr, Pn, g)
    im = ax[k].imshow(zimg[:nsl], origin='lower', aspect='auto', cmap='RdBu_r', vmin=-3, vmax=3,
                      extent=[-3, 3, slices[0, 0], slices[-1, 1]], interpolation='nearest')
    ax[k].axvline(0, color='#ffffff', lw=0.7, alpha=0.7)
    ax[k].axhspan(v_true[0] - 1.5, v_true[1] + 1.5, color=OIII_C, alpha=0.0, lw=1.2, ec=OIII_C, ls='-')
    ax[k].axhline(v_true[0], color=OIII_C, lw=1.0); ax[k].axhline(v_true[1], color=OIII_C, lw=1.0, ls='-')
    ax[k].axhspan(11, 16, color=WARN, alpha=0.0, ec=WARN, lw=1.0)
    label(ax[k], -2.9, v_true[1] + 0.8, '[O III] revised', color=OIII_C, size=7.5)
    label(ax[k], -2.9, 16.4, 'wave-1 window', color=WARN, size=7.5)
    ax[k].axvline(0.86, color=WARN, lw=0.8, ls=':')
    ax[k].set_xlabel('offset d from E spine (deg, + = east/outside)'); ax[k].set_ylabel('V_LSR (km/s)')
    ax[k].set_title(f'({"ab"[k]}) {ttl}: row-normalised profile / null σ', loc='left')
    cb = fig.colorbar(im, ax=ax[k], fraction=0.04, pad=0.01); cb.set_label('σ (coherent nulls)', fontsize=7.5)

# (c) z vs velocity
kA, kS = STATS.index('A_fil'), STATS.index('S_loc')
for g, col, ls in (('all', FG, '-'), ('south', OIII_C, '-'), ('north', WARN, '-')):
    gi = groups.index(g)
    zA, zS = [], []
    for j in range(nsl):
        for kk, store in ((kA, zA), (kS, zS)):
            nl = STc[:, gi, j, kk]; nl = nl[np.isfinite(nl)]
            sr = 0.5 * (np.percentile(nl, 84) - np.percentile(nl, 16))
            store.append((STr[gi, j, kk] - np.median(nl)) / sr)
    ax[2].plot(vc, zA, color=col, lw=1.4, ls='-', marker='o', ms=3, label=f'A_fil {g}')
    ax[2].plot(vc, zS, color=col, lw=1.0, ls='--', label=f'S_loc {g}')
ax[2].axhspan(-2, 2, color='#2a323c', alpha=0.5, lw=0)
ax[2].axvspan(*v_true, color=OIII_C, alpha=0.18, lw=0)
ax[2].axvspan(11, 16, color=WARN, alpha=0.12, lw=0)
ax[2].axvspan(16 - 12.0 + 0.5, 21 - 12.0 + 0.5, color=OK, alpha=0.10, lw=0)
label(ax[2], v_true[0] - 0.5, -2.6, 'revised\n[O III]', color=OIII_C, size=7.5)
label(ax[2], 11.2, -2.6, 'wave-1', color=WARN, size=7.5)
label(ax[2], 4.6, -3.5, 'W strand', color=OK, size=7.5)
ax[2].set_ylim(-3.8, 3.8); ax[2].set_xlabel('slice centre V_LSR (km/s)'); ax[2].set_ylabel('z vs coherent nulls')
ax[2].set_title('(c) z: A_fil (solid), S_loc (dashed)', loc='left')
ax[2].legend(ncol=2, fontsize=6.5, loc='upper left')

# (d) profiles with null envelopes
for idx, (j, col, lab) in enumerate(((nsl, OIII_C, 'V_helio +10.5…+16.5 (revised)'), (nsl + 1, WARN, 'V_LSR +11…+16 (wave-1)'))):
    for g, lw in ((0, 1.6),):
        x = Pr[g, :, j] - np.nanmean(Pr[g, :, j])
        xn = Pn[:, g, :, j] - np.nanmean(Pn[:, g, :, j], axis=1, keepdims=True)
        lo, hi = np.nanpercentile(xn, [16, 84], axis=0)
        ax[3].fill_between(CEN, lo / 1e19, hi / 1e19, color=col, alpha=0.15, lw=0)
        ax[3].plot(CEN, x / 1e19, color=col, lw=lw, label=lab)
ax[3].axvline(0, color='#ffffff', lw=0.7); ax[3].axvline(0.86, color=WARN, lw=0.8, ls=':')
label(ax[3], 0.9, ax[3].get_ylim()[1] * 0.85 if ax[3].get_ylim()[1] > 0 else 1, 'wave-1 wall', color=WARN, size=7)
ax[3].set_xlabel('offset d (deg)'); ax[3].set_ylabel('N(HI) − mean (10¹⁹ cm⁻²)')
ax[3].set_title('(d) whole-spine profiles with 16–84 % coherent-null band', loc='left')
ax[3].legend(fontsize=7, loc='lower left')

# (e) wave-1 wall finder vs nulls
j = names.index('old_oiii'); kw = STATS.index('w1_wall')
nl = STs[:, j, kw]; nl = nl[np.isfinite(nl)]
ax[4].hist(nl, bins=np.arange(-1.0, 3.01, 0.1), color=FGDIM, alpha=0.6, density=True, label='single random segment (null)')
rw = STr[:11, j, kw]
for x in rw[np.isfinite(rw)]:
    ax[4].axvline(x, color=WARN, lw=1.2, alpha=0.8)
rng = np.random.default_rng(3)
by = [np.where(seg_of == i)[0] for i in range(11)]
means = np.array([np.nanmean(STs[[rng.choice(b) for b in by], j, kw]) for _ in range(4000)])
ax4b = ax[4].twinx()
ax4b.hist(means, bins=np.arange(-0.5, 2.01, 0.05), histtype='step', color=ACCENT, lw=1.4, density=True)
ax4b.axvline(np.nanmean(rw), color=WARN, lw=2.2, ls='--')
ax4b.set_yticks([])
label(ax[4], 1.85, ax[4].get_ylim()[1] * 0.75, 'orange: real segments\ndashed: real 11-seg mean\nblue: null 11-seg mean', color=FG, size=7)
ax[4].set_xlabel('wave-1 "wall" offset (deg), V_LSR +11…+16'); ax[4].set_ylabel('density')
ax[4].set_title('(e) wave-1 wall finder returns ~0.8° on random sky', loc='left')
ax[4].legend(fontsize=6.5, loc='upper left')

# (f) 2-sigma upper limits on an HI filament vs velocity (whole spine, coherent nulls)
gi = groups.index('all')
ul = []
for jj in range(nsl):
    nl = STc[:, gi, jj, kA]; nl = nl[np.isfinite(nl)]
    sr = 0.5 * (np.percentile(nl, 84) - np.percentile(nl, 16))
    ul.append(max(STr[gi, jj, kA] - np.median(nl), 0) + 2 * sr)
gi2 = groups.index('south')
ul2 = []
for jj in range(nsl):
    nl = STc[:, gi2, jj, kA]; nl = nl[np.isfinite(nl)]
    sr = 0.5 * (np.percentile(nl, 84) - np.percentile(nl, 16))
    ul2.append(max(STr[gi2, jj, kA] - np.median(nl), 0) + 2 * sr)
ax[5].step(vc, np.array(ul) / 1e19, where='mid', color=FG, lw=1.5, label='whole spine')
ax[5].step(vc, np.array(ul2) / 1e19, where='mid', color=OIII_C, lw=1.2, label='south')
ax[5].axvspan(*v_true, color=OIII_C, alpha=0.18, lw=0)
ax[5].set_xlabel('slice centre V_LSR (km/s, 2 km/s slices)'); ax[5].set_ylabel('2σ UL ΔN(HI) peak (10¹⁹ cm⁻²)')
ax[5].set_title('(f) limit on an HI filament on the spine', loc='left')
ax[5].legend(fontsize=7)
finish(fig, FIG / 'verify_hiwham_hi.png')
print('saved')
