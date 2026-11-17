"""spectra_ism: empirical telluric H2O correction of the Na I D region (5880-5905 A) for the archival spectra.

For a target spectrum, the telluric template is the median optical depth of all OTHER spectra in the
topocentric frame (each high-pass normalised with a 2 A running median, smoothed to the target resolution and
scaled to a common telluric strength), so interstellar lines -- fixed in the heliocentric frame and spread over
~60 km/s of BERV in the topocentric frame -- are suppressed in the median. The template is scaled to the target by
least squares outside the target's interstellar window (V_helio -40..+70 km/s around D2 and D1) and divided out.
"""
import numpy as np
from scipy.ndimage import median_filter, gaussian_filter1d

C = 299792.458
L1, L2, DL = 5878.0, 5907.0, 0.01
LG = np.arange(L1, L2, DL)
D2, D1 = 5889.951, 5895.924


def topo_wave(sp):
    w = sp['wave']
    return w if sp['frame'] == 'topo' else w / (1 + sp['berv'] / C)


def tau_on_grid(sp, R_target=None):
    w = topo_wave(sp)
    m = (w > L1 - 3) & (w < L2 + 3)
    if m.sum() < 100:
        return None
    f = np.interp(LG, w[m], sp['flux'][m], left=np.nan, right=np.nan)
    if np.isnan(f).mean() > 0.05:
        return None
    f = np.where(np.isfinite(f), f, np.nanmedian(f))
    cont = median_filter(f, size=int(2.0 / DL) | 1, mode='nearest')
    x = np.clip(f / cont, 0.02, 1.5)
    if R_target is not None and R_target < sp['R']:
        lam = 0.5 * (L1 + L2)
        s = np.sqrt((lam / R_target) ** 2 - (lam / sp['R']) ** 2) / 2.3548 / DL
        x = gaussian_filter1d(x, s, mode='nearest')
    return -np.log(x)


def is_mask(sp):
    """True where the heliocentric interstellar window of THIS star falls (in topocentric wavelengths)."""
    shift = 0.0 if sp['frame'] == 'topo' else 0.0
    m = np.zeros_like(LG, bool)
    for lam0 in (D2, D1):
        v_helio_of_grid = (LG / lam0 - 1) * C + sp['berv']   # topocentric grid -> helio velocity
        m |= (v_helio_of_grid > -40 + shift) & (v_helio_of_grid < 70)
    return m


def correct(target, pool):
    """Return (wave_topo_grid, corrected_flux_ratio, alpha, template_tau) for the target spectrum, or None."""
    tt = tau_on_grid(target)
    if tt is None:
        return None
    taus = []
    for sp in pool:
        if sp is target:
            continue
        t = tau_on_grid(sp, R_target=target['R'])
        if t is not None:
            mk = ~is_mask(sp)
            taus.append((t, mk))
    if len(taus) < 4:
        return None
    # common strength: scale each tau to the median of the strongest telluric pixels
    ref = np.nanmedian(np.array([t for t, _ in taus]), 0)
    strong = ref > np.percentile(ref, 97)
    scaled = []
    for t, mk in taus:
        mm = strong & mk
        s = np.sum(t[mm] * ref[mm]) / np.sum(ref[mm] ** 2) if mm.sum() > 5 else np.nan
        if np.isfinite(s) and s > 0.05:
            tt_ = np.where(mk, t / s, np.nan)          # blank each star's own IS window
            scaled.append(tt_)
    tmpl = np.nanmedian(np.array(scaled), 0)
    tmpl = np.where(np.isfinite(tmpl), tmpl, 0.0)
    tmpl = np.clip(tmpl - np.median(tmpl), 0, None)
    out = ~is_mask(target)
    alpha = np.sum(tt[out] * tmpl[out]) / np.sum(tmpl[out] ** 2) if np.sum(tmpl[out] ** 2) > 0 else 0.0
    alpha = max(alpha, 0.0)
    return LG, np.exp(alpha * tmpl), alpha, tmpl


def apply(target, corr):
    """Divide the target's native flux by the telluric transmission (interpolated to its wavelengths)."""
    lg, trans_inv, alpha, _ = corr
    w = topo_wave(target)
    f = target['flux'].copy()
    m = (w > L1) & (w < L2)
    f[m] = f[m] * np.interp(w[m], lg, trans_inv)
    out = dict(target)
    out['flux'] = f
    return out
