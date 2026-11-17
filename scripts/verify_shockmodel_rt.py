"""verify_shockmodel_rt: Monte-Carlo estimate of resonance-line (C IV 1548/1551) suppression for an edge-on sheet.

Geometry: a ribbon of thickness D (x in [0, D]) and line-of-sight length L = G*D (z in [-L/2, L/2]), infinite along the
filament (y).  Uniform emissivity.  Pure resonance scattering (no destruction: no dust worth mentioning at E(B-V) < 0.01
in front of / inside the sheet), isotropic redistribution in angle, complete redistribution over a Gaussian Doppler
profile (phi(x) = exp(-x^2), x in Doppler units), line-centre opacity kappa0 = tau_face / D.
Optional linear velocity gradient along the line of sight, w(z) = a*z (Doppler units per unit length), representing the
projected radial flow of a curved (radius R_c) sheet: a = u / (R_c b).
Estimator: peel-off toward +z at every emission and scattering event; S = <sum exp(-tau_peel)> per emitted photon is the
edge-on intensity relative to the optically thin case (S = 1 for kappa0 -> 0).
"""
import numpy as np
from scipy.special import erf, erfinv

SQPI = np.sqrt(np.pi)


def _tau_path(k0, c, a_n, l):
    """Optical depth along a path of length l from a point where (x_f - a z) = c, moving with d(az)/dl = a_n."""
    small = np.abs(a_n) < 1e-9
    out = np.empty_like(l)
    out[small] = k0 * np.exp(-c[small]**2) * l[small]
    s = ~small
    out[s] = k0 * SQPI / (2 * a_n[s]) * (erf(c[s]) - erf(c[s] - a_n[s] * l[s]))
    return out


def edge_on_suppression(tau_face, G, agrad_total=0.0, nphot=20000, seed=1, maxscat=4000):
    """tau_face: line-centre optical depth across the thickness D (=1 unit); G = L/D; agrad_total = total LOS velocity
    span across L in Doppler units (a = agrad_total / L).  Returns S (edge-on intensity / optically thin value)."""
    rng = np.random.default_rng(seed)
    D, L = 1.0, float(G)
    k0 = tau_face / D
    a = agrad_total / L
    x = rng.uniform(0, D, nphot)
    z = rng.uniform(-L / 2, L / 2, nphot)
    xf = a * z + rng.normal(0, 1 / np.sqrt(2), nphot)        # emitted in the local gas frame
    tally = np.zeros(nphot)
    alive = np.ones(nphot, bool)
    for it in range(maxscat):
        idx = np.nonzero(alive)[0]
        if idx.size == 0:
            break
        # peel-off toward +z
        c = xf[idx] - a * z[idx]
        lz = L / 2 - z[idx]
        tally[idx] += np.exp(-_tau_path(k0, c, np.full(idx.size, a), lz))
        # isotropic direction
        mu = rng.uniform(-1, 1, idx.size)
        phi = rng.uniform(0, 2 * np.pi, idx.size)
        st = np.sqrt(1 - mu**2)
        nx, nz = st * np.cos(phi), mu
        # distance to the boundary
        with np.errstate(divide='ignore'):
            lx = np.where(nx > 0, (D - x[idx]) / nx, np.where(nx < 0, -x[idx] / nx, np.inf))
            lzb = np.where(nz > 0, (L / 2 - z[idx]) / nz, np.where(nz < 0, (-L / 2 - z[idx]) / nz, np.inf))
        lb = np.minimum(lx, lzb)
        # sample interaction distance
        td = -np.log(rng.uniform(size=idx.size))
        a_n = a * nz
        tau_b = _tau_path(k0, c, a_n, np.where(np.isfinite(lb), lb, 1e6))
        esc = td >= tau_b
        alive[idx[esc]] = False
        keep = ~esc
        ii = idx[keep]
        if ii.size == 0:
            continue
        cc, an, tdk = c[keep], a_n[keep], td[keep]
        smallg = np.abs(an) < 1e-9
        l = np.empty(ii.size)
        l[smallg] = tdk[smallg] / (k0 * np.exp(-cc[smallg]**2))
        g = ~smallg
        arg = erf(cc[g]) - 2 * an[g] * tdk[g] / (k0 * SQPI)
        arg = np.clip(arg, -1 + 1e-15, 1 - 1e-15)
        l[g] = (cc[g] - erfinv(arg)) / an[g]
        l = np.minimum(l, lb[keep] * (1 - 1e-9))
        x[ii] += nx[keep] * l
        z[ii] += nz[keep] * l
        xf[ii] = a * z[ii] + rng.normal(0, 1 / np.sqrt(2), ii.size)   # CRD re-emission in the local frame
    return float(tally.mean())


if __name__ == '__main__':
    for tf in (0.01, 0.1, 0.3, 1.0, 3.0):
        for G in (10, 30):
            for ag in (0.0, 3.0):
                print(f'tau_face={tf:5.2f} G={G:3d} grad={ag:3.1f}  S={edge_on_suppression(tf, G, ag, nphot=6000):.3f}')
