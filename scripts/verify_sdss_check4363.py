"""verify_sdss: stack-level check of [O III]4363 (and He II 4686, Hbeta) next to the Hg I 4358 street-lamp line.

ON-OFF group spectra (diffuse-corrected R/A) are stacked in the heliocentric frame with [O III]-error weights; Hg I 4358
is topocentric, so it is fitted at the weighted-mean v_corr with its width broadened by the v_corr scatter.
Writes data/derived/verify_sdss_check4363.json
Usage: python -I verify_sdss_check4363.py <tag> f_sdss f_boss
"""
import sys, json
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_sdss_lib import ROOT, VAC, C, load_windows, flam_to_R_per_A
from verify_sdss_results import load, good_groups


def main(tag, f_sdss, f_boss):
    D = load_windows(); lam = D['lam']
    fib, grp, amp, Z = load(tag)
    gok, _ = good_groups(grp)
    fac = np.where(grp['instrument'] == 'SDSS', f_sdss, f_boss)
    a = amp['OIII:OIII5007']
    ON = fib['dmin'] < 0.15
    OFF = (fib['dmin'] > 0.45) & (fib['dmin'] < 2.0) & ~fib['inpatch30']
    out = {}
    for nm, sel in (('south', ON & (fib['dec'] < 20)), ('all', ON), ('north', ON & (fib['dec'] >= 20))):
        m = (lam > 4320) & (lam < 4400)
        lm = lam[m]
        S, W, V, SIG = [], [], [], []
        for g in np.where(gok)[0]:
            s = np.where(fib['gid'] == g)[0]
            so = s[sel[s]]; sf = s[(OFF & ~sel)[s]]
            if len(so) < 1 or len(sf) < 3:
                continue
            X = flam_to_R_per_A(D['tot'][s][:, m].astype(float), lm, fib['diam'][s[0]]) * fac[g]
            X[D['ivar'][s][:, m] <= 0] = np.nan
            d = np.nanmean(X[np.isin(s, so)], 0) - np.nanmean(X[np.isin(s, sf)], 0)
            af = a[sf][np.isfinite(a[sf])]
            if len(af) < 3:
                continue
            sd = 1.4826 * np.median(np.abs(af - np.median(af))) * fac[g] * np.sqrt(1 / len(so) + 1 / len(sf))
            S.append(d); W.append(1 / sd ** 2); V.append(-float(grp['helio_rv'][g]))
            SIG.append(float(np.nanmedian(D['wdisp'][s][:, m])) * 4360 * np.log(10) * 1e-4)
        S = np.array(S); W = np.array(W); V = np.array(V)
        good = np.isfinite(S)
        sp = np.nansum(np.where(good, S * W[:, None], 0), 0) / np.sum(np.where(good, W[:, None], 0), 0)
        sig = float(np.sum(W * np.array(SIG)) / np.sum(W))
        vm = float(np.sum(W * V) / np.sum(W)); vs = float(np.sqrt(np.sum(W * (V - vm) ** 2) / np.sum(W)))
        sHg = np.sqrt(sig ** 2 + (4360 * vs / C) ** 2)
        G = lambda mu, s_: np.exp(-0.5 * ((lm - mu) / s_) ** 2) / (np.sqrt(2 * np.pi) * s_)
        cols = [np.ones_like(lm), (lm - 4360) / 40, G(VAC['Hgam'] * (1 + 13 / C), sig), G(VAC['Hg4358'] * (1 + vm / C), sHg),
                G(VAC['OIII4363'] * (1 + 13 / C), sig)]
        A = np.array(cols).T
        ok = np.isfinite(sp)
        # noise from residual scatter
        c, *_ = np.linalg.lstsq(A[ok], sp[ok], rcond=None)
        r = sp[ok] - A[ok] @ c
        cov = np.linalg.pinv(A[ok].T @ A[ok]) * np.var(r)
        # group bootstrap
        rng = np.random.default_rng(3)
        bs = []
        for _ in range(500):
            ii = rng.integers(0, len(W), len(W))
            gg = np.isfinite(S[ii])
            spb = np.nansum(np.where(gg, S[ii] * W[ii][:, None], 0), 0) / np.sum(np.where(gg, W[ii][:, None], 0), 0)
            okb = np.isfinite(spb)
            cb, *_ = np.linalg.lstsq(A[okb], spb[okb], rcond=None)
            bs.append(cb[2:])
        bs = np.array(bs)
        out[nm] = {'Hgam': [float(c[2]), float(np.std(bs[:, 0]))], 'Hg4358': [float(c[3]), float(np.std(bs[:, 1]))],
                   'OIII4363': [float(c[4]), float(np.std(bs[:, 2]))], 'n_groups': int(len(W)), 'vcorr_mean': vm, 'vcorr_rms': vs}
        print(nm, out[nm])
    json.dump(out, open(ROOT / 'data/derived/verify_sdss_check4363.json', 'w'), indent=1)


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], float(a[1]), float(a[2]))
