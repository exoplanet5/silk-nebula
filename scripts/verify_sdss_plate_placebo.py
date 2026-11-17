"""verify_sdss: per-plate placebo test for the precovery table.

For each plate-MJD, the [O III] ON-OFF is recomputed around the trace shifted by dRA = +-0.6, 0.8, 1.0, 1.2 deg (ON within
0.15 deg of the shifted trace and > 0.45 deg from the real strands; OFF 0.45-2 deg from both).  The rms of these placebo
values is a per-plate systematic; S/N_placebo = ON-OFF / sqrt(err^2 + rms^2).  Updates nothing on disk except
data/derived/verify_sdss_precovery_placebo.csv.
Usage: python -I verify_sdss_plate_placebo.py <tag> f_sdss f_boss
"""
import sys, csv
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_sdss_lib import ROOT, load_strands, strand_distance
from verify_sdss_results import load, good_groups


def onoff(a, gid, groups, on, off, fac):
    vals, errs = [], []
    for g in groups:
        ss = (gid == g) & np.isfinite(a)
        ao = a[ss & on]; af = a[ss & off]
        if len(ao) >= 1 and len(af) >= 3:
            sd = 1.4826 * np.median(np.abs(af - np.median(af)))
            vals.append((np.mean(ao) - np.mean(af)) * fac[g]); errs.append(sd * np.sqrt(1 / len(ao) + 1 / len(af)) * fac[g])
    if not vals:
        return np.nan, np.nan
    w = 1 / np.array(errs) ** 2
    return float(np.sum(w * np.array(vals)) / np.sum(w)), float(1 / np.sqrt(np.sum(w)))


def main(tag, f_sdss, f_boss):
    fib, grp, amp, Z = load(tag)
    gok, _ = good_groups(grp)
    fac = np.where(grp['instrument'] == 'SDSS', f_sdss, f_boss)
    S, patch, tr = load_strands()
    a = amp['OIII:OIII5007']
    shifts = [-1.2, -1.0, -0.8, -0.6, 0.6, 0.8, 1.0, 1.2]
    dsh = {}
    for dra in shifts:
        dd = []
        for k in ('E', 'W', 'W2'):
            P = S[k].copy(); P[:, 0] += dra
            d, _, _ = strand_distance(fib['ra'], fib['dec'], P)
            dd.append(np.abs(d))
        dsh[dra] = np.min(dd, 0)
    pk = np.array([f"{r}|{int(p)}|{int(m)}" for r, p, m in zip(fib['release'], fib['plate'], fib['mjd'])])
    ON = fib['dmin'] < 0.15
    OFF = (fib['dmin'] > 0.45) & (fib['dmin'] < 2.0) & ~fib['inpatch30']
    rows = []
    for key in sorted(set(pk)):
        s = pk == key
        gs = [g for g in np.unique(fib['gid'][s]) if gok[g]]
        if not gs:
            continue
        v, e = onoff(a, fib['gid'], gs, ON, OFF, fac)
        if not np.isfinite(v):
            continue
        pv = []
        for dra in shifts:
            on_p = (dsh[dra] < 0.15) & (fib['dmin'] > 0.45) & ~fib['inpatch30']
            off_p = (dsh[dra] > 0.45) & (dsh[dra] < 2.0) & (fib['dmin'] > 0.45) & ~fib['inpatch30']
            x, _ = onoff(a, fib['gid'], gs, on_p, off_p, fac)
            if np.isfinite(x):
                pv.append(x)
        rms = float(np.sqrt(np.mean(np.square(pv)))) if len(pv) >= 3 else np.nan
        rows.append(dict(plate=int(key.split('|')[1]), mjd=int(key.split('|')[2]), mjdfirst=int(np.min(fib['mjdfirst'][s])),
                         onoff_R=v, onoff_err=e, n_placebo=len(pv), placebo_rms_R=rms,
                         snr_placebo=v / np.sqrt(e ** 2 + rms ** 2) if np.isfinite(rms) else np.nan))
    rows.sort(key=lambda r: r['mjdfirst'])
    with open(ROOT / 'data/derived/verify_sdss_precovery_placebo.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader()
        for r in rows:
            w.writerow(r)
    for r in rows:
        print(r['plate'], r['mjd'], r['mjdfirst'], f"{r['onoff_R']:+.2f}±{r['onoff_err']:.2f} placebo rms {r['placebo_rms_R']:.2f} (n={r['n_placebo']}) S/N_pl {r['snr_placebo']:.1f}")


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], float(a[1]), float(a[2]))
