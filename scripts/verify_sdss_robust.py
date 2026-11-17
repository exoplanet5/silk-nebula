"""verify_sdss stage 5: systematics and precovery.

(1) Placebo strands: the trace shifted by dRA = +-1.0 ... +-3.5 deg; ON/OFF re-defined around the shifted trace (fibres
    within 0.6 deg of the real strands removed); the scatter of placebo ON-OFF amplitudes per line is an empirical
    systematic floor for each region.
(2) Gradient-robust estimate: per group and line, a_f = c0 + c1 dx + c2 dy + A * ON_f fitted over all fibres with
    dmin < 2 deg (Halpha patch excluded unless it is the target); A combined over groups.
(3) Precovery: per plate-MJD, the [O III] ON-OFF (ON < 0.15 deg) and a plane + Gaussian-profile (sigma 0.12 deg) fit;
    S/N from the OFF-fibre scatter.  Dates from MJD (SDSS MJD = UT date of the night's observations) and the earliest
    MJD in the coadd (MJDLIST).
Usage: python -I verify_sdss_robust.py <diff_tag> <f_sdss> <f_boss>
Writes data/derived/verify_sdss_robust_<tag>.json and verify_sdss_precovery.csv
"""
import sys, csv, json
import numpy as np
from astropy.time import Time
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_sdss_lib import ROOT, load_strands, strand_distance
from verify_sdss_results import load, good_groups, group_stats, combine

LINES = ['OIII:OIII5007', 'Hb:Hb', 'HaNII:Ha', 'HaNII:NII6583', 'SII:SII6716', 'SII:SII6731', 'OII:OII3726', 'OII:OII3729',
         'B4363:OIII4363', 'HeII:HeII4686']


def region_point(fib, amp, gok, on, off, fac):
    st = {L: group_stats(fib, amp, L, on, off, gok) for L in LINES}
    mo, eo, non, noff = st['OIII:OIII5007']
    use = np.isfinite(mo) & np.isfinite(eo) & (eo > 0)
    gids = np.where(use)[0]
    w = np.zeros_like(mo); w[use] = 1 / eo[use] ** 2
    return {L: float(combine(st[L][0] * fac, w, gids)) for L in LINES}, int(len(gids)), int(non[use].sum())


def gradient_robust(fib, amp, gok, on, pool, fac):
    out = {}
    ng = int(fib['gid'].max()) + 1
    for L in LINES:
        a = amp[L]
        A = np.full(ng, np.nan); E = np.full(ng, np.nan)
        for g in np.where(gok)[0]:
            s = (fib['gid'] == g) & pool & np.isfinite(a)
            if (s & on).sum() < 1 or (s & ~on).sum() < 6:
                continue
            x0 = np.mean(fib['ra'][s]); y0 = np.mean(fib['dec'][s])
            dx = (fib['ra'][s] - x0) * np.cos(np.radians(y0)); dy = fib['dec'][s] - y0
            M = np.c_[np.ones(s.sum()), dx, dy, on[s].astype(float)]
            y = a[s]
            c, *_ = np.linalg.lstsq(M, y, rcond=None)
            r = y - M @ c
            sd = 1.4826 * np.median(np.abs(r - np.median(r)))
            cov = np.linalg.pinv(M.T @ M) * sd ** 2
            A[g] = c[3] * fac[g]; E[g] = np.sqrt(cov[3, 3]) * fac[g]
        ok = np.isfinite(A) & np.isfinite(E) & (E > 0)
        if ok.sum() < 2:
            out[L] = None
            continue
        w = 1 / E[ok] ** 2
        out[L] = {'I_R': float(np.sum(w * A[ok]) / np.sum(w)), 'err': float(1 / np.sqrt(np.sum(w))), 'n_groups': int(ok.sum())}
    return out


def main(tag, f_sdss, f_boss):
    fib, grp, amp, Z = load(tag)
    gok, _ = good_groups(grp)
    fac = np.where(grp['instrument'] == 'SDSS', f_sdss, f_boss)
    S, patch, tr = load_strands()
    dec = fib['dec']
    out = {'tag': tag, 'f_sdss': f_sdss, 'f_boss': f_boss}
    # ---------------- placebo
    shifts = [-3.5, -3.0, -2.5, -2.0, -1.5, -1.0, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5]
    pl = {'south': [], 'north': [], 'all': []}
    for dra in shifts:
        dd = []
        for k in ('E', 'W', 'W2'):
            P = S[k].copy(); P[:, 0] += dra
            d, end, _ = strand_distance(fib['ra'], fib['dec'], P)
            dd.append(np.abs(d))
        ds = np.min(dd, 0)
        on = (ds < 0.15) & (fib['dmin'] > 0.6) & ~fib['inpatch30']
        off = (ds > 0.45) & (ds < 2.0) & (fib['dmin'] > 0.45) & ~fib['inpatch30']
        for nm, sel in (('south', dec < 20), ('north', dec >= 20), ('all', np.ones_like(dec, bool))):
            r, ng, non = region_point(fib, amp, gok, on & sel, off & ~(on & sel), fac)
            r['shift'] = dra; r['n_groups'] = ng; r['n_on'] = non
            pl[nm].append(r)
    summ = {}
    for nm, lst in pl.items():
        summ[nm] = {L: {'mean': float(np.nanmean([x[L] for x in lst])), 'rms': float(np.sqrt(np.nanmean(np.array([x[L] for x in lst]) ** 2))),
                        'std': float(np.nanstd([x[L] for x in lst]))} for L in LINES}
        summ[nm]['n_on_median'] = float(np.median([x['n_on'] for x in lst]))
    out['placebo'] = {'shifts_deg': shifts, 'per_shift': pl, 'summary': summ}
    # ---------------- gradient-robust
    ON = fib['dmin'] < 0.15
    pool = (fib['dmin'] < 2.0) & ~fib['inpatch30']
    gr = {}
    for nm, sel in (('south_dec_lt20', dec < 20), ('north_dec_ge20', dec >= 20), ('W_strand_15_20', (fib['nearest'] == 'W') & (dec >= 15) & (dec < 20)),
                    ('all', np.ones_like(dec, bool))):
        # ON fibres outside the region are removed from the pool so they do not enter the plane fit
        gr[nm] = gradient_robust(fib, amp, gok, ON & sel, pool & ~(ON & ~sel), fac)
    poolp = (fib['dmin'] < 2.5)
    gr['halpha_patch_30pct'] = gradient_robust(fib, amp, gok, fib['inpatch30'], poolp & ~(ON & ~fib['inpatch30']), fac)
    out['gradient_robust'] = gr
    # ---------------- precovery per plate-MJD
    a = amp['OIII:OIII5007']
    off = (fib['dmin'] > 0.45) & (fib['dmin'] < 2.0) & ~fib['inpatch30']
    pk = np.array([f"{r}|{int(p)}|{int(m)}" for r, p, m in zip(fib['release'], fib['plate'], fib['mjd'])])
    rows = []
    for key in sorted(set(pk)):
        s = pk == key
        gs = np.unique(fib['gid'][s])
        if not np.all(gok[gs]):
            gs = gs[gok[gs]]
        if len(gs) == 0:
            continue
        vals, errs, nons, Ap, Ep = [], [], [], [], []
        for g in gs:
            ss = (fib['gid'] == g) & np.isfinite(a)
            ao = a[ss & ON]; af = a[ss & off]
            if len(ao) >= 1 and len(af) >= 3:
                sd = 1.4826 * np.median(np.abs(af - np.median(af)))
                vals.append((np.mean(ao) - np.mean(af)) * fac[g]); errs.append(sd * np.sqrt(1 / len(ao) + 1 / len(af)) * fac[g]); nons.append(len(ao))
            # profile fit
            sp_ = ss & (fib['dmin'] < 2.0) & ~fib['inpatch30']
            if sp_.sum() >= 8 and (sp_ & (fib['dmin'] < 0.3)).sum() >= 1:
                x0 = np.mean(fib['ra'][sp_]); y0 = np.mean(fib['dec'][sp_])
                dx = (fib['ra'][sp_] - x0) * np.cos(np.radians(y0)); dy = fib['dec'][sp_] - y0
                prof = np.exp(-0.5 * (fib['dmin'][sp_] / 0.12) ** 2)
                M = np.c_[np.ones(sp_.sum()), dx, dy, prof]
                c, *_ = np.linalg.lstsq(M, a[sp_], rcond=None)
                r = a[sp_] - M @ c
                sd = 1.4826 * np.median(np.abs(r - np.median(r)))
                cov = np.linalg.pinv(M.T @ M) * sd ** 2
                Ap.append(c[3] * fac[g]); Ep.append(np.sqrt(cov[3, 3]) * fac[g])
        if not vals and not Ap:
            continue
        row = dict(release=key.split('|')[0], plate=int(key.split('|')[1]), mjd=int(key.split('|')[2]),
                   mjdfirst=int(np.min(fib['mjdfirst'][s])), instrument=str(grp['instrument'][gs[0]]),
                   ra=float(np.mean(fib['ra'][s])), dec=float(np.mean(fib['dec'][s])))
        if vals:
            w = 1 / np.array(errs) ** 2
            row.update(n_on=int(np.sum(nons)), onoff_R=float(np.sum(w * np.array(vals)) / np.sum(w)), onoff_err=float(1 / np.sqrt(np.sum(w))))
            row['onoff_snr'] = row['onoff_R'] / row['onoff_err']
        else:
            row.update(n_on=0, onoff_R=np.nan, onoff_err=np.nan, onoff_snr=np.nan)
        if Ap:
            w = 1 / np.array(Ep) ** 2
            row.update(prof_R=float(np.sum(w * np.array(Ap)) / np.sum(w)), prof_err=float(1 / np.sqrt(np.sum(w))))
            row['prof_snr'] = row['prof_R'] / row['prof_err']
        else:
            row.update(prof_R=np.nan, prof_err=np.nan, prof_snr=np.nan)
        row['date_ut'] = Time(row['mjd'], format='mjd').iso[:10]
        row['date_first_ut'] = Time(row['mjdfirst'], format='mjd').iso[:10]
        rows.append(row)
    rows.sort(key=lambda r: r['mjdfirst'])
    with open(ROOT / 'data/derived/verify_sdss_precovery.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader()
        for r in rows:
            w.writerow(r)
    det = [r for r in rows if (np.isfinite(r['onoff_snr']) and r['onoff_snr'] >= 3) or (np.isfinite(r['prof_snr']) and r['prof_snr'] >= 3)]
    out['precovery'] = {'n_plates': len(rows), 'n_detected_3sigma': len(det),
                        'detections': [{k: r[k] for k in ('plate', 'mjd', 'mjdfirst', 'date_first_ut', 'instrument', 'ra', 'dec', 'n_on',
                                                         'onoff_R', 'onoff_err', 'onoff_snr', 'prof_R', 'prof_err', 'prof_snr')} for r in det]}
    json.dump(out, open(ROOT / f'data/derived/verify_sdss_robust_{tag}.json', 'w'), indent=1, default=float)
    print('placebo summary (rms, R):')
    for nm in summ:
        print(' ', nm, {L.split(':')[1]: round(summ[nm][L]['rms'], 3) for L in LINES}, 'n_on~', summ[nm]['n_on_median'])
    print('gradient-robust:')
    for nm, v in gr.items():
        print(' ', nm, {L.split(':')[1]: (round(v[L]['I_R'], 3), round(v[L]['err'], 3)) for L in LINES if v.get(L)})
    print('precovery: plates', len(rows), 'detected', len(det))
    for r in det[:12]:
        print('  ', r['plate'], r['mjd'], r['mjdfirst'], r['date_first_ut'], r['instrument'], f"{r['onoff_R']:.2f}+-{r['onoff_err']:.2f} ({r['onoff_snr']:.1f}) prof {r['prof_R']:.2f}+-{r['prof_err']:.2f} ({r['prof_snr']:.1f})")


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], float(a[1]), float(a[2]))
