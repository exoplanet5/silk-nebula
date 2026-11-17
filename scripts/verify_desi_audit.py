"""verify_desi_audit.py -- second-pass (resume) audit of the wave-2 DESI verification (tag desi).

Independent of the Gaussian line fitter in verify_desi_linefit.py: line fluxes here are BOX INTEGRALS of the
continuum-subtracted residual (flux - Redrock model), with an empirical sky-line mask built from the median ivar
spectrum. Re-uses only the per-coadd metadata (geometry template Z, chromatic fibre acceptance f(lambda)) from
data/derived/verify_desi_coadd_meta.npz.

Tests:
  1. IV ratios rho_L = sum u z (L - <L>) / sum u z (O - <O>), z = Z - <Z>, HEALPix block bootstrap (nside 16 & 32);
     empirical null boxes (line-free offsets) to calibrate the bootstrap.
  2. Sign split Z > 0 (on-ridge) vs Z < 0 (fibres whose petal sky sampled the ridge) for Hb, Ha, [N II], [S II].
  3. [O III]4363 vs Hg I 4358: pixel IV regression spectrum 4340-4390 A, free-centre two-Gaussian fit, bootstrap.
  4. [O III] velocity, E vs W strand: joint pixel regression on Z_E and Z_W around 5008 A, Gaussian centroid.
  5. Calibration bookkeeping: GFA fibre fraction by program; ratio to MEAN_PSF_TO_FIBER_SPECFLUX.
Outputs: data/derived/verify_desi_audit.json, data/derived/verify_desi_audit_regspec.npz
"""
import csv
import json
import sys

import numpy as np
from scipy.optimize import least_squares
from scipy.special import erf

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_desi_common import D, ROOT, WAVE, C, L, GAL_REST, load_spectra, lsf_sigma  # noqa: E402

RNG = np.random.default_rng(20261116)
V0 = 13.0                       # km/s, nominal nebular velocity for box centres (boxes are +-2.4 A = +-145 km/s)
DLAM = 0.8
HW = 2.4


def healpix_ring(ra, dec, nside):
    import healpy as hp
    return hp.ang2pix(nside, ra, dec, lonlat=True)


# ----------------------------------------------------------------------------- empirical sky-line mask
def sky_mask(I):
    s = 1.0 / np.nanmedian(np.where(I > 0, I, np.nan), 0)
    loc = np.array([np.nanmedian(s[max(0, k - 12):k + 13]) for k in range(s.size)])
    return s / loc > 1.25, s / loc


# ----------------------------------------------------------------------------- box fluxes
def box_flux(R, I, z, lam_lo, lam_hi, other_lines, skym, side=(6.0, 30.0), shift=None):
    """Box integral of R over [lam_lo, lam_hi] (fractional pixel overlap) minus a weighted linear continuum fitted
    to side bands 6-30 A from the box edges. Masks: pixels with ivar<=0, empirical sky lines, other nebular lines
    +-3 A, the object's own redshifted lines +-8 A (box itself invalid if one falls within +-8 A).
    shift (per fibre, in A) moves the box (used for topocentric sky lines)."""
    n = R.shape[0]
    out = np.full(n, np.nan)
    err = np.full(n, np.nan)
    edges_lo = WAVE - DLAM / 2
    edges_hi = WAVE + DLAM / 2
    for i in range(n):
        lo, hi = lam_lo, lam_hi
        if shift is not None:
            lo, hi = lo + shift[i], hi + shift[i]
        zl = GAL_REST * (1 + z[i])
        if np.any((zl > lo - 8) & (zl < hi + 8)):
            continue
        w = np.clip(np.minimum(edges_hi, hi) - np.maximum(edges_lo, lo), 0, None) / DLAM
        inb = w > 0
        sb = (((WAVE > lo - side[1]) & (WAVE < lo - side[0])) | ((WAVE > hi + side[0]) & (WAVE < hi + side[1])))
        sb &= I[i] > 0
        sb &= ~skym
        for lw in other_lines:
            sb &= np.abs(WAVE - lw) > 3.0
        for lw in zl:
            sb &= np.abs(WAVE - lw) > 8.0
        if sb.sum() < 15 or np.any(I[i, inb] <= 0):
            continue
        x = WAVE[sb] - 0.5 * (lo + hi)
        y = R[i, sb]
        wt = I[i, sb]
        keep = np.ones(x.size, bool)
        p = None
        for _ in range(3):
            A = np.c_[np.ones(keep.sum()), x[keep]]
            Wm = wt[keep]
            try:
                cov = np.linalg.inv(A.T @ (A * Wm[:, None]))
            except np.linalg.LinAlgError:
                p = None
                break
            p = cov @ (A.T @ (Wm * y[keep]))
            res = (y - (p[0] + p[1] * x)) * np.sqrt(wt)
            newk = np.abs(res) < 4.0
            if (newk == keep).all() or newk.sum() < 12:
                break
            keep = newk
        if p is None:
            continue
        xb = WAVE[inb] - 0.5 * (lo + hi)
        cb = p[0] + p[1] * xb
        wb = w[inb]
        out[i] = np.sum(wb * (R[i, inb] - cb)) * DLAM
        g = np.array([wb.sum(), np.sum(wb * xb)]) * DLAM
        err[i] = np.sqrt(np.sum(wb ** 2 / I[i, inb]) * DLAM ** 2 + g @ cov @ g)
    return out, err


def vac(lam_rest, v=V0):
    return lam_rest * (1 + v / C)


def box(lam_rest, hw=HW, v=V0):
    c = vac(lam_rest, v)
    return c - hw, c + hw


def build_lines():
    """name: (lo, hi, rest-wavelength for acceptance, kind)."""
    Ls = {}
    for nm in ('Hb', 'OIII4959', 'OIII5007', 'HeII4686', 'NeIII3869', 'Ha', 'NII6583', 'NII6548', 'SII6716',
               'SII6731', 'ArIII7136', 'SIII9069', 'NII5755', 'HeI5876', 'Hg', 'Hd'):
        lo, hi = box(L[nm])
        Ls[nm] = (lo, hi, L[nm], 'neb')
    lo = vac(L['OII3726']) - HW
    hi = vac(L['OII3729']) + HW
    Ls['OII3727sum'] = (lo, hi, 0.5 * (L['OII3726'] + L['OII3729']), 'neb')
    lo = vac(L['OII7320']) - HW
    hi = vac(L['OII7330']) + HW
    Ls['OII7325sum'] = (lo, hi, 7325.0, 'neb')
    lo, hi = box(L['OIII4363'], hw=2.0)
    Ls['OIII4363'] = (lo, hi, L['OIII4363'], 'neb')
    # null boxes (line-free offsets), same half width
    for nm, off in (('OIII5007', 14.0), ('OIII5007', -14.0), ('Hb', 11.0), ('Hb', -11.0), ('OIII4363', 9.0),
                    ('OIII4363', 17.0), ('HeII4686', 10.0), ('NII6583', 11.0), ('SII6731', 12.0),
                    ('NeIII3869', 12.0), ('OII3726', -12.0)):
        c = vac(L[nm]) + off
        Ls[f'null_{nm}_{off:+.0f}'] = (c - HW, c + HW, L[nm] + off, 'null')
    return Ls


def accept_interp(meta, lam):
    """Chromatic point-source acceptance f(lam) per coadd, log-lambda interpolation of the stored values."""
    keys = [k for k in meta.files if k.startswith('f_chrom_') and k[8:].isdigit()]
    lams = np.array([float(k[8:]) for k in keys])
    o = np.argsort(lams)
    lams = lams[o]
    F = np.array([meta[keys[j]] for j in o])
    return np.array([np.interp(np.log(lam), np.log(lams), F[:, i]) for i in range(F.shape[1])])


# ----------------------------------------------------------------------------- IV statistics
def iv_ratio(Z, Lf, Of, u, idx=None):
    if idx is not None:
        Z, Lf, Of, u = Z[idx], Lf[idx], Of[idx], u[idx]
    sw = u.sum()
    z = Z - np.sum(u * Z) / sw
    num = np.sum(u * z * (Lf - np.sum(u * Lf) / sw))
    den = np.sum(u * z * (Of - np.sum(u * Of) / sw))
    return num / den


def block_boot(fun, blocks, nboot=500):
    ub, inv = np.unique(blocks, return_inverse=True)
    members = [np.where(inv == b)[0] for b in range(ub.size)]
    vals = []
    for _ in range(nboot):
        pick = RNG.integers(0, ub.size, ub.size)
        idx = np.concatenate([members[b] for b in pick])
        vals.append(fun(idx))
    return np.array(vals)


def summarize(est, boots):
    b = boots[np.isfinite(boots)]
    p = np.percentile(b, [2.5, 16, 50, 84, 97.5])
    return dict(ratio=float(est), p2p5=float(p[0]), p16=float(p[1]), p50=float(p[2]), p84=float(p[3]),
                p97p5=float(p[4]), sig_boot=float(0.5 * (p[3] - p[1])), sd_boot=float(np.std(b)))


# ----------------------------------------------------------------------------- Gaussian fits for regspec
def gauss_fit(x, y, s, centres, sig, free_centre=True, width_free=False):
    k = len(centres)

    def model(p):
        m = np.full_like(x, p[0]) + p[1] * (x - x.mean())
        for j in range(k):
            a = p[2 + j]
            c = p[2 + k + j] if free_centre else centres[j]
            sg = p[-1] if width_free else sig
            m += a * np.exp(-0.5 * ((x - c) / sg) ** 2) / (np.sqrt(2 * np.pi) * sg)
        return m

    p0 = [0.0, 0.0] + [0.05] * k + (list(centres) if free_centre else []) + ([sig] if width_free else [])
    lb = [-np.inf, -np.inf] + [-np.inf] * k + ([c - 3.0 for c in centres] if free_centre else []) + \
         ([0.3] if width_free else [])
    ub = [np.inf, np.inf] + [np.inf] * k + ([c + 3.0 for c in centres] if free_centre else []) + \
         ([3.0] if width_free else [])
    r = least_squares(lambda p: (model(p) - y) / s, p0, bounds=(lb, ub))
    return r.x, model(r.x)


def main():
    meta = np.load(D / 'verify_desi_coadd_meta.npz')
    ids, F, I, M, src = load_spectra()
    pos = {s: i for i, s in enumerate(ids)}
    take = np.array([pos[s] for s in meta['sparcl_id']])
    F, I, M, src = F[take], I[take], M[take], src[take]
    R = F - M
    z = meta['redshift']
    ra, dec = meta['ra'], meta['dec']
    vcorr = meta['vcorr']
    ok = meta['ok_contrib'] & np.isfinite(meta['Z']) & np.isfinite(meta['f_chrom_5008'])
    print('coadds', len(ids), 'usable', ok.sum())
    skym, skyrel = sky_mask(I)
    lsf, nlsf = lsf_sigma()

    Ls = build_lines()
    neb_rest = [L[k] for k in L if not (k.startswith('Hg') or k.startswith('Na') or k.startswith('sky'))]
    flux, ferr = {}, {}
    for nm, (lo, hi, lr, kind) in Ls.items():
        others = [vac(x) for x in neb_rest if not (lo - 1 < vac(x) < hi + 1)]
        flux[nm], ferr[nm] = box_flux(R, I, z, lo, hi, others, skym)
    # topocentric Hg I 4358 box, per-fibre barycentric shift
    hgc = L['HgI4358']
    shift = hgc * vcorr / C
    others = [vac(x) for x in neb_rest if abs(vac(x) - hgc) > 2]
    flux['HgI4358'], ferr['HgI4358'] = box_flux(R, I, z, hgc - 2.0, hgc + 2.0, others, skym, shift=shift)
    sk = L['sky5577']
    flux['sky5577'], ferr['sky5577'] = box_flux(R, I, z, sk - HW, sk + HW, [], skym, shift=sk * vcorr / C)
    Ls['HgI4358'] = (hgc - 2, hgc + 2, hgc, 'sky')
    Ls['sky5577'] = (sk - HW, sk + HW, sk, 'sky')

    # capture fraction of each box for the LSF at that wavelength (centred)
    capt = {}
    for nm, (lo, hi, lr, kind) in Ls.items():
        c = 0.5 * (lo + hi)
        sg = np.interp(c, WAVE[np.isfinite(lsf)] if lsf.size == WAVE.size else
                       np.linspace(WAVE[0], WAVE[-1], lsf.size), lsf) if lsf.size else 0.7
        sg = float(np.sqrt(sg ** 2 + (12.0 / C * c) ** 2))
        if nm in ('OII3727sum', 'OII7325sum'):
            capt[nm] = 1.0       # sum boxes extend HW beyond both components
        else:
            capt[nm] = float(erf((0.5 * (hi - lo)) / (np.sqrt(2) * sg)))

    # energy surface brightness proxy: flux * f(lambda) (Omega common) -> chromatic-corrected energy ratios
    E, U = {}, {}
    for nm, (lo, hi, lr, kind) in Ls.items():
        f = accept_interp(meta, lr) if kind != 'sky' else meta['f_chrom_5008']
        E[nm] = flux[nm] * f / capt[nm]
        U[nm] = 1.0 / (ferr[nm] * f / capt[nm]) ** 2
    # robust per-line quality cut
    good = {}
    for nm in E:
        e = 1 / np.sqrt(U[nm])
        m = ok & np.isfinite(E[nm]) & np.isfinite(e) & (e > 0)
        m &= e < 4 * np.nanmedian(e[m])
        m &= np.abs(E[nm] / e) < 30
        good[nm] = m
    O = E['OIII5007']
    Z = meta['Z']
    hp16 = healpix_ring(ra, dec, 16)
    hp32 = healpix_ring(ra, dec, 32)

    # matched filter detection on box fluxes (sanity)
    m0 = good['OIII5007']
    uo = U['OIII5007'][m0]
    A = np.sum(uo * Z[m0] * O[m0]) / np.sum(uo * Z[m0] ** 2)
    snr = A * np.sqrt(np.sum(uo * Z[m0] ** 2))
    chi_null = {nm: float(np.std(E[nm][good[nm]] * np.sqrt(U[nm][good[nm]])))
                for nm in E if nm.startswith('null')}
    res = dict(n=int(ok.sum()), box_halfwidth_A=HW, capture=capt, detection=dict(A_flux_units=float(A),
               snr_formal=float(snr), chi_std_O=float(np.std(O[m0] * np.sqrt(uo)))), null_chi_std=chi_null)
    print('box [O III] matched-filter S/N', round(snr, 1), 'chi std O', res['detection']['chi_std_O'])

    subsets = {
        'all': ok,
        'south_dec<19.5': ok & (dec < 19.5),
        'north_dec>=19.5': ok & (dec >= 19.5),
        'Zpos': ok & (Z > 0),
        'Zneg': ok & (Z < 0),
        'south_Zpos': ok & (dec < 19.5) & (Z > 0),
        'south_Zneg': ok & (dec < 19.5) & (Z < 0),
        'QSO_or_z>0.5': ok & ((meta['spectype'] == 'QSO') | (z > 0.5)),
        'w1_fibres': ok & (src == 'w1'),
        'new_fibres': ok & (src == 'ext'),
    }
    lines_out = [nm for nm in E if nm != 'OIII5007']
    ratios = {}
    for sname, smask in subsets.items():
        ratios[sname] = {}
        for nm in lines_out:
            if sname != 'all' and (nm.startswith('null') or nm in ('Hd', 'NII6548', 'HeI5876')):
                continue
            m = smask & good[nm] & good['OIII5007']
            idx = np.where(m)[0]
            if idx.size < 50:
                continue
            u = U[nm][idx]
            fun = (lambda ii, idx=idx, u=u, nm=nm:
                   iv_ratio(Z[idx][ii], E[nm][idx][ii], O[idx][ii], u[ii]))
            est = fun(np.arange(idx.size))
            b16 = block_boot(fun, hp16[idx], 400)
            b32 = block_boot(fun, hp32[idx], 400)
            r = summarize(est, b32)
            r['nside16'] = summarize(est, b16)
            r['n'] = int(idx.size)
            ratios[sname][nm] = r
        print(sname, int(smask.sum()), {k: round(v['ratio'], 3) for k, v in ratios[sname].items()
                                         if k in ('Hb', 'Ha', 'NII6583', 'OII3727sum', 'OIII4363', 'HeII4686')})
    res['ratios'] = ratios
    json.dump(res, open(D / 'verify_desi_audit.json', 'w'), indent=1)
    # empirical scatter of null-box ratios vs their bootstrap sigma
    nl = [k for k in ratios['all'] if k.startswith('null')]
    res['null_summary'] = dict(
        names=nl, ratios=[ratios['all'][k]['ratio'] for k in nl],
        sig_boot32=[ratios['all'][k]['sig_boot'] for k in nl], sig_boot16=[ratios['all'][k]['nside16']['sig_boot']
                                                                         for k in nl],
        rms_ratio=float(np.sqrt(np.mean([ratios['all'][k]['ratio'] ** 2 for k in nl]))),
        rms_over_sig32=float(np.sqrt(np.mean([(ratios['all'][k]['ratio'] / ratios['all'][k]['sig_boot']) ** 2
                                               for k in nl]))),
        rms_over_sig16=float(np.sqrt(np.mean([(ratios['all'][k]['ratio'] / ratios['all'][k]['nside16']['sig_boot'])
                                               ** 2 for k in nl]))))
    print('null rms ratio', res['null_summary']['rms_ratio'], 'rms/sig32', res['null_summary']['rms_over_sig32'],
          'rms/sig16', res['null_summary']['rms_over_sig16'])

    # ------------------------------------------------------------------------- pixel regression spectra
    def regspec(lam_lo, lam_hi, instr, sel, nboot=300, mask_lines=()):
        """r_c(lambda) = a + sum_k b_k(lambda) Z_k ; returns b_k(lambda) per unit of the [O III] box slope."""
        pix = np.where((WAVE >= lam_lo) & (WAVE <= lam_hi))[0]
        idx = np.where(sel & good['OIII5007'])[0]
        # continuum per fibre: weighted linear fit over the window excluding +-3 A of lines and sky mask
        Rc = np.full((idx.size, pix.size), np.nan)
        Wt = np.zeros((idx.size, pix.size))
        x = WAVE[pix] - WAVE[pix].mean()
        for jj, i in enumerate(idx):
            m = (I[i, pix] > 0) & ~skym[pix]
            for lw in mask_lines:
                m &= np.abs(WAVE[pix] - lw) > 3.0
            for lw in GAL_REST * (1 + z[i]):
                m &= np.abs(WAVE[pix] - lw) > 8.0
            if m.sum() < 20 or np.any(I[i, pix] <= 0):
                continue
            Am = np.c_[np.ones(m.sum()), x[m]]
            p = np.linalg.lstsq(Am * np.sqrt(I[i, pix][m])[:, None], R[i, pix][m] * np.sqrt(I[i, pix][m]),
                                rcond=None)[0]
            zl = GAL_REST * (1 + z[i])
            bad = np.zeros(pix.size, bool)
            for lw in zl:
                bad |= np.abs(WAVE[pix] - lw) < 8.0
            Rc[jj] = np.where(bad, np.nan, (R[i, pix] - (p[0] + p[1] * x)) * meta['f_chrom_5008'][i])
            Wt[jj] = np.where(bad, 0.0, I[i, pix] / meta['f_chrom_5008'][i] ** 2)
        okr = np.isfinite(Rc).any(1)
        idx, Rc, Wt = idx[okr], Rc[okr], Wt[okr]
        Rc = np.nan_to_num(Rc)
        cols = [ins[idx] for ins in instr]
        keepc = [j for j, cc in enumerate(cols) if np.std(cc) > 1e-4]
        X = np.c_[np.ones(idx.size), *[cols[j] for j in keepc]]
        regspec.kept = keepc
        # per-fibre weight: median pixel weight (keeps the design common across pixels)
        wf = np.median(np.where(Wt > 0, Wt, np.nan), 1)
        wf = np.nan_to_num(wf)

        def solve(ii):
            Xi = X[ii]
            w = wf[ii]
            XtW = Xi.T * w
            B = np.linalg.solve(XtW @ Xi, XtW @ Rc[ii])
            # [O III] box slope with the same design
            Bo = np.linalg.solve(XtW @ Xi, XtW @ O[idx][ii])
            return B, Bo

        B, Bo = solve(np.arange(idx.size))
        ub, inv = np.unique(hp32[idx], return_inverse=True)
        members = [np.where(inv == b)[0] for b in range(ub.size)]
        boots = []
        for _ in range(nboot):
            pick = RNG.integers(0, ub.size, ub.size)
            ii = np.concatenate([members[b] for b in pick])
            boots.append(solve(ii))
        return WAVE[pix], B, Bo, boots, idx.size

    # 3. [O III]4363 vs Hg I 4358 (instrument: Z; spectrum per unit of [O III] slope)
    lam, B, Bo, boots, n43 = regspec(4335.0, 4395.0, [Z], ok,
                                     mask_lines=[vac(L['Hg']), vac(L['OIII4363']), L['HgI4358'], vac(L['FeII4359'])])
    spec = B[1] / Bo[1]
    bspec = np.array([b[0][1] / b[1][1] for b in boots])
    s_sp = np.std(bspec, 0)
    sg = float(np.interp(4364.6, WAVE, lsf)) if lsf.size == WAVE.size else 0.72
    win = (lam > 4345) & (lam < 4385) & (np.abs(lam - vac(L['Hg'])) > 3.0)
    p, mod = gauss_fit(lam[win], spec[win], np.maximum(s_sp[win], 1e-6), [vac(L['OIII4363']), L['HgI4358']], sg)
    cb, ab = [], []
    for bs in bspec:
        pb, _ = gauss_fit(lam[win], bs[win], np.maximum(s_sp[win], 1e-6), [vac(L['OIII4363']), L['HgI4358']], sg)
        cb.append(pb[4:6])
        ab.append(pb[2:4])
    cb, ab = np.array(cb), np.array(ab)
    # fixed-centre fit (amplitudes only)
    pf, _ = gauss_fit(lam[win], spec[win], np.maximum(s_sp[win], 1e-6), [vac(L['OIII4363']), L['HgI4358']], sg,
                      free_centre=False)
    abf = []
    for bs in bspec:
        pb, _ = gauss_fit(lam[win], bs[win], np.maximum(s_sp[win], 1e-6), [vac(L['OIII4363']), L['HgI4358']], sg,
                          free_centre=False)
        abf.append(pb[2:4])
    abf = np.array(abf)
    acc4363 = float(np.nanmedian(meta['f_chrom_4363'] / meta['f_chrom_5008']))
    res['regspec_4363'] = dict(
        n=int(n43), lsf_sigma_A=sg, expected_OIII4363=float(vac(L['OIII4363'])), expected_Hg=float(L['HgI4358']),
        centre_4363=float(p[4]), centre_4363_p16_p84=np.percentile(cb[:, 0], [16, 84]).tolist(),
        centre_Hg=float(p[5]), centre_Hg_p16_p84=np.percentile(cb[:, 1], [16, 84]).tolist(),
        amp_4363_free=float(p[2]), amp_Hg_free=float(p[3]),
        amp_4363_fixed=float(pf[2]), amp_4363_fixed_p16_p84=np.percentile(abf[:, 0], [2.5, 16, 50, 84, 97.5]).tolist(),
        amp_Hg_fixed=float(pf[3]), amp_Hg_fixed_p16_p84=np.percentile(abf[:, 1], [2.5, 16, 50, 84, 97.5]).tolist(),
        note='amplitudes = line flux per unit [O III]5007 flux (photon-count like, no chromatic correction; '
             f'f(4363)/f(5008) median {acc4363:.3f})')
    print('4363 regspec', {k: (round(v, 3) if isinstance(v, float) else v) for k, v in res['regspec_4363'].items()
                           if k != 'note'})
    np.savez_compressed(D / 'verify_desi_audit_regspec.npz', lam4363=lam, spec4363=spec, sig4363=s_sp,
                        model4363=np.interp(lam, lam[win], mod))

    # 4. velocity, E vs W strand, joint pixel regression around 5008
    vel = {}
    lam0 = L['OIII5007']
    for sname, sel in (('all', ok), ('dec15.5-17.5', ok & (dec > 15.5) & (dec < 17.5)),
                       ('dec<15.5', ok & (dec < 15.5)), ('dec>=17.5', ok & (dec >= 17.5))):
        lamv, B, Bo, boots, nv = regspec(4985.0, 5031.0, [meta['Z_E'], meta['Z_W'], meta['Z_W2']], sel,
                                         mask_lines=[vac(lam0, 13.0)])
        sgv = float(np.interp(5008.2, WAVE, lsf)) if lsf.size == WAVE.size else 0.70
        sgv = float(np.sqrt(sgv ** 2 + (12.0 / C * 5008.2) ** 2))
        wv = (lamv > 4995) & (lamv < 5021)
        out = {}
        names_k = [('E', 'W', 'W2')[j] for j in regspec.kept]
        for k, nmk in enumerate(names_k, start=1):
            y = B[k]
            sb = np.std(np.array([b[0][k] for b in boots]), 0)
            try:
                pk, _ = gauss_fit(lamv[wv], y[wv], np.maximum(sb[wv], 1e-9), [vac(lam0, 13.0)], sgv)
            except Exception:
                continue
            vb = []
            for b in boots:
                try:
                    pb, _ = gauss_fit(lamv[wv], b[0][k][wv], np.maximum(sb[wv], 1e-9), [vac(lam0, 13.0)], sgv)
                    if abs(pb[2]) > 0:
                        vb.append((pb[3] / lam0 - 1) * C)
                except Exception:
                    pass
            vb = np.array(vb)
            amp_snr = float(pk[2] / np.std([b[0][k][wv].sum() * DLAM for b in boots])) if len(boots) else np.nan
            out[nmk] = dict(v=float((pk[3] / lam0 - 1) * C), v_p16_p84=np.percentile(vb, [16, 84]).tolist(),
                            v_sd=float(np.std(vb)), amp=float(pk[2]), amp_snr=amp_snr)
        lamt, Bt, Bot, bootst, _ = regspec(4985.0, 5031.0, [Z], sel, mask_lines=[vac(lam0, 13.0)])
        sbt = np.std(np.array([b[0][1] for b in bootst]), 0)
        pt, _ = gauss_fit(lamt[wv], Bt[1][wv], np.maximum(sbt[wv], 1e-9), [vac(lam0, 13.0)], sgv)
        vbt = []
        for b in bootst:
            pb, _ = gauss_fit(lamt[wv], b[0][1][wv], np.maximum(sbt[wv], 1e-9), [vac(lam0, 13.0)], sgv)
            vbt.append((pb[3] / lam0 - 1) * C)
        out['total_Z'] = dict(v=float((pt[3] / lam0 - 1) * C), v_p16_p84=np.percentile(vbt, [16, 84]).tolist(),
                              v_sd=float(np.std(vbt)), amp=float(pt[2]), amp_snr=np.nan)
        cE = np.corrcoef(meta['Z_E'][sel], meta['Z_W'][sel])[0, 1]
        out['corr_ZE_ZW'] = float(cE)
        if 'E' in out and 'W' in out:
            out['E_minus_W'] = out['E']['v'] - out['W']['v']
        out['n'] = int(nv)
        vel[sname] = out
        print('velocity', sname, out.get('corr_ZE_ZW'), {k: (round(v['v'], 1), round(v['v_sd'], 1), round(v['amp'], 4))
                                  for k, v in out.items() if isinstance(v, dict)})
    res['velocity_joint_regression'] = vel

    # 5. calibration bookkeeping straight from the GFA exposure table
    rows = list(csv.DictReader(open(ROOT / 'data/raw/verify_desi_tap/exposure_gfa_box.csv')))
    ff = np.array([float(r['fiber_fracflux_gfa']) if r['fiber_fracflux_gfa'] not in ('', 'nan') else np.nan
                   for r in rows])
    prog = np.array([r['program'] for r in rows])
    g = np.isfinite(ff) & (ff > 0.05)
    cal = dict(n_exp=len(rows), n_with_ff=int(g.sum()), median_ff=float(np.median(ff[g])),
               by_program={p: float(np.median(ff[g & (prog == p)])) for p in set(prog) if (g & (prog == p)).any()},
               legacy_1arcsec_gaussian_1p5arcsec=float(1 - np.exp(-0.75 ** 2 / (2 * (1 / 2.3548) ** 2))),
               coadd_f5008_median=float(np.nanmedian(meta['f_chrom_5008'])),
               ratio_f5008_over_psf2fib=np.nanpercentile(meta['f_chrom_5008'] /
                                                         np.clip(meta['mean_psf_to_fiber_specflux'], 0.3, 1),
                                                         [16, 50, 84]).tolist())
    res['calibration'] = cal
    print('calibration', cal)

    # 6. tiles 41485/41487 among strong-[O III] fibres
    sn = O / np.sqrt(1 / U['OIII5007'])
    strong = good['OIII5007'] & (sn > 4)
    res['tiles_41485_41487'] = dict(n_strong=int(strong.sum()),
                                    n_strong_on_bad_tiles=int((strong & meta['has_41485_41487']).sum()),
                                    n_any_on_bad_tiles=int((ok & meta['has_41485_41487']).sum()))
    print(res['tiles_41485_41487'])
    json.dump(res, open(D / 'verify_desi_audit.json', 'w'), indent=1)
    np.savez_compressed(D / 'verify_desi_audit_boxfluxes.npz', sparcl_id=meta['sparcl_id'],
                        **{f'E_{k}': v for k, v in E.items()}, **{f'U_{k}': v for k, v in U.items()},
                        skyrel=skyrel)


if __name__ == '__main__':
    main()
