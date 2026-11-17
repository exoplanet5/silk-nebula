"""spectra_desi_petalmodel.py -- forward model of DESI per-petal sky subtraction for a structured nebula.

DESI subtracts, per exposure and per spectrograph (= petal), a sky spectrum fitted to that petal's sky fibres
(desispec.sky.compute_sky; plus low-order skygradpca/tpcorr terms with fixed spectral templates). A Galactic
emission line present in the sky fibres is therefore removed at its petal-mean level, and a fibre's coadded
residual is   res_i = T(x_i) - sum_t w_t <T>_{tile t, petal p(i,t)}   with w_t the coadd weights.
Sky-fibre positions are not in the public desi_dr1.fiberassign table, so <T> is approximated by the mean of
T over all science/standard fibres of that tile+petal (they sample the petal area; sky fibres are spread
across it too). Weights w_t ~ summed EFFTIME_SPEC of the tile's exposures (coadds are ivar-weighted).

Provides
  PetalModel(): loads fiberassign_box.csv + exposure_box.csv
  .design(ra, dec, targetid, survey, program, basis) -> matrix [N, K] of predicted residuals per unit basis amp
Basis functions: Gaussian ridges along the trace strands in Dec segments (see ridge_basis()).
"""
import csv
from collections import defaultdict
import numpy as np

import sys
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_geometry import ROOT, load_trace  # noqa: E402

FA = ROOT / 'data/raw/datalab_tap/fiberassign_box.csv'
EXP = ROOT / 'data/raw/datalab_tap/exposure_box.csv'
TILES = ROOT / 'data/raw/datalab_tap/tiles_box.csv'


def tangent(ra, dec, ra0, dec0):
    """Gnomonic offsets (deg) about (ra0, dec0)."""
    ra, dec, ra0, dec0 = map(np.radians, (ra, dec, ra0, dec0))
    cosc = np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0)
    x = np.cos(dec) * np.sin(ra - ra0) / cosc
    y = (np.cos(dec0) * np.sin(dec) - np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)) / cosc
    return np.degrees(x), np.degrees(y)


class RidgeBasis:
    """Gaussian cross-section ridges along polylines, split into Dec segments with triangular (linear-interp)
    amplitude weights so the fitted amplitude varies smoothly along each strand."""

    def __init__(self, strands=('E', 'W', 'W2'), seg_deg=1.0, hwhm=None, ra_shift=0.0):
        t = load_trace()
        self.lines = {}
        self.ra_shift = ra_shift          # deg of RA added to every strand (null tests with a displaced trace)
        sp = t['spine']
        self.lines['E'] = (np.array([p['ra'] for p in sp]), np.array([p['dec'] for p in sp]),
                           np.array([p['hwhm_deg'] for p in sp]))
        for k in ('W', 'W2'):
            if k in strands and k in t.get('strands', {}):
                pts = t['strands'][k]['points']
                self.lines[k] = (np.array([p['ra'] for p in pts]), np.array([p['dec'] for p in pts]),
                                 np.array([p.get('hwhm_deg', t['half_width_deg']) for p in pts]))
        self.hwhm_override = hwhm
        self.names, self.nodes = [], []
        for k, (ra, dec, hw) in self.lines.items():
            if k not in strands:
                continue
            nodes = np.arange(np.ceil(dec.min() / seg_deg) * seg_deg, dec.max() + 1e-6, seg_deg)
            nodes = np.unique(np.concatenate([[dec.min()], nodes, [dec.max()]]))
            for j, d in enumerate(nodes):
                self.names.append(f'{k}@{d:.2f}')
                self.nodes.append((k, j, nodes))

    def evaluate(self, ra, dec):
        """[N, K] basis values at positions (unit peak amplitude)."""
        ra, dec = np.asarray(ra, float), np.asarray(dec, float)
        out = np.zeros((ra.size, len(self.names)))
        for kk, (k, j, nodes) in enumerate(self.nodes):
            lra, ldec, lhw = self.lines[k]
            lra = lra + self.ra_shift
            o = np.argsort(ldec)
            # ridge RA and width at each Dec (within strand Dec range)
            inside = (dec >= ldec.min() - 0.3) & (dec <= ldec.max() + 0.3)
            rra = np.interp(dec, ldec[o], lra[o])
            hw = np.interp(dec, ldec[o], lhw[o]) if self.hwhm_override is None else self.hwhm_override
            hw = np.clip(hw, 0.05, 0.4)
            dx = (ra - rra) * np.cos(np.radians(dec))
            prof = np.exp(-0.5 * (dx / (hw / 1.1774)) ** 2)
            # taper beyond strand ends (0.3 deg Gaussian roll-off)
            dd = np.where(dec < ldec.min(), ldec.min() - dec, np.where(dec > ldec.max(), dec - ldec.max(), 0))
            prof *= np.exp(-0.5 * (dd / 0.15) ** 2) * inside
            # triangular weight in Dec around node j
            wts = np.zeros_like(dec)
            d0 = nodes[j]
            if j > 0:
                lo = nodes[j - 1]
                m = (dec >= lo) & (dec <= d0)
                wts[m] = (dec[m] - lo) / (d0 - lo)
            else:
                wts[dec <= d0] = 1.0
            if j < len(nodes) - 1:
                hi = nodes[j + 1]
                m = (dec > d0) & (dec <= hi)
                wts[m] = (hi - dec[m]) / (hi - d0)
            else:
                wts[dec > d0] = 1.0
            out[:, kk] = prof * wts
        return out


class PetalModel:
    def __init__(self):
        eff = defaultdict(float)
        for r in csv.DictReader(open(EXP)):
            try:
                eff[int(r['tileid'])] += float(r['efftime_spec'])
            except ValueError:
                pass
        self.tile_eff = dict(eff)
        self.tile_prog = {int(r['tileid']): (r['survey'], r['program']) for r in csv.DictReader(open(TILES))}
        self.by_target = defaultdict(list)       # targetid -> [(tileid, petal)]
        pet_ra, pet_dec = defaultdict(list), defaultdict(list)
        for r in csv.DictReader(open(FA)):
            t, p = int(r['tileid']), int(r['petal_loc'])
            self.by_target[int(r['targetid'])].append((t, p))
            pet_ra[(t, p)].append(float(r['target_ra']))
            pet_dec[(t, p)].append(float(r['target_dec']))
        self.pet_pos = {k: (np.array(pet_ra[k]), np.array(pet_dec[k])) for k in pet_ra}

    def petal_means(self, basis):
        self.pmean = {k: basis.evaluate(*v).mean(0) for k, v in self.pet_pos.items()}
        return self.pmean

    def contributions(self, targetid, survey, program):
        """List of (tileid, petal, weight) for the coadd of this target in (survey, program)."""
        out = []
        for t, p in self.by_target.get(int(targetid), []):
            if self.tile_prog.get(t) == (survey, program) and self.tile_eff.get(t, 0) > 0:
                out.append((t, p, self.tile_eff[t]))
        s = sum(w for _, _, w in out)
        return [(t, p, w / s) for t, p, w in out] if s > 0 else []

    def design(self, ra, dec, targetid, survey, program, basis):
        if not hasattr(self, 'pmean'):
            self.petal_means(basis)
        X = basis.evaluate(ra, dec)
        sub = np.zeros_like(X)
        ok = np.ones(len(ra), bool)
        info = []
        for i in range(len(ra)):
            c = self.contributions(targetid[i], survey[i], program[i])
            if not c:
                ok[i] = False
                info.append([])
                continue
            for t, p, w in c:
                sub[i] += w * self.pmean[(t, p)]
            info.append(c)
        return X - sub, X, sub, ok, info
