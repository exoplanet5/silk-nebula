"""Build the target list for the LAMOST MRS-N plate MN074209N160100N01 (2024-02-01, the paper's Silk plate):
its 83 public star spectra (80 F-Std + 3 Obj; DR12 v2.0 med_catalogue) become 'on' targets of
spectra_sdss_lamost_lamost_resid.py (mode mrsmn); stellar templates come from the 900 OFF-strand stars of the
other MRS plates (data/raw/lamost_dr12/mrs_strand_sample.csv).  On this plate each spectrograph's super-sky
is built from its ~165 nebular fibres (header NSKIES), so the residual at 5007 A in a star is
I([O III], star position) - <I([O III])> over that spectrograph's nebular fibres.  Topic tag spectra_sdss_lamost.

Output: data/derived/spectra_sdss_lamost_mnplate_sample.csv
"""
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import spectra_sdss_lamost_trace as TR

ROOT = Path('/Users/mickey/astronomy/silk-nebula')


def main():
    t = Table.read(ROOT / 'data/raw/lamost_dr12/med_catalogue_MN074209N160100N01.csv', format='ascii.csv')
    rows = []
    for o in sorted(set(t['obsid'])):
        r = t[t['obsid'] == o]
        rv = np.array(r['rv_br1'], float); rv = rv[rv > -999]
        rl = np.array(r['rv_lasp1'], float); rl = rl[rl > -999]
        rows.append(dict(obsid=int(o), planid=str(r['planid'][0]), obsdate=str(r['obsdate'][0]), spid=int(r['spid'][0]),
                         fiberid=int(r['fiberid'][0]), ra=float(r['ra'][0]), dec=float(r['dec'][0]),
                         snr=float(np.max(r['snr'])), fibertype=str(r['fibertype'][0]),
                         gaia_g_mean_mag=float(r['gaia_g_mean_mag'][0]),
                         rv_br1=float(np.median(rv)) if len(rv) else -9999.0,
                         rv_br1_err=-9999.0, rv_lasp1=float(np.median(rl)) if len(rl) else -9999.0))
    M = Table(rows=rows)
    geo = TR.load()
    d, ins, _ = TR.nearest_strand_distance(np.array(M['ra']), np.array(M['dec']), geo)
    M['dist'] = d
    M['sample'] = 'on'
    S = Table.read(ROOT / 'data/raw/lamost_dr12/mrs_strand_sample.csv', format='ascii.csv')
    S = S[S['sample'] == 'off']
    keep = [c for c in M.colnames if c in S.colnames]
    out = vstack([M[keep], S[keep]])
    out.write(ROOT / 'data/derived/spectra_sdss_lamost_mnplate_sample.csv', overwrite=True)
    print(len(M), 'MN-plate stars;', len(S), 'template stars; MN |dist| range %.2f-%.2f deg' % (np.min(np.abs(d)), np.max(np.abs(d))))
    print('MN stars with valid RV:', int(np.sum((M['rv_br1'] > -999) | (M['rv_lasp1'] > -999))))


if __name__ == '__main__':
    main()
