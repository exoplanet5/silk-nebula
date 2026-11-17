"""spectra_desi_lsf.py -- median DESI LSF sigma (A) in the saved windows from SPARCL 'wave_sigma' (pilot).
Output: data/derived/spectra_desi_lsf.json {wave, sigma_A, sigma_p16, sigma_p84, per_line{...}}"""
import json
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_desi_lines import WIN_WAVE, LINES, C_KMS  # noqa: E402

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
S = []
for f in sorted((ROOT / 'data/derived/spectra_desi_windows').glob('pilot_*.npz')):
    with np.load(f, allow_pickle=False) as d:
        if 'wave_sigma' in d.files:
            S.append(d['wave_sigma'])
S = np.concatenate(S)
S = np.where(S > 0, S, np.nan)
med = np.nanmedian(S, 0)
p16, p84 = np.nanpercentile(S, [16, 84], 0)
per = {}
for k, l0 in LINES.items():
    i = np.argmin(np.abs(WIN_WAVE - l0))
    per[k] = dict(lam=round(l0, 2), sigma_A=round(float(med[i]), 4), p16=round(float(p16[i]), 4),
                  p84=round(float(p84[i]), 4), fwhm_kms=round(float(2.3548 * med[i] / l0 * C_KMS), 1))
out = dict(n_spectra=int(S.shape[0]), wave=WIN_WAVE.round(2).tolist(), sigma_A=np.round(med, 4).tolist(),
           sigma_p16=np.round(p16, 4).tolist(), sigma_p84=np.round(p84, 4).tolist(), per_line=per,
           source='SPARCL DESI-DR1 wave_sigma (per-pixel LSF Gaussian sigma), pilot spectra')
(ROOT / 'data/derived/spectra_desi_lsf.json').write_text(json.dumps(out))
for k, v in per.items():
    print(k, v)
