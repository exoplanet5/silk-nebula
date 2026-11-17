"""Quick-look mosaic of every EQ grid with the trace overlaid (diagnostic only; scratch output)."""
import sys, json
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import new_sky_figure, finish
from imaging_mw_common import DER, load_grid, overlay_eq, eq_axes, gal_axes, EQ_EXTENT, GAL_EXTENT, stretch, get_spine, spine_gal
which = sys.argv[1] if len(sys.argv) > 1 else 'eq'
out = sys.argv[2]
cat = json.loads((DER / 'imaging_mw_grids' / 'grid_catalog.json').read_text())
names = sorted(k.split(':')[1] for k in cat if k.startswith(which + ':') and not k.endswith('hdu1'))
n = len(names); ncol = 8; nrow = int(np.ceil(n / ncol))
fig, axes = new_sky_figure(2.6 * ncol, 3.3 * nrow, 'quicklook', which, nrows=nrow, ncols=ncol)
sp = get_spine()
for ax, nm in zip(axes.ravel(), names):
    a = load_grid(nm, which)
    lo, hi = stretch(a, 2, 98)
    ext = EQ_EXTENT if which == "eq" else GAL_EXTENT
    ax.imshow(a, origin='lower', extent=ext, vmin=lo, vmax=hi, cmap='gray', interpolation='nearest')
    if which == 'eq':
        overlay_eq(ax, sp, label_regions=False, color='#ff5252', lw=0.6, lamost=False); eq_axes(ax)
    else:
        l, b = spine_gal(sp); ax.plot(l, b, color='#ff5252', lw=0.6); gal_axes(ax)
    ax.set_title(nm, fontsize=7); ax.tick_params(labelsize=5); ax.set_xlabel(''); ax.set_ylabel('')
for ax in axes.ravel()[n:]:
    ax.axis('off')
finish(fig, out, dpi=80)
