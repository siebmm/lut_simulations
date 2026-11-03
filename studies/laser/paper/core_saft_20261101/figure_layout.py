"""Reserve legend space and check that keys cannot hide scientific data."""
from pathlib import Path
import json
from matplotlib.legend import Legend
from matplotlib.text import Text

HERE = Path(__file__).resolve().parent


def legend_strip(fig, spec, handles, labels=None, **kwargs):
    """A layout-managed, data-free row: legends never float over a plot."""
    ax = fig.add_subplot(spec)
    ax.set_gid('legend-strip')
    ax.set_axis_off()
    if labels is None:
        labels = [h.get_label() for h in handles]
    return ax.legend(handles, labels, loc='center', frameon=False,
                     borderaxespad=0, **kwargs)


def verify_layout(fig, name):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    legends = fig.findobj(Legend)
    data_axes = [ax for ax in fig.axes if ax.get_gid() != 'legend-strip']
    records = []
    for legend in legends:
        box = legend.get_window_extent(renderer)
        assert fig.bbox.contains(box.x0, box.y0) and fig.bbox.contains(box.x1, box.y1), (name, 'clipped legend')
        collisions = []
        for ax in data_axes:
            if box.overlaps(ax.get_window_extent(renderer)):
                collisions.append(ax.get_title() or ax.get_ylabel() or 'data axes')
            for item in ax.findobj(Text):
                if item.get_visible() and item.get_text() and box.overlaps(item.get_window_extent(renderer)):
                    collisions.append(item.get_text())
        assert not collisions, (name, 'legend covers plot or label', collisions)
        records.append({'labels': [t.get_text() for t in legend.get_texts()],
                        'bounds_figure_fraction': list(box.transformed(fig.transFigure.inverted()).bounds),
                        'data_or_label_collisions': collisions})
    folder = HERE / 'verification/figure_layout'
    folder.mkdir(parents=True, exist_ok=True)
    (folder / (name + '.json')).write_text(json.dumps({
        'status': 'passed', 'figure': name, 'legend_count': len(legends),
        'legends': records,
        'check': 'Rendered legend bounds inside canvas and outside all data axes and their text.'
    }, indent=2) + '\n', encoding='utf-8')
