r"""
+---------------------+-----------+
|    ^                | OQMD refs |
|    |  *       *     | ...       |
|    |   \   * /      |           |
| ΔH,|    \   /       | C2DB refs |
| eV/|     \ /        | ...       |
| atm|      *         |           |
|    |                |           |
|     ------------>   |           |
|         A   B       |           |
|          1-x x      |           |
+---------------------+-----------+
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from typing import Iterable

import plotly
import plotly.graph_objs as go
from ase.formula import Formula
from ase.phasediagram import PhaseDiagram

from camdweb.c2db.asr_panel import read_result_file
from camdweb.html import table
from camdweb.material import Material
from camdweb.panels.panel import Panel, PanelData
from camdweb.utils import html_format_formula

HTML = """
<div class="row">
  <div class="col-6">
    {table}
    <div id='{id}' class='{id}'></div>
  </div>
  <div class="col-6">
    {tables}
  </div>
</div>
"""

SCRIPT = """
<script type='text/javascript'>
var graphs = {chull_json};
Plotly.newPlot('{id}', graphs, {{}});
</script>
"""

# D3 colors from plotly.express.colors.qualitative
colors = ['#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD',
          '#8C564B', '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF']


class ConvexHullPanel(Panel):
    def __init__(self, sources: dict[str, tuple[str, str]] | None = None):
        super().__init__()
        self.sources = sources

    def get_data(self,
                 material: Material) -> PanelData:
        tbl = table(
            [f'{html_format_formula(material.formula)} ({material.uid})', ''],
            [['Heat of formation [eV/atom]', f'{material.hform:.2f}'],
             ['Energy above convex hull [eV/atom]', f'{material.ehull:.2f}']])
        root = material.folder.parents[2] / 'convex-hulls'
        name = ''.join(sorted(material.count))
        ch_file = root / f'{name}.json'
        refs = read_result_file(ch_file)
        chull, tables = make_figure_and_tables(refs,
                                               higlight_uid=material.uid,
                                               sources=self.sources,
                                               verbose=False)
        html = HTML.format(table=tbl, tables=tables, id='chull')
        if chull:
            script = SCRIPT.format(chull_json=chull, id='chull')
        else:
            script = ''
        return PanelData(html,
                         title='Convex hull',
                         script=script)


def make_figure_and_tables(refs: dict[str, tuple[dict[str, int],
                                                 float,
                                                 str]],
                           higlight_uid: str | None = None,
                           sources: dict[str, tuple[str, str]] | None = None,
                           verbose: bool = True) -> tuple[str, str]:
    """Make convex-hull figure and tables.

    >>> refs = {'u1': ({'B': 1}, 0.0, 'OQMD'),
    ...         'u2': ({'B': 1, 'N': 1}, -0.5, 'OQMD'),
    ...         'u3': ({'N': 1}, 0.0, 'OQMD'),
    ...         '11': ({'B': 1, 'N': 1}, -0.2, 'C2DB')}
    >>> ch, tables = make_figure_and_tables(refs)
    Species: B, N
    References: 4
    0    B              0.000
    1    BN            -0.500
    2    N              0.000
    3    BN            -0.200
    Simplices: 2
    """
    if sources is None:
        sources = {}
    tables = defaultdict(list)
    source_names = []
    uids = []

    for uid, (count, e, source) in refs.items():
        f = Formula.from_dict(count)
        hform = e / len(f)
        if uid == higlight_uid or source not in sources:
            frmt = '{formula:html}, ({uid})'
        else:
            _, frmt = sources[source]
        tables[source].append((hform, frmt.format(uid=uid, formula=f)))
        source_names.append(source)
        uids.append(uid)

    pd = PhaseDiagram([(count, e) for count, e, source in refs.values()],
                      verbose=verbose)
    if 2 <= len(pd.symbols) <= 3:
        if len(pd.symbols) == 2:
            fig = plot_2d(pd, uids, source_names, uid=higlight_uid)
        else:
            fig = plot_3d(pd, uids, source_names, uid=higlight_uid)
        chull = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    else:
        chull = ''

    html = '\n'.join(
        table([sources.get(source, ('References', ''))[0], ''],
              [[link, f'{h:.2f} eV/atom'] for h, link in sorted(tbl)])
        for source, tbl in tables.items())
    return chull, html


def plot_2d(pd: PhaseDiagram,
            uids: list[str] | None = None,
            sources: list[str] | None = None,
            uid: str | None = None) -> go.Figure:

    if uids is None:
        uids = [r[2] for r in pd.references]

    if sources is None:
        sources = ['Materials'] * len(uids)

    x, y = pd.points[:, 1:].T

    X = []
    Y = []
    for i, j in pd.simplices:
        X += [x[i], x[j], None]
        Y += [y[i], y[j], None]

    data = [go.Scatter(x=X, y=Y, mode='lines',
                       line=dict(color=colors[2]),
                       showlegend=False)]

    names = [format(Formula(ref[2]).reduce()[0], 'html')
             for ref in pd.references]

    for i, source in enumerate(set(sources)):
        mask = [True if source in label else False
                for label in sources]
        hovertemplate = '%{customdata} <br> ΔH: %{y:.2f} eV/atom'

        data.append(go.Scatter(
            x=x[mask],
            y=y[mask],
            customdata=[f'{name} ({uid})' for uid, name, x
                        in zip(uids, names, mask) if x],
            name=source,
            hovertemplate=hovertemplate,
            mode='markers',
            marker=dict(color=colors[i], size=8)))

    delta = y.ptp() / 30
    ymin = y.min() - 2.5 * delta
    fig = go.Figure(data=data, layout_yaxis_range=[ymin, 0.1])

    # Highlight selected material:
    if uid is not None:
        this_idx = uids.index(uid)
        names[this_idx] = '<b>' + names[this_idx] + '</b>'

    # Add annotations for materials on the hull and the selected material:
    annotate_idx = [i for i, x in enumerate(pd.hull) if x]
    if uid is not None:
        if this_idx not in annotate_idx:
            annotate_idx.append(this_idx)  # pragma: no cover
    for i in annotate_idx:
        fig.add_annotation(x=x[i], y=y[i],
                           text=names[i],
                           xanchor='left',
                           showarrow=False,
                           xshift=10)

    A, B = pd.symbols
    fig.update_layout(
        xaxis_title=f'{A}<sub>1-x</sub>{B}<sub>x</sub>',
        yaxis_title='ΔH [eV/atom]',
        template='simple_white')

    return fig


def plot_3d(pd: PhaseDiagram,
            uids: list[str] | None = None,
            sources: list[str] | None = None,
            uid: str | None = None) -> go.Figure:

    if uids is None:
        uids = [r[2] for r in pd.references]

    if sources is None:
        sources = ['Materials'] * len(uids)

    x, y, z = pd.points[:, 1:].T
    i, j, k = pd.simplices.T

    data = [go.Mesh3d(x=x, y=y, z=z, i=i, j=j, k=k,
                      opacity=0.3, hoverinfo='skip', color=colors[2])]

    names = [format(Formula(ref[2]).reduce()[0], 'html')
             for ref in pd.references]

    for i, source in enumerate(set(sources)):
        mask = [True if source in label else False
                for label in sources]
        hovertemplate = '%{customdata} <br> ΔH: %{z:.2f} eV/atom'

        data.append(
            go.Scatter3d(
                x=x[mask], y=y[mask], z=z[mask],
                customdata=[f'{name} ({uid})' for uid, name, x
                            in zip(uids, names, mask) if x],
                name=source,
                hovertemplate=hovertemplate,
                mode='markers',
                marker=dict(color=colors[i], size=5),))

    fig = go.Figure(data=data)

    # Highlight selected material:
    if uid is not None:
        this_idx = uids.index(uid)
        names[this_idx] = '<b>' + names[this_idx] + '</b>'

    annotate_idx = [i for i, x in enumerate(pd.hull) if x]
    if uid is not None:
        if this_idx not in annotate_idx:
            annotate_idx.append(this_idx)  # pragma: no cover

    annotations = []
    for i in annotate_idx:
        annotations.append(dict(showarrow=False,
                                x=x[i],
                                y=y[i],
                                z=z[i],
                                text=names[i],
                                xanchor='left',
                                xshift=10,
                                opacity=0.7))

    A, B, C = pd.symbols

    fig.update_layout(scene=dict(xaxis_title=B,
                                 yaxis_title=C,
                                 zaxis_title='ΔH [eV/atom]',
                                 annotations=annotations,))
    return fig


def group_references(references: dict[str, tuple[str, ...]],
                     uids: Iterable[str],
                     check=True) -> dict[tuple[str, ...], list[str]]:
    """Group references into sets of convex hull candidates.

    >>> refs = {'1': ('A',),
    ...         '2': ('B',),
    ...         '3': ('A', 'B'),
    ...         'u1': ('A',),
    ...         'u2': ('A', 'B')}
    >>> group_references(refs, ['u1', 'u2'])
    {('A',): ['1', 'u1'], ('A', 'B'): ['1', '2', '3', 'u1', 'u2']}
    """
    index = defaultdict(set)
    for uid, symbols in references.items():
        if check and sorted(symbols) != list(symbols):
            print(symbols)
            raise ValueError
        for symbol in symbols:
            index[symbol].add(uid)
    chulls = {}
    for uid in uids:
        symbols = references[uid]
        if symbols in chulls:
            continue
        chull = set()
        for symbol in symbols:
            for uid2 in index[symbol]:
                if all(s in symbols for s in references[uid2]):
                    chull.add(uid2)
        chulls[symbols] = sorted(chull)
    return chulls


def calculate_ehull_energies(refs: dict[str, tuple[dict[str, int], float]],
                             uids: set[str], verbose=0) -> dict[str, float]:
    """Calculate energies above hull.

    >>> calculate_ehull_energies(
    ...     {'i1': ({'A': 1}, 0.0),
    ...      'i2': ({'B': 1}, 0.0),
    ...      'i3': ({'A': 1, 'B': 1}, 1.0)},
    ...     {'i3'})
    {'i3': 1.0}
    """
    pd = PhaseDiagram([ref for ref in refs.values()], verbose=verbose)
    ehull_energies = {}
    for uid in refs:
        if uid in uids:
            count, hform = refs[uid]
            ehull_energies[uid] = hform - pd.decompose(**count)[0]
    return ehull_energies


if __name__ == '__main__':
    # Example:
    # python -m camdweb.panels.convex_hull A:0 B:0 AB:-0.5
    refs = []
    for arg in sys.argv[1:]:
        formula, energy = arg.split(':')
        refs.append((formula, float(energy)))
    pd = PhaseDiagram(refs)
    if len(pd.symbols) == 2:
        fig = plot_2d(pd)
    else:
        fig = plot_3d(pd)
    fig.show()
