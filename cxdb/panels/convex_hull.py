r"""
+---------------------+-----------+
|    ^                | OQMR refs |
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

import json
import sys


import plotly
import plotly.graph_objs as go
from ase.formula import Formula
from ase.phasediagram import PhaseDiagram

from cxdb.html import table
from cxdb.material import Material, Materials
from cxdb.panels.asr_panel import read_result_file
from cxdb.panels.panel import Panel

HTML = """
<div class="row">
  <div class="col-6">
    <div id='chull' class='chull'></div>
  </div>
  <div class="col-6">
    {tables}
  </div>
</div>
"""

FOOTER = """
<script type='text/javascript'>
var graphs = {chull_json};
Plotly.newPlot('chull', graphs, {{}});
</script>
"""

OQMD = 'https://cmrdb.fysik.dtu.dk/oqmd123/row'


class ConvexHullPanel(Panel):
    title = 'Convex hull'

    def get_html(self,
                 material: Material,
                 materials: Materials) -> tuple[str, str]:
        result_file = material.folder / 'results-asr.convex_hull.json'
        if not result_file.is_file():
            return '', ''
        data = read_result_file(result_file)
        chull, tbls = make_figure_and_tables(data['references'])
        html = HTML.format(tables=tbls)
        if chull:
            return (html, FOOTER.format(chull_json=chull))
        return html, ''  # pragma: no cover


def make_figure_and_tables(references: list[dict]) -> tuple[str, str]:
    """Make convex-hull figure and tables.

    >>> refs = [
    ...     {'title': 'OQMD', 'hform': 0.0, 'formula': 'B', 'uid': 'u1'},
    ...     {'title': 'OQMD', 'hform': -0.5, 'formula': 'BN', 'uid': 'u2'},
    ...     {'title': 'OQMD', 'hform': 0.0, 'formula': 'N', 'uid': 'u3'},
    ...     {'title': 'C2DB', 'hform': -0.2, 'formula': 'BN', 'uid': 'a1'}]
    >>> chullhtml, tableshtml = make_figure_and_tables(refs)
    Species: B, N
    References: 4
    0    B              0.000
    1    BN            -1.000
    2    N              0.000
    3    BN            -0.400
    Simplices: 2
    """
    tbl1 = []
    tbl2 = []
    pdrefs = []
    labels = []
    for ref in references:
        e = ref['hform']
        f = Formula(ref['formula'])
        uid = ref['uid']
        if 'natoms' in ref:
            assert ref['natoms'] == len(f)
        pdrefs.append((f.count(), e * len(f)))
        if 'OQMD' in ref['title']:
            source = 'OQMD'
            tbl1.append(
                [f'<a href={OQMD}/{uid}>{f:html}</a>', f'{e:.2f} eV/atom'])
        else:
            assert 'C2DB' in ref['title']
            source = 'C2DB'
            tbl2.append(
                [f'<a href={uid}>{f:html}</a>', f'{e:.2f} eV/atom'])
        labels.append(f'{source}({uid})')

    try:
        pd = PhaseDiagram(pdrefs)
    except ValueError:
        chull = ''  # only one species
    else:
        if len(pd.symbols) < 4:
            if len(pd.symbols) == 2:
                fig = plot_2d(pd, labels)
            else:
                fig = plot_3d(pd, labels)
            chull = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
        else:
            chull = ''

    tbls = (
        table(['Bulk crystals from OQMD123', ''], tbl1) +
        table(['Monolayers from C2DB', ''], tbl2))

    return chull, tbls


def plot_2d(pd: PhaseDiagram,
            labels: list[str] | None = None) -> go.Figure:
    if labels is None:
        labels = [r[2] for r in pd.references]

    x, y = pd.points[:, 1:].T

    X = []
    Y = []
    for i, j in pd.simplices:
        X += [x[i], x[j], None]
        Y += [y[i], y[j], None]
    data = [go.Scatter(x=X, y=Y, mode='lines')]

    data.append(go.Scatter(
        x=x,
        y=y,
        text=labels,
        hovertemplate='%{text}: %{y} eV/atom',
        mode='markers'))

    delta = y.ptp() / 30
    ymin = y.min() - 2.5 * delta
    fig = go.Figure(data=data, layout_yaxis_range=[ymin, 0.1])

    A, B = pd.symbols
    fig.update_layout(
        xaxis_title=f'{A}<sub>1-x</sub>{B}<sub>x</sub>',
        yaxis_title='ΔH [eV/atom]',
        template='simple_white')

    return fig


def plot_3d(pd: PhaseDiagram,
            labels: list[str] | None = None) -> go.Figure:
    if labels is None:
        labels = [r[2] for r in pd.references]
    x, y, z = pd.points[:, 1:].T
    i, j, k = pd.simplices.T
    data = [go.Mesh3d(x=x, y=y, z=z, i=i, j=j, k=k, opacity=0.5)]
    data.append(
        go.Scatter3d(
            x=x, y=y, z=z,
            text=labels,
            hovertemplate='%{text}: %{z} eV/atom',
            mode='markers'))

    fig = go.Figure(data=data)
    A, B, C = pd.symbols
    fig.update_layout(scene=dict(xaxis_title=B,
                                 yaxis_title=C,
                                 zaxis_title='ΔH [eV/atom]'),
                      template='simple_white')

    return fig


if __name__ == '__main__':
    # Example:
    # pyhton -m cxdb.panels.convex_hull A:0 B:0 AB:-0.5
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
