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
from pathlib import Path

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
        chull, tbls = self.make_figure_and_tables(result_file)
        html = HTML.format(tables=tbls)
        if chull:
            return (html, FOOTER.format(chull_json=chull))
        return html, ''

    def make_figure_and_tables(self,
                               result_file: Path) -> tuple[str, str]:
        data = read_result_file(result_file)

        tbl1 = []
        tbl2 = []
        references = []
        labels = []
        for ref in data['references']:
            e = ref['hform']
            f = Formula(ref['formula'])
            uid = ref['uid']
            references.append((f.count(), e * ref['natoms']))
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

        pd = PhaseDiagram(references)
        if 2 <= len(pd.symbols) <= 3:
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
    """
    x, y, e = pd.points[:, 1:].T
    df_ref['x'] = x
    df_ref['y'] = y
    df_ref['e'] = e

    hull = np.array(pd.hull)
    df_ref['hull'] = hull

    figs = []
    for i, j, k in pd.simplices:
        fig_temp = go.Figure(
            data=[
                go.Mesh3d(
                    x=x[[i, j, k, i]],
                    y=y[[i, j, k, i]],
                    z=e[[i, j, k, i]],
                    color=colors[2],
                    opacity=0.5,
                ),
            ]
        )
        fig_temp.update_traces(hoverinfo='skip')
        figs.append(fig_temp)

    # Plot materials
    hover_data = {
        'x': False,
        'y': False,
        'e': ':.2f',
        'legend': False,
        'latexname': True,
        'uid': True,
        'name': True,
    }
    fig_temp = go.scatter_3d(
        df_ref,
        x='x',
        y='y',
        z='e',
        hover_data=hover_data,
        color='legend',
        custom_data=['link'],
        color_discrete_sequence=colors,
        labels={
            'x': pd.symbols[1],
            'y': pd.symbols[2],
            'e': '\u0394H [eV/atom]',
            'latexname': 'Formula',
        },
    )
    fig_temp.update_traces(marker={'size': 6})
    figs.append(fig_temp)

    delta = e.ptp() / 30
    ymin = e.min() - 2.5 * delta
    fig = go.Figure(data=sum([fig.data for fig in figs], ()))

    #  Highlight materials on the hull with formula and thisrow
    materials_with_text = df_ref[df_ref.hull | df_ref.thisrow]
    annotations = []
    for row in materials_with_text.itertuples(index=False):
        annotations.append(
            dict(
                showarrow=False,
                x=row.x,
                y=row.y,
                z=row.e,
                text=row.latexname,
                xanchor='left',
                xshift=10,
                opacity=0.7,
            )
        )

    fig.update_layout(
        scene=dict(
            xaxis_title=pd.symbols[1],
            yaxis_title=pd.symbols[2],
            zaxis_title='\u0394H [eV/atom]',
            zaxis=dict(range=[ymin, 0.1]),
            annotations=annotations,
            aspectratio={'x': 1, 'y': 1, 'z': 1},
        ),
        margin=dict(l=0, r=0, b=0, t=0),
        legend=dict(
            orientation='h',
            entrywidth=100,
            yanchor='bottom',
            y=0.9,
            xanchor='left',
            x=0.01,
            font=dict(size=14),
        ),
    )
    return fig
    """


if __name__ == '__main__':
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
