import numpy as np
from pathlib import Path

from cxdb.panels.asr_panel import read_result_file
from cxdb.material import Material, Materials
from cxdb.panels.panel import Panel

from ase.phasediagram import PhaseDiagram
from ase.formula import Formula
import plotly.graph_objs as go
import plotly.express as px
import pandas

HTML = """
<img alt="DOS for {uid}" src="/png/{uid}/dos.png" />
"""


class ConvexHullPanel(Panel):
    title = 'Convex hull'

    def get_html(self, material: Material, materials: Materials) -> tuple[str, str]:
        result_file = material.folder / 'results-asr.convex_hull.json'
        n_elements = len(set(material.atoms.get_atomic_numbers()))
        uid = material.uid
        if not result_file.is_file():
            return ('', '')
        self.make_figures(result_file, n_elements, uid)
        return (HTML.format(uid=material.uid), '')

    def make_figures(self, result_file: Path, n_elements: int, uid: str):
        data = read_result_file(result_file)
        plot_convex_hull(data)


def plot_convex_hull(data: dict, n_elements: int, uid: str) -> go.Figure:
    colors = px.colors.qualitative.D3

    if not (2 <= n_elements <= 3):
        return

    references = data['references']

    df_ref = pandas.DataFrame(references)

    df_ref['thisrow'] = df_ref.apply(
        lambda x: True if x['uid'] == uid else False, axis=1
    )

    names = [ref['label'] for ref in references]
    latexnames = [
        format(Formula(name.split(' ')[0]).reduce()[0], 'html') for name in names
    ]

    df_ref['latexname'] = latexnames

    # Highlight this material by making it bold
    name_column_to_plot = 'latexname'
    try:
        df_ref.loc[df_ref['thisrow'], name_column_to_plot] = (
            '<b>' + df_ref[df_ref['thisrow']][name_column_to_plot].values[0] + '</b>'
        )
    except IndexError:
        pass

    pdrefs = []

    for reference in references:
        h = reference['natoms'] * reference['hform']
        pdrefs.append((reference['formula'], h))

    pd = PhaseDiagram(pdrefs, verbose=False)

    if n_elements == 2:
        fig = plot_2D(pd)
    else:
        fig = plot_3D(df_ref, pd, colors)

    return fig


def plot_2D(pd):
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
        text=[name for _, _, name, _ in pd.references],
        hovertemplate=' %{text}: %{y} eV/atom',
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


def plot_3D(df_ref, pd, colors):
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
    fig_temp = px.scatter_3d(
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


if __name__ == '__main__':
    # Example data
    data = {
        'hform': -0.920927544752896,
        'references': [
            {
                'hform': 0.0,
                'formula': 'S48',
                'uid': 'S48',
                'natoms': 48,
                'title': 'Bulk crystals (from OQMD123)',
                'legend': 'Bulk crystals',
                'name': 'S48',
                'label': 'S48',
                'link': 'https://cmrdb.fysik.dtu.dk/oqmd123/row/S48',
                'method': 'DFT',
            },
            {
                'hform': 0.0,
                'formula': 'Mo',
                'uid': 'Mo',
                'natoms': 1,
                'title': 'Bulk crystals (from OQMD123)',
                'legend': 'Bulk crystals',
                'name': 'Mo',
                'label': 'Mo',
                'link': 'https://cmrdb.fysik.dtu.dk/oqmd123/row/Mo',
                'method': 'DFT',
            },
            {
                'hform': -0.18018421551178587,
                'formula': 'Mo2S2',
                'uid': 'Mo2S2-925d20f42e31',
                'natoms': 4,
                'title': 'Monolayers (from C2DB)',
                'legend': 'Monolayers',
                'name': 'Mo2S2 (AB-187-hi)',
                'label': 'Mo2S2 (AB-187-hi)',
                'link': '/c2db/row/Mo2S2-925d20f42e31',
                'method': 'DFT',
            },
            {
                'hform': -0.920927544752896,
                'formula': 'MoS2',
                'uid': 'MoS2-b3b4685fb6e1',
                'natoms': 3,
                'title': 'Monolayers (from C2DB)',
                'legend': 'Monolayers',
                'name': 'MoS2 (AB2-187-bi)',
                'label': 'MoS2 (AB2-187-bi)',
                'link': '/c2db/row/MoS2-b3b4685fb6e1',
                'method': 'DFT',
            },
        ],
        'indices': [7, 2],
        'coefs': [1.0, 0.0],
        'ehull': 0.0,
        'thermodynamic_stability_level': 3,
    }
    print(data)
    fig = plot_convex_hull(data, 2, 'MoS2-b3b4685fb6e1')
    fig.show()
