import json

import numpy as np
import plotly
import plotly.graph_objs as go

from camdweb.panels.panel import Panel


HTML = """
<div class="row">
  <div class="col-6">
    <div id='{plot_name}' class='{plot_name}'></div>
  </div>
</div>

<script type='text/javascript'>
var graphs = {plot_data};
Plotly.newPlot('{plot_name}', graphs, {{}});
</script>
"""


OQMD = 'https://cmrdb.fysik.dtu.dk/oqmd123/row'


class BandStructurePanel(Panel):
    title = 'Band structure'

    def get_html(self, material, materials):
        from camdweb.c2db.asr_panel import Row
        row = Row(material)
        fig = plot_bs_html(row)
        bandstructure_json = json.dumps(
            fig, cls=plotly.utils.PlotlyJSONEncoder)
        yield HTML.format(plot_data=bandstructure_json,
                          plot_name='bandstructure')


def make_plot():
    import numpy as np
    x = np.linspace(-5, 5, 100)
    y = np.sin(x)

    data = [go.Scatter(x=x, y=y, mode='lines')]
    fig = go.Figure(data=data)
    bandstructure_json = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    return bandstructure_json



def plot_bs_html(row):  # , filename=bs_html):
    import plotly
    import plotly.graph_objs as go

    from ase.dft.kpoints import labels_from_kpts

    traces = []
    d = row.data.get('results-asr.bandstructure.json')
    xcname = 'PBE'
    # xcname = gs_xcname_from_row(row)

    path = d['bs_nosoc']['path']
    kpts = path.kpts
    ef = d['bs_nosoc']['efermi']

    if row.get('evac') is not None:
        label = '<i>E</i> - <i>E</i><sub>vac</sub> [eV]'
        reference = row.get('evac')
    else:
        label = '<i>E</i> - <i>E</i><sub>F</sub> [eV]'
        reference = ef

    gaps = row.data.get('results-asr.gs.json', {}).get('gaps_nosoc', {})
    if gaps.get('vbm'):
        emin = gaps.get('vbm', ef) - 3
    else:
        emin = ef - 3
    if gaps.get('cbm'):
        emax = gaps.get('cbm', ef) + 3
    else:
        emax = ef + 3
    e_skn = d['bs_nosoc']['energies']
    shape = e_skn.shape
    xcoords, label_xcoords, orig_labels = labels_from_kpts(
        kpts, row.cell, special_points=path.special_points
    )
    xcoords = np.vstack([xcoords] * shape[0] * shape[2])
    # colors_s = plt.get_cmap('viridis')([0, 1])  # color for sz = 0
    e_kn = np.hstack([e_skn[x] for x in range(shape[0])])
    trace = go.Scattergl(
        x=xcoords.ravel(),
        y=e_kn.T.ravel() - reference,
        mode='markers',
        name=f'{xcname} no SOC',
        showlegend=True,
        marker=dict(size=4, color='#999999'),
        hovertemplate='%{y:.3f} eV'
    )
    traces.append(trace)

    e_mk = d['bs_soc']['energies']
    path = d['bs_soc']['path']
    kpts = path.kpts
    ef = d['bs_soc']['efermi']
    sz_mk = d['bs_soc']['sz_mk']

    xcoords, label_xcoords, orig_labels = labels_from_kpts(
        kpts, row.cell, special_points=path.special_points
    )

    shape = e_mk.shape
    perm = (-sz_mk).argsort(axis=None)
    e_mk = e_mk.ravel()[perm].reshape(shape)
    sz_mk = sz_mk.ravel()[perm].reshape(shape)
    xcoords = np.vstack([xcoords] * shape[0])
    xcoords = xcoords.ravel()[perm].reshape(shape)

    # Unicode for <S_z>
    sdir = row.get('spin_axis', 'z')
    cbtitle = '&#x3008; <i><b>S</b></i><sub>{}</sub> &#x3009;'.format(sdir)
    trace = go.Scattergl(
        x=xcoords.ravel(),
        y=e_mk.ravel() - reference,
        mode='markers',
        name=xcname,
        showlegend=True,
        marker=dict(
            size=4,
            color=sz_mk.ravel(),
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(
                tickmode='array',
                tickvals=[-1, 0, 1],
                ticktext=['-1', '0', '1'],
                title=cbtitle,
                titleside='right',
            ),
        ),
        hovertemplate='%{y:.3f} eV'
    )
    traces.append(trace)

    linetrace = go.Scatter(
        x=[np.min(xcoords), np.max(xcoords)],
        y=[ef - reference, ef - reference],
        mode='lines',
        line=dict(color=('rgb(0, 0, 0)'), width=2, dash='dash'),
        name='Fermi level',
    )
    traces.append(linetrace)

    def pretty(kpt):
        if kpt == 'G':
            kpt = '&#x393;'  # Gamma in unicode
        elif len(kpt) == 2:
            kpt = kpt[0] + '$_' + kpt[1] + '$'
        return kpt

    labels = [pretty(name) for name in orig_labels]
    i = 1
    while i < len(labels):
        if label_xcoords[i - 1] == label_xcoords[i]:
            labels[i - 1] = labels[i - 1][:-1] + ',' + labels[i][1:]
            labels[i] = ''
        i += 1

    bandxaxis = go.layout.XAxis(
        title='k-points',
        range=[0, np.max(xcoords)],
        showgrid=True,
        showline=True,
        ticks='',
        showticklabels=True,
        mirror=True,
        linewidth=2,
        ticktext=labels,
        tickvals=label_xcoords,
        gridcolor='lightgrey',
        linecolor='black',
    )

    bandyaxis = go.layout.YAxis(
        title=label,
        range=[emin - reference, emax - reference],
        showgrid=True,
        showline=True,
        zeroline=False,
        mirror='ticks',
        ticks='inside',
        linewidth=2,
        tickwidth=2,
        zerolinewidth=2,
        gridcolor='lightgrey',
        linecolor='black',
    )

    bandlayout = go.Layout(
        xaxis=bandxaxis,
        yaxis=bandyaxis,
        plot_bgcolor='white',
        hovermode='closest',
        margin=dict(t=20, r=20),
        font=dict(size=14),
        legend=dict(
            orientation='h',
            yanchor='bottom',
            y=1.01,
            xanchor='left',
            x=0.0,
            font=dict(size=14),
            itemsizing='constant',
            itemwidth=35
        ),
    )

    return go.Figure(data=traces, layout=bandlayout)
