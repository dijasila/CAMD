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


def plot_bs_html(row):
    return PlotUtil(row).plot()


class PlotUtil:
    def __init__(self, row):
        self.row = row
        self.dct = self.row.data.get('results-asr.bandstructure.json')
        self.path = self.dct['bs_nosoc']['path']
        self.kpts = self.path.kpts
        self.evac = row.get('evac')
        assert self.evac is not None
        self.fermilevel_nosoc = self.dct['bs_nosoc']['efermi']
        self.fermilevel_soc = self.dct['bs_soc']['efermi']
        self.e_skn = self.dct['bs_nosoc']['energies']
        self.xcname = 'PBE'

        self.gaps = row.data.get(
            'results-asr.gs.json', {}).get('gaps_nosoc', {})

        self.emin = self.gaps.get('vbm', self.fermilevel_nosoc) - 3
        self.emax = self.gaps.get('cbm', self.fermilevel_nosoc) + 3

        assert np.allclose(self.dct['bs_soc']['path'].kpts,
                           self.dct['bs_nosoc']['path'].kpts)

        (self.xcoords, self.label_xcoords,
         self.orig_labels) = self.path.get_linear_kpoint_axis()

        self.xmax = self.xcoords[-1]
        assert all(self.xmax >= self.xcoords)

        self.esoc_mk = self.dct['bs_soc']['energies']
        self.zsoc_mk = self.dct['bs_soc']['sz_mk']

        self.axisargs = dict(
            showgrid=True,
            showline=True,
            linewidth=2,
            gridcolor='lightgrey',
            linecolor='black',
        )
        self.boringmarker = dict(size=4, color='#999999')

    def fancymarker_and_also_colorbar(self, color):
        sdir = self.row.get('spin_axis', 'z')
        cbtitle = f'〈<i><b>S</b></i><sub>{sdir}</sub>〉'

        return dict(
            size=4,
            color=color,
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(
                tickmode='array',
                tickvals=[-1, 0, 1],
                ticktext=['-1', '0', '1'],
                title=cbtitle,
                titleside='right',
            ),
        )

    def plot_bands(self, xcoords_k, energies_xk, name, marker):
        assert len(xcoords_k) == energies_xk.shape[-1]

        ndatasets = energies_xk.size // len(xcoords_k)
        xcoords_xk = np.tile(xcoords_k, ndatasets)

        return go.Scattergl(
            x=xcoords_xk.ravel(),
            y=energies_xk.ravel() - self.evac,
            name=name,
            marker=marker,
            mode='markers',
            showlegend=True,
            hovertemplate='%{y:.3f} eV')

    def plot_bands_boring(self):
        return self.plot_bands(
            self.xcoords, self.e_skn.transpose(0, 2, 1),
            name=f'{self.xcname} no SOC',
            marker=self.boringmarker)

    def plot_bands_fancy(self):
        perm = (-self.zsoc_mk).ravel().argsort()
        esoc_mk = self.esoc_mk.ravel()[perm]
        zsoc_mk = self.zsoc_mk.ravel()[perm]
        ndatasets = esoc_mk.size // len(self.xcoords)
        xcoords_mk = np.tile(self.xcoords, ndatasets)[perm]

        return self.plot_bands(
            xcoords_mk, esoc_mk, name=self.xcname,
            marker=self.fancymarker_and_also_colorbar(color=zsoc_mk))

    def plot(self):
        return go.Figure(
            data=[self.plot_bands_boring(),
                  self.plot_bands_fancy(),
                  self.plot_reference_energy_as_line()],
            layout=self.bandlayout())

    def plot_reference_energy_as_line(self):
        line_position = self.fermilevel_soc - self.evac
        return go.Scatter(
            x=[0, self.xmax],
            y=[line_position, line_position],
            mode='lines',
            line=dict(color=('rgb(0, 0, 0)'), width=2, dash='dash'),
            name='Fermi level',
        )

    def bandlayout(self):
        return go.Layout(
            xaxis=self.bandxaxis(),
            yaxis=self.bandyaxis(),
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

    def bandxaxis(self):
        return go.layout.XAxis(
            title='k-points',
            range=[0, self.xmax],
            ticks='',
            showticklabels=True,
            mirror=True,
            ticktext=prettify_labels(self.orig_labels, self.label_xcoords),
            tickvals=self.label_xcoords,
            **self.axisargs)

    def bandyaxis(self):
        return go.layout.YAxis(
            title='<i>E</i> - <i>E</i><sub>vac</sub> [eV]',
            range=[self.emin - self.evac, self.emax - self.evac],
            zeroline=False,
            mirror='ticks',
            ticks='inside',
            tickwidth=2,
            zerolinewidth=2,
            **self.axisargs)


def prettify_labels(orig_labels, label_xcoords):

    def pretty(kpt):
        if kpt == 'G':
            kpt = 'Γ'
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

    return labels
