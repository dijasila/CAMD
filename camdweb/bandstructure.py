import dataclasses

import numpy as np
import plotly.graph_objs as go
from ase.dft.kpoints import BandPath


@dataclasses.dataclass
class PlotUtil:
    path: BandPath

    energy_nosoc_skn: np.ndarray
    energy_soc_mk: np.ndarray
    spin_zprojection_soc_mk: np.ndarray

    fermilevel: float
    emin: float
    emax: float
    spin_axisname: str

    xcname: str = 'PBE'  # XXX Should not take a default

    def __post_init__(self):
        xcoords, label_xcoords, labels = self.path.get_linear_kpoint_axis()

        self.xcoords = xcoords
        self.kpoint_labels, self.label_xcoords = prettify_labels(
            labels, [*label_xcoords])

        self.axisargs = dict(
            showgrid=True,
            showline=True,
            linewidth=2,
            gridcolor='lightgrey',
            linecolor='black')

    def plot(self) -> go.Figure:
        """Plot band structure.

        This plots:
         * The ordinary (non-spin–orbit-coupled) band structure
         * The spin–orbit coupled band structure coloured by spin projection
         * Reference energy (Fermi level).
        """
        return go.Figure(
            data=[self.plot_bands_boring(),
                  self.plot_bands_fancy(),
                  self.plot_reference_energy_as_line()],
            layout=self.layout())

    def subtract_reference_energy(self, reference):
        return dataclasses.replace(
            self,
            energy_nosoc_skn=self.energy_nosoc_skn - reference,
            energy_soc_mk=self.energy_soc_mk - reference,
            fermilevel=self.fermilevel - reference,
            emin=self.emin - reference,
            emax=self.emax - reference)

    @property
    def xmax(self):
        return self.xcoords[-1]

    def fancymarker_and_also_colorbar(self, color):
        cbtitle = f'〈<i><b>S</b></i><sub>{self.spin_axisname}</sub>〉'

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
                titleside='right'))

    def boringmarker(self):
        return dict(size=4, color='#999999')

    def plot_bands(self, xcoords_k, energies_xk, name, marker):
        assert len(xcoords_k) == energies_xk.shape[-1]

        ndatasets = energies_xk.size // len(xcoords_k)
        xcoords_xk = np.tile(xcoords_k, ndatasets)

        return go.Scattergl(
            x=xcoords_xk.ravel(),
            y=energies_xk.ravel(),
            name=name,
            marker=marker,
            mode='markers',
            showlegend=True,
            hovertemplate='%{y:.3f} eV')

    def plot_bands_boring(self):
        return self.plot_bands(
            self.xcoords, self.energy_nosoc_skn.transpose(0, 2, 1),
            name=f'{self.xcname} no SOC',
            marker=self.boringmarker())

    def plot_bands_fancy(self):
        perm = (-self.spin_zprojection_soc_mk).ravel().argsort()
        esoc_mk = self.energy_soc_mk.ravel()[perm]
        zsoc_mk = self.spin_zprojection_soc_mk.ravel()[perm]
        ndatasets = esoc_mk.size // len(self.xcoords)
        xcoords_mk = np.tile(self.xcoords, ndatasets)[perm]

        return self.plot_bands(
            xcoords_mk, esoc_mk, name=self.xcname,
            marker=self.fancymarker_and_also_colorbar(color=zsoc_mk))

    def plot_reference_energy_as_line(self):
        return go.Scatter(
            x=[0, self.xmax],
            y=[self.fermilevel, self.fermilevel],
            mode='lines',
            line=dict(color=('rgb(0, 0, 0)'), width=2, dash='dash'),
            name='Fermi level')

    def layout(self):
        return go.Layout(
            xaxis=self.xaxis(),
            yaxis=self.yaxis(),
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
                itemwidth=35))

    def xaxis(self):
        return go.layout.XAxis(
            title='k-points',
            range=[0, self.xmax],
            ticks='',
            showticklabels=True,
            mirror=True,
            ticktext=self.kpoint_labels,
            tickvals=self.label_xcoords,
            **self.axisargs)

    def yaxis(self):
        return go.layout.YAxis(
            title='<i>E</i> - <i>E</i><sub>vac</sub> [eV]',
            range=[self.emin, self.emax],
            zeroline=False,
            mirror='ticks',
            ticks='inside',
            tickwidth=2,
            zerolinewidth=2,
            **self.axisargs)


def prettify_labels(orig_labels, label_xcoords):
    import re

    label_xcoords = [*label_xcoords]

    def pretty(kpt):
        if kpt == 'G':
            return 'Γ'

        # Convert Abc123 ----> Abc<sub>123</sub>:
        return re.sub('[0-9]+', lambda match:
                      rf'<sub>{match.group()}</sub>', kpt)

    labels = [pretty(name) for name in orig_labels]

    assert len(labels) == len(label_xcoords)

    i = 1
    while i < len(labels):
        if abs(label_xcoords[i - 1] - label_xcoords[i]) < 1e-6:
            # Merge two special points A and B at same kpoint into a composite
            # label "A,B" and shorten the list of labels:
            labels[i - 1] = f'{labels[i - 1]},{labels[i]}'
            labels.pop(i)
            label_xcoords.pop(i)
        else:
            i += 1

    return labels, label_xcoords
