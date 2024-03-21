import json
from pathlib import Path
from camdweb.panels.panel import Panel, PanelData
from camdweb.material import Material
from camdweb.html import table
from camdweb.html import image
import numpy as np

HTML = """
<div class="row">
  <div class="col-6">
   {vbmfig}
   {vbmtable}
  </div>
  <div class="col-6">
   {cbmfig}
   {cbmtable}
  </div>
</div>
"""


class EmassPanel(Panel):
    datafiles = ['emass.json']

    def get_data(self,
                 material: Material) -> PanelData:
        path = material.folder / self.datafiles[0]
        with open(path, 'r') as file:
            band_data = json.load(file)

        def emass2html(emass):
            if emass == np.inf:
                emass_str = '∞'
            else:
                emass_str = f'{emass:.2f} m<sub>0</sub>'
            return emass_str

        tables = []
        # generate data tables
        for band_name in ['vbm', 'cbm']:
            data = band_data[band_name]
            min_emass = ('Min eff. mass', emass2html(data['min_emass']))
            max_emass = ('Max eff. mass', emass2html(data['max_emass']))
            dos_emass = ('DOS eff. mass', emass2html(data['m_dos']))
            coords = ('Crystal coordinates', '[%.3f, %.3f]' % (
                data['coords'][0], data['coords'][1]))
            warping = ('Warping parameter', '%.3f' % data['warping'])
            if data['barrier_found']:
                barrier_height = ('Barrier height', '%.1f meV' %
                                  data['extremum_depth'])
                dist_to_barrier = ('Distance to barrier',
                                   '%.3g Å<sup>-1</sup>' %
                                   data['dist_to_barrier'])
            else:
                barrier_height = ('Barrier height', '> %.1f meV' %
                                  data['extremum_depth'])
                dist_to_barrier = ('Distance to barrier',
                                   '> %.3g Å<sup>-1</sup>' %
                                   data['dist_to_barrier'])
            new_table = table(['Property (' + band_name.upper() + ')',
                               'Value'],
                              [min_emass,
                               max_emass,
                               dos_emass,
                               coords,
                               warping,
                               barrier_height,
                               dist_to_barrier])
            tables.append(new_table)
        # columns = p['columns']
        vbmtable = tables[0]
        cbmtable = tables[1]
        vbmfig = material.folder / 'emass_vbm.png'
        cbmfig = material.folder / 'emass_cbm.png'
        if not (vbmfig.is_file() and cbmfig.is_file()):
            make_figure(band_data, folder=material.folder)

        html = HTML.format(vbmfig=image(vbmfig, 'VBM'),
                           vbmtable=vbmtable,
                           cbmfig=image(cbmfig, 'CBM'),
                           cbmtable=cbmtable)
        return PanelData(html, title='Effective masses (PBE)')


def make_figure(data, folder: Path):
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.colors import Normalize as pltnorm
    band_names = ['vbm', 'cbm']  # hard coded band-names. change in future
    for i, band_name in enumerate(band_names):
        band_data = data[band_name]
        for key in band_data:
            band_data[key] = np.asarray(band_data[key])

        X = band_data['X']
        Y = band_data['Y']
        Z = band_data['Z']
        min_emass_direction = band_data['min_emass_direction']
        max_emass_direction = band_data['max_emass_direction']

        line_radius = np.sqrt(X**2 + Y**2).max() / 2
        line = np.linspace(0, line_radius, 101)

        fig, ax = plt.subplots()
        ax.contourf(X, Y, Z, cmap='viridis',
                    levels=80, vmin=Z.min(), vmax=Z.max())

        diff = Z.max() - Z.min()
        if band_name == 'cbm':
            extra_contours = Z.min()\
                + (np.linspace(0, np.sqrt(diff), 10)[1:-1])**2
            ax.contour(X, Y, Z, cmap='viridis', levels=extra_contours,
                       vmin=Z.min() - 0.1 * diff, vmax=Z.max() + 0.1 * diff)
        elif band_name == 'vbm':
            extra_contours = np.flip(
                Z.max() - (np.linspace(0, np.sqrt(diff), 10)[1:-1])**2)
            ax.contour(X, Y, Z, cmap='viridis', levels=extra_contours,
                       vmin=Z.min() - 0.1 * diff, vmax=Z.max() + 0.1 * diff)

        ax.set_xlabel(r'$k_x$ / Å$^{-1}$')
        ax.set_ylabel(r'$k_y$ / Å$^{-1}$')
        cbar = fig.colorbar(cm.ScalarMappable(cmap='viridis', norm=pltnorm(
            1000 * Z.min(), 1000 * Z.max())), ax=ax)
        cbar.set_label(r'$(E - E_0)$ / meV')

        # add circles to contour plot
        ax.plot(line * max_emass_direction[0],
                line * max_emass_direction[1], color='tab:green', ls='dashed',
                label='Max eff. mass direction')
        ax.plot(line * min_emass_direction[0],
                line * min_emass_direction[1], color='tab:orange', ls='dashed',
                label='Min eff. mass direction')

        ax.set_xlim(X.min(), X.max())
        ax.set_ylim(Y.min(), Y.max())
        if band_name == 'cbm':
            plt.title('Conduction band minimum (CBM)')
        elif band_name == 'vbm':
            plt.title('Valence band maximum (VBM)')
        else:
            plt.title(band_name)

        plt.tight_layout()

        if len(ax.get_xticks()) >= 7:
            ax.set_xticklabels(np.round(ax.get_xticks(), 4), rotation=15)
        ax.legend()
        filename = 'emass_' + band_name + '.png'
        plt.savefig(folder / filename)
