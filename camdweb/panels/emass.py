import json
from typing import Generator
from pathlib import Path
from camdweb.panels.panel import Panel
from camdweb.material import Material
from camdweb.html import table
from camdweb.html import image
from asr.effective_masses import make_figure

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
    title = 'Effective masses (PBE)'
    datafiles = ['emass.json']

    def get_html(self,
                 material: Material) -> Generator[str, None, None]:
        path = material.folder / self.datafiles[0]
        with open(path, 'r') as file:
            band_data = json.load(file)
        vbm_data = band_data['vbm']
        cbm_data = band_data['cbm']
        tables = []
        # generate data tables
        for data in [vbm_data, cbm_data]:
            min_emass = ('Min eff. mass', '%.2f m<sub>0</sub>' %
                         data['min_emass'])
            max_emass = ('Max eff. mass', '%.2f m<sub>0</sub>' %
                         data['max_emass'])
            dos_emass = ('DOS eff. mass', '%.2f m<sub>0</sub>' %
                         data['m_dos'])
            coords = ('Crystal coordinates', '[%.3f, %.3f]' % (
                data['coords'][0], data['coords'][1]))
            warping = ('Warping parameter', '%.3f' % data['warping'])
            if data['barrier_found']:
                barrier_height = ('Barrier height', '%.1f meV' %
                                  data['extremum_depth'])
                dist_to_barrier = ('Distance to barrier', '%.3g Å<sup>-1</sup>' %
                                   data['dist_to_barrier'])
            else:
                barrier_height = ('Barrier height', '> %.1f meV' %
                                  data['extremum_depth'])
                dist_to_barrier = ('Distance to barrier', '> %.3g Å<sup>-1</sup>' %
                                   data['dist_to_barrier'])
            new_table = table(['Property (VBM)', 'Value'],
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
            make_figure(data, folder=material.folder)

        html = HTML.format(vbmfig=image(vbmfig, 'VBM'),
                           vbmtable=vbmtable,
                           cbmfig=image(cbmfig, 'CBM'),
                           cbmtable=cbmtable,)
        yield html


"""
            table=table(['#', 'Chemical symbol', 'Charges [|e|]'],
                        [(n, s, f'{c:.2f}') for n, (s, c)
                         in enumerate(zip(material.atoms.symbols,
                                          charges))]))
"""
