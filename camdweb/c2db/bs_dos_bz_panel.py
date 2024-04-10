import json

import plotly
from asr.pdos import plot_pdos
from asr.gs import bz_with_band_extremums
import numpy as np

from camdweb.panels.panel import Panel, PanelData, SkipPanel
from camdweb.c2db.asr_panel import Row
from camdweb.html import image
from camdweb.bandstructure import PlotUtil

HTML = """
<div class="row">
  <div class="col-6">
    <div id='bandstructure'></div>
  </div>
  <div class="col-6">
    {dos}
    {bz}
  </div>
</div>

<script type='text/javascript'>
var graphs = {bs_data};
Plotly.newPlot('bandstructure', graphs, {{}});
</script>
"""


def plotter_from_row(row):
    dct = row.data.get('results-asr.bandstructure.json')
    gaps = row.data['results-asr.gs.json']['gaps_nosoc']
    fermilevel_soc = dct['bs_soc']['efermi']

    assert np.allclose(dct['bs_soc']['path'].kpts,
                       dct['bs_nosoc']['path'].kpts)

    vbm = gaps['vbm']
    cbm = gaps['cbm']
    return PlotUtil(
        energy_soc_mk=dct['bs_soc']['energies'],
        energy_nosoc_skn=dct['bs_nosoc']['energies'],
        spin_zprojection_soc_mk=dct['bs_soc']['sz_mk'],
        path=dct['bs_nosoc']['path'],
        fermilevel=fermilevel_soc,
        emin=(fermilevel_soc if vbm is None else vbm) - 3,
        emax=(fermilevel_soc if cbm is None else cbm) + 3,
        spin_axisname=row.get('spin_axis', 'z')  # XXX crazy to have a default
    ).subtract_reference_energy(row.get('evac'))


class BSDOSBZPanel(Panel):
    def get_data(self, material):
        bs = material.folder / 'results-asr.bandstructure.json'
        if not (bs.is_file() or bs.with_suffix('.json.gz').is_file()):
            raise SkipPanel
        row = Row(material)
        plotter = plotter_from_row(row)
        fig = plotter.plot()

        dos_file = material.folder / 'dos.png'
        if not dos_file.is_file():
            plot_pdos(row, dos_file, soc=False)

        bz_file = material.folder / 'bz-with-gaps.png'
        if not bz_file.is_file():
            bz_with_band_extremums(row, bz_file)

        bs_json = json.dumps(
            fig, cls=plotly.utils.PlotlyJSONEncoder)
        return PanelData(
            HTML.format(bs_data=bs_json,
                        dos=image(dos_file, 'DOS'),
                        bz=image(bz_file, 'BZ')),
            title='Electronic band structure and density of states')
