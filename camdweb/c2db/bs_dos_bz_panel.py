import json

import plotly
from asr.pdos import plot_pdos
from asr.gs import bz_with_band_extremums

from camdweb.panels.panel import Panel
from camdweb.c2db.asr_panel import Row
from camdweb.panels.bandstructure import plotter_from_row
from camdweb.html import image

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


class BSDOSBZPanel(Panel):
    title = 'Band structure and density of states'

    def get_html(self, material):
        bs = material.folder / 'results-asr.bandstructure.json'
        if not (bs.is_file() or bs.with_suffix('.json.gz').is_file()):
            return
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
        yield HTML.format(bs_data=bs_json,
                          dos=image(dos_file, 'DOS'),
                          bz=image(bz_file, 'BZ'))
