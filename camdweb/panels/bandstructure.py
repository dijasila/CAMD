import json

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
        yield HTML.format(plot_data=make_plot(),
                          plot_name='bandstructure')


def make_plot():
    import numpy as np
    x = np.linspace(-5, 5, 100)
    y = np.sin(x)

    data = [go.Scatter(x=x, y=y, mode='lines')]
    fig = go.Figure(data=data)
    bandstructure_json = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    return bandstructure_json
