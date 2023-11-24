import json

import plotly
import plotly.graph_objects as go

from cxdb.section import Section
from cxdb.material import Material

HTML = """
<h4>{formula}</h4>
<input
  type="text"
  name="repeat"
  onchange="cb(this.value, 'atoms', '{id}')"
  placeholder="repeat">
<div id='atoms' class='atoms'></div>
"""

FOOTER = """
<script type='text/javascript'>
var graphs = {atoms_json};
Plotly.newPlot('atoms', graphs, {{}});
</script>
"""


class AtomsSection(Section):
    title = 'Atoms'

    def __init__(self):
        self.callbacks = {'atoms': self.plot}

    def get_html(self, material: Material) -> tuple[str, str]:
        return (HTML.format(id=material.id,
                            formula=material.formula_html),
                FOOTER.format(atoms_json=self.plot(material, 1)))

    def plot(self, material: Material, repeat: int = 1) -> str:
        assert repeat < 5, 'DOS!'
        atoms = material.atoms * repeat
        x, y, z = atoms.positions.T
        fig = go.Figure(data=[go.Scatter3d(x=x, y=y, z=z,
                                           mode='markers')])
        return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
