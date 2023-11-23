import plotly.graph_objects as go
from cxdb.section import Section

HTML = """
<input type="text" name="repeat" onchange="cb(this.value, 'atoms')">
<div id='atoms' class='atoms'></div>
"""

FOOTER = """
<script type='text/javascript'>
var graphs = {atoms_json};
Plotly.newPlot('atoms', graphs, {});
</script>
"""


class Atoms(Section):
    def __init__(self, atoms):
        self.atoms = atoms
        super().__init__(
            html=HTML, footer=FOOTER)

    def plot(self, repeat: int = 1):
        atoms = self.atoms * repeat
        x, y, z = atoms.positions.T
        fig = go.Figure(
            data=[
                go.Scatter3d(
                    x=x,
                    y=y,
                    z=z,
                    mode='markers')])
        return fig


def plot_atoms(n=5):
    df = pd.DataFrame(
        {'Fruit': ['Apples', 'Oranges', 'Bananas', 'Apples', 'Oranges',
                   'Bananas'],
         'Amount': [4, 1, 2, 2, 4, n],
         'City': ['SF', 'SF', 'SF', 'Montreal', 'Montreal', 'Montreal']})
    fig = px.bar(df, x='Fruit', y='Amount', color='City',
                 barmode='group')
    return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
