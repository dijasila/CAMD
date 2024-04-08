import json

from camdweb.panels.panel import Panel, PanelData
from camdweb.material import Material
from camdweb.html import table


HTML = """
<div class="row">
{table}
</div>
"""


class BaderPanel(Panel):
    datafiles = ['bader.json']

    def get_data(self,
                 material: Material) -> PanelData:
        path = material.folder / self.datafiles[0]
        charges = json.loads(path.read_text())['charges']
        html = HTML.format(
            table=table(['#', 'Chemical symbol', 'Charges [|e|]'],
                        [(n, s, f'{c:.2f}') for n, (s, c)
                         in enumerate(zip(material.atoms.symbols,
                                          charges))]))
        return PanelData(html,
                         title='Bader-charge analysis')
