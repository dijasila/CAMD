import json

from cxdb.panels.panel import Panel
from cxdb.material import Material, Materials
from cxdb.html import table


HTML = """
<div class="row">
{table}
</div>
"""


class BaderPanel(Panel):
    title = 'Bader-charge analysis'

    def get_html(self,
                 material: Material,
                 materials: Materials) -> tuple[str, str]:
        path = material.folder / 'bader.json'
        if not path.is_file():
            return ('', '')
        charges = json.loads(path.read_text())['charges']
        return (
            HTML.format(
                table=table(['#', 'Chemical symbol', 'Charges [|e|]'],
                            [(n, s, f'{c:.2f}') for n, (s, c)
                             in enumerate(zip(material.atoms.symbols,
                                              charges))])),
            '')
