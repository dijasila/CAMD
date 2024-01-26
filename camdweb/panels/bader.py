import json
from typing import Generator

from camdweb.panels.panel import Panel
from camdweb.material import Material, Materials
from camdweb.html import table


HTML = """
<div class="row">
{table}
</div>
"""


class BaderPanel(Panel):
    title = 'Bader-charge analysis'

    def get_html(self,
                 material: Material,
                 materials: Materials) -> Generator[str, None, None]:
        path = material.folder / 'bader.json'
        if not path.is_file():
            return
        charges = json.loads(path.read_text())['charges']
        yield HTML.format(
            table=table(['#', 'Chemical symbol', 'Charges [|e|]'],
                        [(n, s, f'{c:.2f}') for n, (s, c)
                         in enumerate(zip(material.atoms.symbols,
                                          charges))]))
