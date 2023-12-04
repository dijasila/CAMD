import json

from cxdb.section import Section
from cxdb.material import Material
from cxdb.html import table


class BaderSection(Section):
    title = 'Bader-charge analysis'

    def get_html(self, material: Material) -> tuple[str, str]:
        path = material.folder / 'bader.json'
        if not path.is_file():
            return ('', '')
        charges = json.loads(path.read_text())['charges']
        return table(['#', 'Chemical symbol', 'Charges [|e|]'],
                     [(n, s, c) for n, (s, c)
                      in enumerate(material.atoms.symbols, charges)]), ''
