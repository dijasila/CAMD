import pandas as pd
import json

from cxdb.section import Section
from cxdb.material import Material


class BaderSection(Section):
    title = 'Bader-charge analysis'

    def get_html(self, material: Material) -> tuple[str, str]:
        path = material.folder / 'bader.json'
        if not path.is_file():
            return ('', '')
        charges = json.loads(path.read_text())['charges']
        df = pd.DataFrame(data={'Chemical symbol': material.atoms.symbols,
                                'Charges': charges})
        return (df.to_html(), '')
