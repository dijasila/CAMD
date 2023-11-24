import pandas as pd
import json


class BaderSection:
    def __init__(self):
        self.callbacks = {}

    def get_html(self, material):
        path = material.folder / 'bader.json'
        if not path.is_file():
            return ('', '', '')
        charges = json.loads(path.read_text())['charges']
        df = pd.DataFrame(data={'Chemical symbol': material.atoms.symbols,
                                'Charges': charges})
        return ('Bader-charge analysis',
                df.to_html(),
                '')
