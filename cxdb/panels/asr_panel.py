"""Hack to use webpanel() functions from ASR."""
import importlib
from pathlib import Path

import matplotlib.pyplot as plt
from ase.db.core import KeyDescription
from ase.io.jsonio import decode
from cxdb.html import table
from cxdb.material import Material, Materials
from cxdb.panels.panel import Panel

HTML = """
<h4>{title}</h4>
<div class="row">
  <div class="col-6">
   {col1}
  </div>
  <div class="col-6">
   {col2}
  </div>
</div>
"""


def read_result_file(path: Path) -> dict:
    dct = decode(path.read_text())
    if 'kwargs' in dct:
        dct = dct['kwargs']['data']
    return dct


class Row:
    """Fake row object."""
    def __init__(self, material: Material):
        self.data = Data(material.folder)
        self.__dict__.update(material.columns)
        self.atoms = material.atoms
        self.cell = self.atoms.cell
        self.pbc = self.atoms.pbc
        self.symbols = self.atoms.symbols

    def toatoms(self):
        return self.atoms

    def get(self, name: str, default=None):
        if hasattr(self, name):
            return getattr(self, name)
        print('MISSING:', name, default)
        return default


class Data:
    def __init__(self, folder: Path):
        self.folder = folder

    def get(self, name, default=None):
        dct = decode((self.folder / name).read_text())
        if 'kwargs' in dct:
            return dct['kwargs']['data']
        return dct

    def __contains__(self, name):
        return False

    def __getitem__(self, name):
        return self.get(name)


class ASRPanel(Panel):
    """Generic ASR-panel."""
    def __init__(self, name: str, keys: set[str]):
        self.name = name
        mod = importlib.import_module(f'asr.{name}')
        self.webpanel = mod.webpanel
        self.result_class = mod.Result
        self.key_descriptions = {}
        self.column_names = {}
        for key, desc in getattr(self.result_class,
                                 'key_descriptions', {}).items():
            self.key_descriptions[key] = KeyDescription(key, desc)
            if key in keys:
                self.column_names[key] = desc

    def get_html(self,
                 material: Material,
                 materials: Materials) -> tuple[str, str]:
        """Create row and result objects and call webpanel() function."""
        print('+++++++', self.name)
        row = Row(material)
        try:
            dct = row.data.get(f'results-asr.{self.name}.json')
        except FileNotFoundError:
            return ('', '')
        result = self.result_class(dct)
        (p,) = self.webpanel(result, row, self.key_descriptions)
        plt.close()

        columns: list[list[str]] = [[], []]
        for i, column in enumerate(p['columns']):
            for thing in column:
                # print(thing)
                if thing is not None:
                    html = thing2html(thing, material.folder)
                    columns[i].append(html)

        for desc in p.get('plot_descriptions', []):
            paths = [material.folder / filename
                     for filename in desc['filenames']]
            for f in paths:
                if not f.is_file():
                    # Call plot-function:
                    desc['function'](row, *paths)
                    break

        self.title = p['title']
        return (HTML.format(title=p['title'],
                            col1='\n'.join(columns[0]),
                            col2='\n'.join(columns[1])),
                '')


def thing2html(thing: dict, path: Path) -> str:
    """Convert webpanel() output to HTML."""
    if thing['type'] == 'figure':
        filename = thing['filename']
        html = f'<img src="/png/{path}/{filename}" />'
    elif thing['type'] == 'table':
        html = table(thing.get('header'), thing['rows'])
    else:
        raise ValueError
    return html
