"""Hack to use webpanel() functions from ASR."""
import importlib
from pathlib import Path

from ase.db.core import KeyDescription
from ase.io.jsonio import decode

from cxdb.material import Material, Materials
from cxdb.panels.panel import Panel
from cxdb.utils import table

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
        self.magstate = material.magstate
        self.has_inversion_symmetry = material.has_inversion_symmetry
        self.cell = material.atoms.cell
        self.minhessianeig = 117.0
        self.evac = material.evac

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
    def __init__(self, name: str):
        self.name = name
        mod = importlib.import_module(f'asr.{name}')
        self.webpanel = mod.webpanel
        self.result_class = mod.Result
        self.key_descriptions = {
            key: KeyDescription(key, desc)
            for key, desc in self.result_class.key_descriptions.items()}

    def get_html(self,
                 material: Material,
                 materials: Materials) -> tuple[str, str]:
        """Create row and result objects and call webpanel() function."""
        uid = material.uid
        row = Row(material)
        try:
            dct = row.data.get(f'results-asr.{self.name}.json')
        except FileNotFoundError:
            return ('', '')
        result = self.result_class(dct)
        (p,) = self.webpanel(result, row, self.key_descriptions)

        columns: list[list[str]] = [[], []]
        for i, column in enumerate(p['columns']):
            for thing in column:
                html = thing2html(thing, uid)
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


def thing2html(thing: dict, uid: str) -> str:
    """Convert webpanel() output to HTML."""
    if thing['type'] == 'figure':
        filename = thing['filename']
        html = f'<img src="/png/{uid}/{filename}" />'
    elif thing['type'] == 'table':
        html = table(thing['header'], thing['rows'])
    else:
        raise ValueError
    return html
