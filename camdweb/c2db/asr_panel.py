"""Hack to use webpanel() functions from ASR."""
import gzip
import importlib
from pathlib import Path
from typing import Generator

from ase.db.core import KeyDescription
from ase.io.jsonio import decode
from camdweb.html import table
from camdweb.material import Material, Materials
from camdweb.panels.panel import Panel

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
    gz = path.with_suffix('.json.gz')
    if gz.is_file():
        with gzip.open(gz, 'rt') as fd:
            txt = fd.read()
    else:
        txt = path.read_text()
    dct = decode(txt)
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
        try:
            dct = read_result_file(self.folder / name)
        except FileNotFoundError:
            return None
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
                 materials: Materials) -> (
            tuple[str, str] |
            Generator[tuple[str, str], None, None]):
        """Create row and result objects and call webpanel() function."""
        row = Row(material)
        dct = row.data.get(f'results-asr.{self.name}.json')
        if dct is None:
            return
        result = self.result_class(dct)
        (p,) = self.webpanel(result, row, self.key_descriptions)
        self.title = p['title']

        columns: list[list[str]] = [[], []]
        for i, column in enumerate(p['columns']):
            for thing in column:
                if thing is not None:  # pragma: no branch
                    html = thing2html(thing, material.folder)
                    columns[i].append(html)

        async_results = []
        for desc in p.get('plot_descriptions', []):
            if desc['filenames'] == ['bs.html']:
                continue
            paths = [material.folder / filename
                     for filename in desc['filenames']]
            for f in paths:
                if not f.is_file():
                    # Call plot-function:
                    if 1:
                        pool = materials.process_pool
                        result = pool.apply_async(
                            desc['function'], (row, *paths))
                        async_results.append(result)
                    else:
                        desc['function'](row, *paths)
                    break

        yield

        for result in async_results:
            result.wait()

        html = HTML.format(title=p['title'],
                           col1='\n'.join(columns[0]),
                           col2='\n'.join(columns[1]))

        yield html


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
