import importlib
from pathlib import Path

from ase.io.jsonio import decode

from cxdb.material import Material, Materials
from cxdb.panel import Panel
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


class Row:
    def __init__(self, material):
        self.data = Data(material.folder)
        self.magstate = material.magstate
        self.has_inversion_symmetry = material.has_inversion_symmetry
        self.cell = material.atoms.cell
        self.minhessianeig = 117.0

    def get(self, name, default=None):
        if hasattr(self, name):
            return getattr(self, name)
        print('MISSING:', name, default)
        return default


class Data:
    def __init__(self, folder):
        self.folder = folder

    def get(self, name, default=None):
        dct = decode((self.folder / name).read_text())
        if 'kwargs' in dct:
            return dct['kwargs']['data']
        return dct

    def __getitem__(self, name):
        return self.get(name)


class ASRPanel(Panel):
    def __init__(self, name, index=None):
        mod = importlib.import_module(f'asr.{name}')
        self.webpanel = mod.webpanel
        self.key_descriptions = mod.Result.key_descriptions

    def get_html(self,
                 material: Material,
                 materials: Materials) -> tuple[str, str]:
        uid = material.uid
        row = Row(material)
        (p,) = self.webpanel(None, row, self.key_descriptions)

        columns: list[list[str]] = []
        for column in p['columns']:
            columns.append([])
            for thing in column:
                html = thing2html(thing, uid)
                columns[-1].append(html)

        for desc in p['plot_descriptions']:
            paths = [material.folder / filename
                     for filename in desc['filenames']]
            for f in paths:
                if not f.is_file():
                    desc['function'](row, *paths)
                    # for filename in desc['filenames']:
                    #     f0 = Path(filename)
                    #     f.write_bytes(f0.read_bytes())
                    #     f0.unlink()
                    break

        self.title = p['title']
        return (HTML.format(title=p['title'],
                            col1='\n'.join(columns[0]),
                            col2='\n'.join(columns[1])),
                '')


def thing2html(thing: dict, uid: str) -> str:
    if thing['type'] == 'figure':
        filename = thing['filename']
        html = f'<img src="/png/{uid}/{filename}" />'
    elif thing['type'] == 'table':
        html = table(thing['header'], thing['rows'])
    else:
        print(thing)
        1 / 0
    return html
