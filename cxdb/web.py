from __future__ import annotations

import sys
from pathlib import Path

from bottle import Bottle, request, template, TEMPLATE_PATH

from cxdb.atoms import AtomsSection
from cxdb.bader import BaderSection
from cxdb.dos import DOSSection
from cxdb.material import Material

TEMPLATE_PATH[:] = [str(Path(__file__).parent)]
# STATIC_PATH ???


class C2DB:
    def __init__(self, materials: dict[str, Material]):
        self.materials = materials

        self.app = Bottle()
        self.app.route('/')(self.index)
        self.app.route('/material/<id>')(self.material)
        self.app.route('/callback')(self.callback)
        self.app.route('/png/<id>/<filename>')(self.png)
        self.app.route('/stop/<code:int>')(self.stop)

        self.sections = [AtomsSection(),
                         DOSSection(),
                         BaderSection()]
        self.callbacks = {}
        for section in self.sections:
            self.callbacks.update(section.callbacks)

    def index(self) -> str:
        query = request.query.get('query', '')
        return template('index.html',
                        query=query,
                        rows=[(id, material.columns)
                              for id, material in self.materials.items()
                              if query in id],
                        header=Material.header)

    def material(self, id: str) -> str:
        material = self.materials[id]
        sections = []
        footer = ''
        for section in self.sections:
            html1, html2 = section.get_html(material)
            if html1:
                sections.append((section.title, html1))
                footer += html2
        return template('material.html',
                        title=id,
                        sections=sections,
                        footer=footer)

    def stop(self, code: int) -> str:
        if code == 117:
            sys.stderr.close()
        return ''

    def callback(self) -> str:
        name = request.query.name
        id = request.query.id
        material = self.materials[id]
        return self.callbacks[name](material, int(request.query.data))

    def png(self, id: str, filename: str) -> bytes:
        material = self.materials[id]
        return (material.folder / filename).read_bytes()


def main() -> None:
    materials = {}
    for arg in sys.argv[1:]:
        folder = Path(arg)
        id = folder.name
        materials[id] = Material(folder, id)

    C2DB(materials).app.run(host='0.0.0.0', port=8081, debug=True)


if __name__ == '__main__':
    main()
