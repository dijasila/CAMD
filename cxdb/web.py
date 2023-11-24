from __future__ import annotations

import sys
from pathlib import Path

from bottle import Bottle, request, template

from cxdb.atoms import AtomsSection
from cxdb.bader import BaderSection
from cxdb.dos import DOSSection
from cxdb.material import Material

TEMPLATES = Path(__file__).parent


class C2DB:
    def __init__(self, materials: dict[str, Material]):
        self.materials = materials

        self.app = Bottle()
        self.app.route('/')(self.index)
        self.app.route('/material/<id>')(self.material)
        self.app.route('/callback')(self.callback)
        self.app.route('/png/<id>/<filename>')(self.png)

        self.sections = [AtomsSection(),
                         DOSSection(),
                         BaderSection()]
        self.callbacks = {}
        for section in self.sections:
            self.callbacks.update(section.callbacks)

    def index(self) -> str:
        query = request.query.get('query', '')
        print(repr(query))
        print(list(self.materials))
        return template(str(TEMPLATES / 'index.html'),
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
        return template(str(TEMPLATES / 'material.html'),
                        title=id,
                        sections=sections,
                        footer=footer)

    def callback(self) -> str:
        name = request.query.name
        id = request.query.id
        material = self.materials[id]
        return self.callbacks[name](material, int(request.query.data))

    def png(self, id: str, filename: str) -> bytes:
        return (Path(id) / filename).read_bytes()


def main():
    materials = {}
    for arg in sys.argv[1:]:
        folder = Path(arg)
        id = folder.name
        materials[id] = Material(folder, id)

    C2DB(materials).app.run(host='localhost', port=8080, debug=True)


if __name__ == '__main__':
    main()
