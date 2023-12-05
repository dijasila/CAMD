from __future__ import annotations

import sys
from pathlib import Path

from bottle import Bottle, request, template, TEMPLATE_PATH

from cxdb.atoms import AtomsSection
from cxdb.bader import BaderSection
from cxdb.dos import DOSSection
from cxdb.material import Material, get_rows
from cxdb.session import Sessions

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
        self.app.route('/help')(self.help)

        self.sections = [AtomsSection(),
                         DOSSection(),
                         BaderSection()]
        self.callbacks = {}
        for section in self.sections:
            self.callbacks.update(section.callbacks)

        self.sessions = Sessions({'id', 'energy', 'volume', 'formula'})

        self.stoichiometries = {'Any'}
        for material in materials.values():
            self.stoichiometries.add(material['s11y'])

    def index(self,
              query: dict | None = None) -> str:
        if query is None:
            query = request.query

        session = self.sessions.get(int(query.get('sid', '-1')))
        session.update(query)

        rows, header, pages, new_columns = get_rows(self.materials, session)
        return template('index.html',
                        query=query,
                        stoichiometries=self.stoichiometries,
                        session=session,
                        pages=pages,
                        rows=rows,
                        header=header,
                        new_columns=new_columns)

    def material(self, id: str) -> str:
        if id == 'stop':
            sys.stderr.close()
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

    def callback(self, query: dict | None = None) -> str:
        if query is None:
            query = request.query
        name = query['name']
        id = query['id']
        material = self.materials[id]
        return self.callbacks[name](material, int(query['data']))

    def help(self):
        return template('help.html')

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
