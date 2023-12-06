from __future__ import annotations

import sys
from pathlib import Path

from bottle import Bottle, request, template, TEMPLATE_PATH, static_file

from cxdb.atoms import AtomsSection
from cxdb.bader import BaderSection
from cxdb.dos import DOSSection
from cxdb.material import Material, Materials
from cxdb.session import Sessions
from cxdb.panel import Panel

TEMPLATE_PATH[:] = [str(Path(__file__).parent)]


class C2DB:
    def __init__(self,
                 materials: Materials,
                 panels: list[Panel],
                 column_names: dict[str, str],
                 initial_columns: set[str],
                 root: Path):
        self.materials = materials
        self.panels = panels
        self.column_names = column_names
        self.root = root

        self.app = Bottle()
        self.app.route('/')(self.index)
        self.app.route('/material/<id>')(self.material)
        self.app.route('/callback')(self.callback)
        self.app.route('/png/<id>/<filename>')(self.png)
        self.app.route('/help')(self.help)

        self.callbacks = {}
        for panel in self.panels:
            self.callbacks.update(panel.callbacks)

        # User sessions (selcted columns, sorting, filter string, ...)
        self.sessions = Sessions(initial_columns)

        # For selecting materials (A, AB, AB2, ...)
        self.stoichiometries = {'Any'} + self.materials.stoichiometries()

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
        panels = []
        footer = ''
        for panel in self.panels:
            html1, html2 = panel.get_html(material)
            if html1:
                panels.append((panel.title, html1))
                footer += html2
        return template('material.html',
                        title=id,
                        panels=panels,
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
        return static_file(material.folder / filename, self.root)


def main() -> None:
    panels = [AtomsSection(),
              DOSSection(),
              BaderSection()]

    materials = Materials(panels)

    for arg in sys.argv[1:]:
        folder = Path(arg)
        id = folder.name
        material = Material(folder, id)
        materials.add(material)

        for panel in panels:
            panel.update_column_data(material)

    initial_columns = {'id', 'energy', 'volume', 'formula'}

    root = Path.cwd()

    C2DB(materials, panels, initial_columns, root).app.run(
        host='0.0.0.0', port=8081, debug=True)


if __name__ == '__main__':
    main()
