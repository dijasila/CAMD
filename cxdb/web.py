from __future__ import annotations

import sys
from pathlib import Path

from bottle import Bottle, request, template, TEMPLATE_PATH, static_file

from cxdb.atoms import AtomsPanel
from cxdb.bader import BaderPanel
from cxdb.dos import DOSPanel
from cxdb.material import Material, Materials
from cxdb.session import Sessions

TEMPLATE_PATH[:] = [str(Path(__file__).parent)]


class CXDBApp:
    def __init__(self,
                 materials: Materials,
                 initial_columns: set[str],
                 root: Path):
        self.materials = materials
        self.root = root

        self.app = Bottle()
        self.app.route('/')(self.index)
        self.app.route('/material/<uid>')(self.material)
        self.app.route('/callback')(self.callback)
        self.app.route('/png/<uid>/<filename>')(self.png)
        self.app.route('/help')(self.help)

        self.callbacks = self.materials.get_callbacks()

        # User sessions (selcted columns, sorting, filter string, ...)
        self.sessions = Sessions(initial_columns)

        # For selecting materials (Any, A, AB, AB2, ...)
        self.stoichiometries = [
            ('Any', 'Any')] + self.materials.stoichiometries()

    def index(self,
              query: dict | None = None) -> str:
        if query is None:
            query = request.query

        session = self.sessions.get(int(query.get('sid', '-1')))
        session.update(query)

        rows, header, pages, new_columns = self.materials.get_rows(session)
        return template('index.html',
                        query=query,
                        stoichiometries=self.stoichiometries,
                        session=session,
                        pages=pages,
                        rows=rows,
                        header=header,
                        new_columns=new_columns)

    def material(self, uid: str) -> str:
        if uid == 'stop':
            sys.stderr.close()
        material = self.materials[uid]
        panels = []
        footer = ''
        for panel in self.materials.panels:
            html1, html2 = panel.get_html(material, self.materials)
            if html1:
                panels.append((panel.title, html1))
                footer += html2
        return template('material.html',
                        title=uid,
                        panels=panels,
                        footer=footer)

    def callback(self, query: dict | None = None) -> str:
        if query is None:
            query = request.query
        name = query['name']
        uid = query['uid']
        material = self.materials[uid]
        return self.callbacks[name](material, int(query['data']))

    def help(self):
        return template('help.html')

    def png(self, uid: str, filename: str) -> bytes:
        material = self.materials[uid]
        return static_file(str(material.folder / filename), self.root)


def main() -> None:
    panels = [AtomsPanel(3),
              DOSPanel(),
              BaderPanel()]

    mlist = []
    for arg in sys.argv[1:]:
        folder = Path(arg)
        uid = folder.name
        mlist.append(Material.from_file(folder / 'rlx.traj', uid))

    materials = Materials(mlist, panels)

    initial_columns = {'uid', 'energy', 'volume', 'formula'}

    root = Path.cwd()

    CXDBApp(materials, initial_columns, root).app.run(
        host='0.0.0.0', port=8081, debug=True)


if __name__ == '__main__':
    main()
