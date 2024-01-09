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
                 initial_columns: list[str],
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
        self.stoichiometries = ['Any'] + self.materials.stoichiometries()
        self.maxnspecies = max(material.nspecies
                               for material
                               in self.materials._materials.values())

    def index(self,
              query: dict | None = None) -> str:
        if query is None:
            query = request.query

        session = self.sessions.get(int(query.get('sid', '-1')))
        session.update(self.get_filter_string(query), query)

        rows, header, pages, new_columns = self.materials.get_rows(session)
        return template('index.html',
                        query=query,
                        stoichiometries=self.stoichiometries,
                        maxnspecies=self.maxnspecies,
                        session=session,
                        pages=pages,
                        rows=rows,
                        header=header,
                        new_columns=new_columns)

    def get_filter_string(self, query: dict) -> str:
        """Generate filter string from URL query.

        Example::

            {'filter': Cu=1,gap>1.5',
             'stoichiometry': 'AB2',
             'nspecies': ''}

        will give the string "Cu=1,gap>1.5,stoichiometry=AB2".
        """
        filters = []
        filter = query.get('filter', '')
        if filter:
            filters.append(filter)
        s11y = query.get('stoichiometry', 'Any')
        if s11y != 'Any':
            filters.append(f'stoichiometry={s11y}')
        nspecies = query.get('nspecies', '')
        if nspecies:
            filters.append(f'nspecies={nspecies}')
        return ','.join(filters)

    def material(self, uid: str) -> str:
        if uid == 'stop':  # pragma: no cover
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
            query = request.query  # pragma: no cover
        name = query['name']
        uid = query['uid']
        material = self.materials[uid]
        return self.callbacks[name](material, int(query['data']))

    def help(self):
        return template('help.html')

    def png(self, uid: str, filename: str) -> bytes:
        material = self.materials[uid]
        return static_file(str(material.folder / filename), self.root)


def main() -> None:  # pragma: no cover
    panels = [AtomsPanel(3),
              DOSPanel(),
              BaderPanel()]

    mlist = []
    for arg in sys.argv[1:]:
        folder = Path(arg)
        uid = folder.name
        mlist.append(Material.from_file(folder / 'rlx.traj', uid))

    materials = Materials(mlist, panels)

    initial_columns = ['uid', 'volume', 'formula']

    root = Path.cwd()

    CXDBApp(materials, initial_columns, root).app.run(
        host='0.0.0.0', port=8081, debug=True)


if __name__ == '__main__':
    main()
