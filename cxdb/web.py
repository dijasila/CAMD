"""Base web-app class."""
from __future__ import annotations

import sys
from pathlib import Path

from bottle import Bottle, request, template, TEMPLATE_PATH, static_file

from cxdb.material import Materials
from cxdb.session import Sessions
from cxdb.utils import Select

TEMPLATE_PATH[:] = [str(Path(__file__).parent)]


class CXDBApp:
    title = 'CXDB'

    def __init__(self,
                 materials: Materials,
                 initial_columns: list[str],
                 root: Path | None = None):
        self.materials = materials
        self.root = root or Path()

        self.route()

        # For updating plots:
        self.callbacks = self.materials.get_callbacks()

        # User sessions (selected columns, sorting, filter string, ...)
        self.sessions = Sessions(initial_columns)

        self.search = []

        # For selecting materials (A, AB, AB2, ...):
        self.search.append(
            Select('Stoichiometry', 'stoichiometry',
                   [''] + self.materials.stoichiometries()))

        # For nspecies selection:
        maxnspecies = max(material.nspecies for material in self.materials)
        self.search.append(
            Select('Number of chemical species', 'nspecies',
                   [''] + [str(i) for i in range(1, maxnspecies)]))

    def route(self):
        self.app = Bottle()
        self.app.route('/')(self.index)
        self.app.route('/material/<uid>')(self.material)
        self.app.route('/callback')(self.callback)
        self.app.route('/png/<uid>/<filename>')(self.png)
        self.app.route('/help')(self.help)

    def index(self,
              query: dict | None = None) -> str:
        """Page showing table of selected materials."""
        if query is None:
            query = request.query

        filter_string = self.get_filter_string(query)
        session = self.sessions.get(int(query.get('sid', '-1')))
        session.update(filter_string, query)
        search = '\n'.join(s.render(query) for s in self.search)
        rows, header, pages, new_columns = self.materials.get_rows(session)

        return template('index.html',
                        title=self.title,
                        query=query,
                        search=search,
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
        for s in self.search:
            filters += s.get_filter_strings(query)
        return ','.join(filters)

    def material(self, uid: str) -> str:
        """Page showing one selected material."""
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
