"""Base web-app class."""
from __future__ import annotations

from functools import partial
from io import BytesIO, StringIO
from pathlib import Path

from ase.io import write
from bottle import TEMPLATE_PATH, Bottle, request, static_file, template

from camdweb.html import FormPart, Select, StoichiometryInput
from camdweb.materials import Materials
from camdweb.session import Sessions
from camdweb.panels.panel import SkipPanel

TEMPLATE_PATH[:] = [str(Path(__file__).parent)]


class CAMDApp:
    title = 'CAMD'
    logo = ''
    links = [
        ('CMR', 'https://cmr.fysik.dtu.dk')]

    def __init__(self,
                 materials: Materials,
                 initial_columns: list[str],
                 *,
                 initial_filter_string: str = '',
                 root: Path | None = None):
        self.materials = materials
        self.root = root or Path()

        self.route()

        # For updating plots:
        self.callbacks = self.materials.get_callbacks()

        # User sessions (selected columns, sorting, filter string, ...)
        self.sessions = Sessions(initial_columns,
                                 filter_string=initial_filter_string)

        self.form_parts: list[FormPart] = []

        # For selecting materials (A, AB, AB2, ...):
        stoichiometries = self.materials.stoichiometries()
        if len(stoichiometries) > 20:
            # too many to list them all
            self.form_parts.append(StoichiometryInput())
        else:
            self.form_parts.append(
                Select('Stoichiometry:', 'stoichiometry',
                       [''] + stoichiometries))

        # For nspecies selection:
        maxnspecies = max(len(material.count) for material in self.materials)
        self.form_parts.append(
            Select('Number of chemical species:', 'nspecies',
                   [''] + [str(i) for i in range(1, maxnspecies + 1)]))

    def route(self):
        self.app = Bottle()
        self.app.route('/')(self.index_page)
        self.app.route('/table')(self.table_html)
        self.app.route('/material/<uid>')(self.material_page)
        self.app.route('/callback')(self.callback)
        self.app.route('/png/<path:path>')(self.png)
        self.app.route('/favicon.ico')(self.favicon)

        for fmt in ['xyz', 'cif', 'json']:
            self.app.route(f'/material/<uid>/download/{fmt}')(
                partial(self.download, fmt=fmt))

    def download(self, uid: str, fmt: str) -> bytes | str:
        ase_fmt = fmt

        if fmt == 'xyz':
            # Only the extxyz writer includes cell, pbc etc.
            ase_fmt = 'extxyz'

        # (Can also query ASE's IOFormat for whether bytes or str,
        # in fact, ASE should make this easier.)
        buf: BytesIO | StringIO = BytesIO() if fmt == 'cif' else StringIO()

        atoms = self.materials[uid].atoms
        write(buf, atoms, format=ase_fmt)
        return buf.getvalue()

    def index_page(self) -> str:
        """Page showing table of selected materials."""
        query = request.query
        session = self.sessions.get(int(query.get('sid', '-1')))
        search = '\n'.join(fp.render() for fp in self.form_parts)
        table = self.table_html(session)
        sidebar = self.persistent_sidebar()
        return template('index.html',
                        title=self.title,
                        search=search,
                        session=session,
                        table=table,
                        sidebar=sidebar)

    def table_html(self, session=None) -> str:
        """Get HTML for table."""
        if session is None:
            query = request.query
            session = self.sessions.get(int(query.get('sid')))
            if 'filter' in query:
                filter_string = self.get_filter_string(query)
                session.update(filter=filter_string)
            else:
                session.update(query=query)
        rows, header, pages, new_columns, error = self.materials.get_rows(
            session)
        summary_string = pages.summary(len(self.materials))
        return template('table.html',
                        session=session,
                        pages=pages,
                        rows=rows,
                        summary_string=summary_string,
                        header=header,
                        new_columns=new_columns,
                        error=error)

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
        for s in self.form_parts:
            filters += s.get_filter_strings(query)
        return ','.join(filters)

    def material_page(self, uid: str) -> str:
        """Page showing one selected material."""
        sidebar = self.persistent_sidebar()
        material = self.materials[uid]
        webpanels = []
        for panel in self.materials.panels:
            if not all((material.folder / datafile).is_file()
                       for datafile in panel.datafiles):
                continue
            try:
                data = panel.get_data(material)
            except SkipPanel:
                continue
            webpanels.append(data)

        return template('material.html',
                        title=uid,
                        panels=webpanels,
                        sidebar=sidebar)

    def persistent_sidebar(self):
        """Provide persistent sidebar for all pages."""
        return template('sidebar.html',
                        logo=self.logo,
                        sidebar_links=self.links)

    def callback(self) -> str:
        """Send new json data.

        Currently only used for the atoms-plot.
        """
        query = request.query
        name = query['name']
        uid = query['uid']
        material = self.materials[uid]
        return self.callbacks[name](material, int(query['data']))

    def png(self, path: str) -> bytes:
        """Return binary data for png-figures."""
        return static_file(path, self.root)

    def favicon(self) -> bytes:
        path = self.root / 'favicon.ico'
        return static_file(path.name, path.parent)
