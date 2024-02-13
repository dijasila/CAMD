"""Base web-app class."""
from __future__ import annotations

import re
from functools import partial
from io import BytesIO, StringIO
from pathlib import Path
from typing import Iterator

from ase.io import write
from bottle import TEMPLATE_PATH, Bottle, request, static_file, template

from camdweb.html import FormPart, Select
from camdweb.materials import Materials
from camdweb.session import Sessions

TEMPLATE_PATH[:] = [str(Path(__file__).parent)]


class CAMDApp:
    title = 'CAMD'

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

        self.form_parts: list[FormPart] = []

        # For selecting materials (A, AB, AB2, ...):
        self.form_parts.append(
            Select('Stoichiometry', 'stoichiometry',
                   [''] + self.materials.stoichiometries()))

        # For nspecies selection:
        maxnspecies = max(len(material.count) for material in self.materials)
        self.form_parts.append(
            Select('Number of chemical species', 'nspecies',
                   [''] + [str(i) for i in range(1, maxnspecies + 1)]))

    def route(self):
        self.app = Bottle()
        self.app.route('/')(self.index_page)
        self.app.route('/table')(self.table_html)
        self.app.route('/material/<uid>')(self.material_page)
        self.app.route('/callback')(self.callback)
        self.app.route('/png/<path:path>')(self.png)

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
        table = self.get_table(session)
        return template('index.html',
                        title=self.title,
                        search=search,
                        session=session,
                        table=table)

    def get_table(self, session) -> tuple[str, str]:
        """Get HTML and error message (if there was one) for table."""
        rows, header, pages, new_columns, error = self.materials.get_rows(
            session)
        return template('table.html',
                        session=session,
                        pages=pages,
                        rows=rows,
                        header=header,
                        new_columns=new_columns,
                        error=error)

    def table_html(self) -> str:
        """Get HTML for table."""
        query = request.query
        filter_string = self.get_filter_string(query)
        session = self.sessions.get(int(query.get('sid', '-1')))
        session.update(filter_string, query)
        return self.get_table(session)

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
        material = self.materials[uid]
        titles = []
        info_strings = []
        generators: list[Iterator[str]] = []
        for panel in self.materials.panels:
            if not all((material.folder / datafile).is_file()
                       for datafile in panel.datafiles):
                continue
            generator = panel.get_html(material)
            try:
                html = next(generator)
            except StopIteration:
                continue
            if html == '':  # result will come next time
                generators.append(generator)
            else:
                generators.append(iter([html]))
            titles.append(panel.title)
            info_strings.append(panel.info)

        panels = []
        scripts = []
        for gen, title, info in zip(generators, titles, info_strings):
            html = next(gen)
            html, script = cut_out_script(html)
            panels.append((title, info, html))
            scripts.append(script)

        return template('material.html',
                        title=uid,
                        panels=panels,
                        footer='\n'.join(scripts))

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


def cut_out_script(html: str) -> tuple[str, str]:
    r"""We need to put the script tags in the footer.

    >>> cut_out_script('''Hello
    ... <script>
    ...   ...
    ... </script>
    ... CAMd''')
    ('Hello\n\nCAMd', '<script>\n  ...\n</script>')
    """
    m = re.search(r'(<script.*</script>)', html, re.MULTILINE | re.DOTALL)
    if m:
        i, j = m.span()
        return html[:i] + html[j:], html[i:j]
    return html, ''
