from __future__ import annotations

from math import nan
from typing import Generator, Sequence

from camdweb.filter import Index
from camdweb.parse import parse
from camdweb.material import Material
from camdweb.paging import get_pages
from camdweb.panels.panel import Panel
from camdweb.session import Session
from camdweb.utils import html_format_formula


class Materials:
    def __init__(self,
                 materials: Sequence[Material],
                 panels: Sequence[Panel]):
        self.panels = panels
        self._materials = {material.uid: material for material in materials}
        for material in self:
            for panel in panels:
                panel.update_material(material)

        self.index = Index([(mat.count, mat.columns) for mat in self])
        self.i2uid = {i: mat.uid for i, mat in enumerate(self)}
        self.column_descriptions = {
            'formula': 'Formula',
            'reduced': 'Reduced formula',
            'stoichiometry': 'Stoichiometry',
            'nspecies': 'Number of species',
            'natoms': 'Number of atoms',
            'uid': 'Unique ID',
            'length': 'Unit cell length [Å]',
            'area': 'Unit cell area [Å<sup>2</sup>]',
            'volume': 'Unit cell volume [Å<sup>3</sup>]'}

        self.html_formaters = {
            name: html_format_formula
            for name in ['formula', 'reduced', 'stoichiometry']}

        self._update_panels()

    def _update_panels(self):
        for panel in self.panels:
            panel.update_column_descriptions(self.column_descriptions)
            panel.update_html_formatters(self.html_formatters)

        keys: set[str] = set()
        for material in self:
            keys.update(material.columns)

        for name in self.column_descriptions - keys:
            del self.column_descriptions[name]

    def __iter__(self) -> Generator[Material, None, None]:
        yield from self._materials.values()

    def __len__(self) -> int:
        return len(self._materials)

    def get_callbacks(self):
        callbacks = {}
        for panel in self.panels:
            callbacks.update(panel.callbacks)
        return callbacks

    def stoichiometries(self) -> list[str]:
        """Construct list of stoichiometries present."""
        s = set()
        for material in self:
            s.add(material.stoichiometry)
        return sorted(s)

    def __getitem__(self, uid: str) -> Material:
        return self._materials[uid]

    def get_rows(self,
                 session: Session) -> tuple[list[tuple[str, list[str]]],
                                            list[tuple[str, str]],
                                            list[tuple[int, str]],
                                            list[tuple[str, str]],
                                            str]:
        """Filter rows for table.

        Example::

            rows, header, pages, new_columns = materials.get_rows(session)

        The returned values are:

        rows:
            list of rows, where each row is a tuple of uid and list of
            HTML-strings.

        header:
            list of (column name, column HTML-string) tuples.

        pages:
            stuff for pagination buttons (see get_pages() function).

        new_columns:
            list of (column name, columns HTML-string) tuples for columns not
            shown.

        error:
            Error message.
        """
        filter = session.filter
        try:
            func = parse(filter)
        except SyntaxError as ex:
            error = ex.args[0]
            rows = []
        else:
            rows = list(func(self.index))
            error = ''

        if rows and session.sort:
            missing = '' if session.sort in self.index.strings else nan

            def key(i):
                return self.index.columns[i].get(session.sort, missing)

            rows = sorted(rows, key=key, reverse=session.direction == -1)

        page = session.page
        n = session.rows_per_page
        pages = get_pages(page, len(rows), n)
        rows = rows[n * page:n * (page + 1)]
        table = []
        for i in rows:
            uid = self.i2uid[i]
            material = self._materials[uid]
            columns = self.index.columns[i]
            table.append(
                (uid,
                 [material.html_format_column(name, columns.get(name, ''))
                  for name in session.columns]))
        return (table,
                [(name, self.column_descriptions.get(name, name))
                 for name in session.columns],
                pages,
                [(name, value)
                 for name, value in self.column_descriptions.items()
                 if name not in session.columns],
                error)
