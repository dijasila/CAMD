from __future__ import annotations

from math import nan
from typing import Callable, Generator, Sequence

from camdweb.filter import Index
from camdweb.material import Material
from camdweb.paging import get_pages
from camdweb.panels.panel import Panel, default_formatter
from camdweb.parse import parse
from camdweb.session import Session
from camdweb.utils import html_format_formula
from camdweb import ColVal


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

        self.html_formatters: dict[str, Callable[[ColVal, bool], str]] = {
            name: html_format_formula
            for name in ['formula', 'reduced', 'stoichiometry']}

        self._update_panels()

    def _update_panels(self) -> None:
        for panel in self.panels:
            panel.update_column_descriptions(self.column_descriptions)
            panel.update_html_formatters(self.html_formatters)

        keys: set[str] = set()
        for material in self:
            keys.update(material.columns)

        for name in self.column_descriptions.keys() - keys:
            del self.column_descriptions[name]

    def __iter__(self) -> Generator[Material, None, None]:
        yield from self._materials.values()

    def __len__(self) -> int:
        return len(self._materials)

    def get_callbacks(self) -> dict[str, Callable[[Material, int], str]]:
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
            col_numbers = []
        else:
            col_numbers = list(func(self.index))
            error = ''

        rows = [self[self.i2uid[i]] for i in col_numbers]

        if rows and session.sort:
            missing = '' if session.sort in self.index.strings else nan

            def key(material: Material) -> ColVal:
                return material.columns.get(session.sort, missing)

            rows = sorted(rows, key=key, reverse=session.direction == -1)

        page = session.page
        n = session.rows_per_page
        pages = get_pages(page, len(rows), n)
        rows = rows[n * page:n * (page + 1)]

        formatters = [self.html_formatters.get(name, default_formatter)
                      for name in session.columns]
        table = []
        for material in rows:
            columns = []
            for name, formatter in zip(session.columns, formatters):
                value = material.columns.get(name)
                if value is None:
                    columns.append('')
                else:
                    columns.append(formatter(value))
            table.append((material.uid, columns))

        headers = [(name, self.column_descriptions.get(name, name))
                   for name in session.columns]

        new_columns = [(name, value)
                       for name, value in self.column_descriptions.items()
                       if name not in session.columns]

        return (table, headers, pages, new_columns, error)
