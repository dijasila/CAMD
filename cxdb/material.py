from __future__ import annotations

from collections.abc import Container
from math import nan
from pathlib import Path
from typing import Any

from ase import Atoms
from ase.io import read

from cxdb.filter import Index, parse
from cxdb.paging import get_pages
from cxdb.panel import Panel
from cxdb.session import Session


class Material:
    def __init__(self, folder: Path, uid: str, atoms: Atoms):
        self.folder = folder
        self.uid = uid
        self.atoms = atoms

        formula = self.atoms.symbols.formula.convert('periodic')
        s11y, _, _ = formula.stoichiometry()

        self._count: dict[str, int] = formula.count()

        self._data: dict[str, Any] = {}
        self._html_reprs: dict[str, str] = {}
        self._values: dict[str, bool | int | float | str] = {}

        self.add_column('formula', formula.format(), formula.format('html'))
        self.add_column('stoichiometry', s11y.format(), s11y.format('html'))
        self.add_column('uid', uid)

    @classmethod
    def from_file(cls, file: Path, uid: str) -> Material:
        atoms = read(file)
        assert isinstance(atoms, Atoms)
        return cls(file.parent, uid, atoms)

    def add(self,
            name: str,
            value) -> None:
        assert name not in self._data, (name, self._data)
        self._data[name] = value

    def add_column(self,
                   name: str,
                   value,
                   html: str | None = None) -> None:
        self.add(name, value)
        self._values[name] = value
        if html is None:
            if isinstance(value, float):
                html = f'{value:.3f}'
            else:
                html = str(value)
        self._html_reprs[name] = html

    def __getattr__(self, name: str) -> Any:
        return self._data[name]

    def __getitem__(self, name: str) -> str:
        return self._html_reprs[name]

    def get(self, name: str, default: str = '') -> str:
        return self._html_reprs.get(name, default)

    def check_columns(self, column_names: Container[str]) -> None:
        for name in self._values:
            assert name in column_names, name


class Materials:
    def __init__(self,
                 materials: list[Material],
                 panels: list[Panel]):
        self.column_names = {
            'formula': 'Formula',
            'stoichiometry': 'Stoichiometry',
            'uid': 'Unique ID'}

        for panel in panels:
            assert panel.column_names.keys().isdisjoint(self.column_names)
            self.column_names.update(panel.column_names)

        self._materials: dict[str, Material] = {}
        for material in materials:
            for panel in panels:
                panel.update_data(material)
            material.check_columns(self.column_names)
            self._materials[material.uid] = material

        self.index = Index([(mat._count, mat._values)
                            for mat in self._materials.values()])
        self.i2uid = {i: mat.uid
                      for i, mat in enumerate(self._materials.values())}

        self.panels = panels

    def get_callbacks(self):
        callbacks = {}
        for panel in self.panels:
            callbacks.update(panel.callbacks)
        return callbacks

    def stoichiometries(self):
        s = []
        for material in self._materials.values():
            s.append((material.stoichiometry, material['stoichiometry']))
        return s

    def __getitem__(self, uid):
        return self._materials[uid]

    def get_rows(self,
                 session: Session) -> tuple[list[tuple[str, list[str]]],
                                            list[tuple[str, str]],
                                            list[tuple[int, str]],
                                            list[tuple[str, str]]]:
        filter = session.filter
        if session.stoichiometry != 'Any':
            if filter:
                filter += ','
            filter += f'stoichiometry={session.stoichiometry}'
        func = parse(filter)
        rows = [self._materials[self.i2uid[i]] for i in func(self.index)]

        if rows and session.sort:
            missing = '' if session.sort in self.index.strings else nan

            def key(material):
                return material._values.get(session.sort, missing)

            rows = sorted(rows, key=key, reverse=session.direction == -1)

        page = session.page
        n = session.rows_per_page
        pages = get_pages(page, len(rows), n)
        rows = rows[n * page:n * (page + 1)]
        table = [(material.uid,
                  [material.get(name, '') for name in session.columns])
                 for material in rows]
        return (table,
                [(name, self.column_names[name]) for name in session.columns],
                pages,
                [(name, value) for name, value in self.column_names.items()
                 if name not in session.columns])
