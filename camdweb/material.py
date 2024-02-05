from __future__ import annotations

import multiprocessing
from functools import cached_property
from math import nan
from pathlib import Path
from typing import Generator, Sequence

import numpy as np
from ase import Atoms
from ase.io import read

from camdweb.filter import Index, parse
from camdweb.paging import get_pages
from camdweb.panels.panel import Panel
from camdweb.session import Session
from camdweb.utils import fft, html_format_formula


class Material:
    def __init__(self, folder: Path, uid: str, atoms: Atoms):
        """Object representing a material and associated data.

        >>> mat = Material(Path(), 'x1', Atoms('H2O'))
        >>> mat.formula
        'OH2'
        >>> mat['formula']
        'OH<sub>2</sub>'
        >>> mat.stoichiometry
        'AB2'
        >>> mat.add_column('energy', -1.23456)
        >>> mat.energy, mat['energy']
        (-1.23456, '-1.235')
        """
        self.folder = folder
        self.uid = uid
        self.atoms = atoms

        # Get number-of-atoms dicts:
        self.count, self.formula, self.reduced, self.stoichiometry = fft(
            atoms.numbers)

    def get_columns(self):
        cols = {'uid': self.uid,
                'natoms': sum(self.count.values()),
                'nspecies': len(self.count),
                'formula': self.formula,
                'reduced': self.reduced,
                'stoichiometry': self.stoichiometry}
        pbc = self.atoms.pbc
        dims = pbc.sum()
        if dims > 0:
            vol = abs(np.linalg.det(self.atoms.cell[pbc][:, pbc]))
            name = ['length', 'area', 'volume'][dims - 1]
            cols[name] = vol
        return cols

    def html_format_column(self,
                           key: str,
                           value: bool | int | float | str) -> str:
        if isinstance(value, float):
            return f'{value:.3f}'
        if key in ['formula', 'reduced', 'stoichiometry']:
            return html_format_formula(value)
        return str(value)

    @classmethod
    def from_file(cls, file: Path, uid: str) -> Material:
        atoms = read(file)
        assert isinstance(atoms, Atoms)
        return cls(file.parent, uid, atoms)


class Materials:
    def __init__(self,
                 materials: list[Material],
                 panels: Sequence[Panel],
                 column_names: dict[str, str]):
        self.column_names = {
            'formula': 'Formula',
            'reduced': 'Reduced formula',
            'stoichiometry': 'Stoichiometry',
            'nspecies': 'Number of species',
            'natoms': 'Number of atoms',
            'uid': 'Unique ID',
            'length': 'Unit cell length [Å]',
            'area': 'Unit cell area [Å<sup>2</sup>]',
            'volume': 'Unit cell volume [Å<sup>3</sup>]'}
        self.column_names.update(column_names)

        self._materials: dict[str, Material] = {}
        for material in materials:
            for panel in panels:
                panel.update_data(material)
            self._materials[material.uid] = material

        for panel in panels:
            if not panel.column_names.keys().isdisjoint(self.column_names):
                overlap = panel.column_names.keys() & self.column_names
                raise ValueError(f'{overlap}')
            self.column_names.update(panel.column_names)

        self.index = Index([(mat.reduced, mat.count, mat.get_columns())
                            for mat in self._materials.values()])
        self.i2uid = {i: mat.uid for i, mat in enumerate(self)}

        self.panels = panels

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
        return list(s)

    def table(self,
              material: Material,
              columns: list[str]) -> list[tuple[str, str]]:
        return [(self.column_names[name], material[name])
                for name in columns
                if name in material.columns]

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
                [(name, self.column_names[name]) for name in session.columns],
                pages,
                [(name, value) for name, value in self.column_names.items()
                 if name not in session.columns],
                error)

    @cached_property
    def process_pool(self):
        return multiprocessing.Pool(maxtasksperchild=100)
