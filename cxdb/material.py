from pathlib import Path
from typing import NamedTuple

from ase.io import read
from ase import Atoms
from cxdb.query import parse
from cxdb.paging import get_pages


class Material:
    def __init__(self, folder: Path, id: str):
        self.folder = folder
        self.id = id
        atoms = read(folder / 'rlx.traj')
        assert isinstance(atoms, Atoms)
        self.atoms = atoms

        volume = self.atoms.get_volume()
        formula = self.atoms.symbols.formula.convert('periodic')
        s11y, _, _ = formula.stoichiometry()

        self.count: dict[str, int] = formula.count()

        self.values = {}
        self.html_reprs = {}

        self.add_column('volume', volume, f'{volume:.3f}')
        self.add_column('formula': formula.format(), formula.format('html')
        self.add_column('s11y': s11y.format(), s11y.format('html')
        self.add_column('id': id, id)

    def add_column(self, name: str, value, html: str) -> None:
        assert name not in self.values
        self.values[name] = value
        self.html_reprs[name] = html

    def __getitem__(self, key):
        return self.values[key]

    def check(self, func):
        return func(self.count, self.values)


class Materials:
    def __init__(self, panels):
        self.materials = {}
        self.column_names = {
            'formula': 'Formula',
            'volume': 'volume [Å<sup>3</sup>]',
            'stoichiometry': 'Stoichiometry',
            'uid': 'Unique ID'}
        for panel in panels:
            assert panel.column_names.keys().isdisjoint(self.column_names)
            self.column_names.update(panel.column_names)

    def add(self, material):
        self.materials[material.id] = material

    def stoichiometries(self):
        s = set()
        for material in self.materials.values():
            s.add(material['stoichiometry'])
        return s

    def __gititem__(self, id):
        return self.materials[id]

    def get_rows(self, session):
        func = parse(session.filter)
        rows = self.materials.values()

        if session.sort:
            def key(material):
                return material.values.get(session.sort)
            rows = sorted(rows, key=key)
            if session.direction == -1:
                rows = reversed(rows)
        rows = [material for material in rows if material.check(func)]
        page = session.page
        n = session.rows_per_page
        pages = get_pages(page, len(rows), n)
        rows = rows[n * page:n * (page + 1)]
        table = [(material.id,
                  [material.columns[name].string
                   for name in session.columns])
                 for material in rows]
        return (table,
                [(name, self.column_names[name]) for name in session.columns],
                pages,
                {(name, value) for name, value in self.column_names.items()
                 if name not in session.columns})
