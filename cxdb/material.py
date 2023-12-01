from pathlib import Path
from typing import NamedTuple

from ase.io import read
from ase import Atoms
from cxdb.query import parse
from cxdb.paging import get_pages


class Column(NamedTuple):
    value: bool | int | float | str
    string: str


class Material:
    headers = {'formula': 'Formula',
               'energy': 'energy [eV]',
               'volume': 'volume [Å<sup>3</sup>]',
               'id': 'Unique ID'}

    def __init__(self, folder: Path, id: str):
        self.folder = folder
        self.id = id
        atoms = read(folder / 'rlx.traj')
        assert isinstance(atoms, Atoms)
        self.atoms = atoms

        energy = self.atoms.get_potential_energy()
        volume = self.atoms.get_volume()

        formula = self.atoms.symbols.formula.convert('periodic')
        s11y, _, _ = formula.stoichiometry()

        self.columns = {
            'energy': Column(energy, f'{energy:.3f}'),
            'volume': Column(volume, f'{volume:.3f}'),
            'formula': Column(formula.format(), formula.format('html')),
            's11y': Column(s11y.format(), s11y.format('html')),
            'id': Column(id, id)}

        self.values = {key: column.value
                       for key, column in self.columns.items()}
        self.count = formula.count()

    def __getitem__(self, key):
        return self.values[key]

    def check(self, func):
        return func(self.count, self.values)


def get_rows(materials, session):
    func = parse(session.filter)
    rows = materials.values()

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
            [(name, Material.headers[name]) for name in session.columns],
            pages,
            {(name, value) for name, value in Material.headers.items()
             if name not in session.columns})
