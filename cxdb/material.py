from pathlib import Path

from ase.io import read
from ase import Atoms
from cxdb.query import parse
from cxdb.paging import get_pages
from cxdb.session import Session


class Material:
    def __init__(self, folder: Path, uid: str):
        self.folder = folder
        self.uid = uid
        atoms = read(folder / 'rlx.traj')
        assert isinstance(atoms, Atoms)
        self.atoms = atoms

        volume = self.atoms.get_volume()
        formula = self.atoms.symbols.formula.convert('periodic')
        s11y, _, _ = formula.stoichiometry()

        self.count: dict[str, int] = formula.count()

        self.values: dict[str, bool | int | float | str] = {}
        self.html_reprs: dict[str, str] = {}

        self.add_column('volume', volume, f'{volume:.3f}')
        self.add_column('formula', formula.format(), formula.format('html'))
        self.add_column('stoichiometry', s11y.format(), s11y.format('html'))
        self.add_column('uid', uid, uid)

    def add_column(self, name: str, value, html: str) -> None:
        assert name not in self.values
        self.values[name] = value
        self.html_reprs[name] = html

    def __getitem__(self, name):
        return self.html_reprs[name]

    def check(self, func):
        return func(self.count, self.values)


class Materials:
    def __init__(self, materials, panels):
        self.column_names = {
            'formula': 'Formula',
            'volume': 'volume [Å<sup>3</sup>]',
            'stoichiometry': 'Stoichiometry',
            'uid': 'Unique ID'}

        for panel in panels:
            assert panel.column_names.keys().isdisjoint(self.column_names)
            self.column_names.update(panel.column_names)

        self._materials = {}
        for material in materials:
            for panel in panels:
                panel.update_column_data(material)
            self._materials[material.uid] = material

        self.panels = panels

    def get_callbacks(self):
        callbacks = {}
        for panel in self.panels:
            callbacks.update(panel.callbacks)
        return callbacks

    def stoichiometries(self):
        s = []
        for material in self._materials.values():
            s.append((material.values['stoichiometry'],
                      material['stoichiometry']))
        return s

    def __getitem__(self, uid):
        return self._materials[uid]

    def get_rows(self, session: Session):
        filter = session.filter
        if session.stoichiometry != 'Any':
            if filter:
                filter += ','
            filter += f'stoichiometry={session.stoichiometry}'
        func = parse(filter)
        rows = self._materials.values()

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
        table = [(material.uid,
                  [material[name] for name in session.columns])
                 for material in rows]
        return (table,
                [(name, self.column_names[name]) for name in session.columns],
                pages,
                {(name, value) for name, value in self.column_names.items()
                 if name not in session.columns})
