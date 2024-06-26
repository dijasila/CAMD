from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from ase import Atoms

from camdweb import ColVal
from camdweb.utils import fft, read_atoms


class Material:
    def __init__(self,
                 uid: str,
                 atoms: Atoms | None = None,
                 folder: Path | None = None):
        """Object representing a material and associated data.

        >>> mat = Material('x1', atoms=Atoms('H2O'))
        >>> mat.formula
        'OH2'
        >>> mat.stoichiometry
        'AB2'
        """
        self.uid = uid
        self.atoms = atoms or Atoms()
        self.folder = folder or Path()

        self.data: dict[str, Any] = {}

        # Get number-of-atoms dicts:
        self.count, formula, reduced, stoichiometry = fft(self.atoms.numbers)

        self.columns: dict[str, ColVal] = {
            'uid': uid,
            'natoms': sum(self.count.values()),
            'nspecies': len(self.count),
            'formula': formula,
            'reduced': reduced,
            'stoichiometry': stoichiometry}
        pbc = self.atoms.pbc
        dims = pbc.sum()
        if dims > 0:
            vol = abs(np.linalg.det(self.atoms.cell[pbc][:, pbc]))
            name = ['length', 'area', 'volume'][dims - 1]
            self.columns[name] = vol

    def __repr__(self):
        return f'Material({self.uid}, {self.atoms}, {self.folder})'

    def __getattr__(self, name):
        if name.startswith('_'):
            raise AttributeError
        try:
            return self.columns[name]
        except KeyError:
            raise AttributeError(f'Material object has no attribute {name!r}')

    @classmethod
    def from_file(cls, file: Path, uid: str) -> Material:
        atoms = read_atoms(file)
        return cls(uid, atoms, file.parent)
