"""Structural information.

We add some utils to oqmd12345 which should at some point be removed
either to ASE or another place since we don't want to depend on spglib or ASE.
"""
import re
import string
import numpy as np
from math import gcd
from collections.abc import Iterable
from ase import Atoms
from functools import reduce
from spglib import get_symmetry_dataset
from spglib.spglib import get_symmetry_layerdataset


class SymmetryAnalyzer:
    def __init__(self, atoms: Atoms, wrap: bool):
        """
        :param atoms: ASE atoms object to analyze symmetry for.
        :param wrap: Wrap atoms outside the unit cell into the cell in
            directions with periodic boundary conditions so the scaled
            coordinates are bounded between [0, 1] when passed to spglib.
        """
        self.atoms = atoms
        self.spglib_cell = self._to_spglib(wrap=wrap)
        self._dataset = {}

    def _to_spglib(self, wrap: bool):
        """
        Return a structure object that spglib understands from an ASE atoms
        object.

        :param wrap:
        :return:
        """
        lattice = self.atoms.get_cell().array
        positions = self.atoms.get_scaled_positions(wrap=wrap)
        numbers = self.atoms.get_atomic_numbers()

        return lattice, positions, numbers

    def generate_spglib_dataset(self, symprec: float, angle_tolerance: float,
                                hall_number: int = 0):
        self._dataset = get_symmetry_dataset(
            cell=self.spglib_cell, symprec=symprec,
            angle_tolerance=angle_tolerance, hall_number=hall_number)

        if self.atoms.pbc.sum() == 2:
            self._dataset.update(self.get_layer_group(symprec=symprec))

    def has_inversion(self):
        """
        Determine if atoms has inversion symmetry.

        :param rotations: list of 3x3 matrices List of point group symmetries.
        :return: Whether rotations contain inversion symmetry.
        """
        rotations = self.rotations
        if isinstance(rotations, Iterable):
            inversion = -np.identity(3, dtype=int)
            self._dataset['has_inversion_symmetry'] = np.any(
                [np.all(rotation == inversion) for rotation in rotations])
        else:
            raise ValueError("Please generate the spglib dataset prior to "
                             "checking if the dataset has inversions. e.g. "
                             "Obj.generate_spglib_dataset(...)")

    def __getitem__(self, item):
        return self._dataset[item]

    def __getattr__(self, item):
        if item in self._dataset:
            return self._dataset[item]
        else:
            raise AttributeError(f"'{type(self).__name__}' object has no attribute '{item}'")

    @property
    def dataset(self):
        return self._dataset

    @property
    def crystal_type(self):
        wyckoffs = ''.join(sorted(set(self.wyckoffs)))
        crystal_type = f'{self.atoms.symbols.formula.stoichiometry()[0]}-' \
                       f'{self.number}-{wyckoffs}'
        return crystal_type

    def get_layer_group(self, symprec: float):
        assert self.atoms.pbc.sum() == 2
        aperiodic_dir = np.where(~self.atoms.pbc)[0][0]
        # Prepare for spglib v3 API change to always have the aperiodic_dir==2
        # See: https://github.com/spglib/spglib/issues/314.
        if aperiodic_dir != 2:
            perm = np.array([0, 1, 2])
            # Swap axes such that aperiodic is always 2
            perm[2], perm[aperiodic_dir] = perm[aperiodic_dir], perm[2]

            atoms = self.atoms.copy()
            atoms.set_pbc(atoms.get_pbc()[perm])

            # Atoms are stored in cartesian coordinates, therefore, we're free
            # to permute the cell vectors, and the system remains invariant
            atoms.set_cell(atoms.get_cell()[perm], scale_atoms=False)
            aperiodic_dir = 2

        assert aperiodic_dir == 2
        lg_dct = get_symmetry_layerdataset(cell=(self.atoms.get_cell(),
                                                 self.atoms.get_scaled_positions(),
                                                 self.atoms.get_atomic_numbers()),
                                           symprec=symprec,
                                           aperiodic_dir=aperiodic_dir)
        layergroup_info = zip(['international', 'number'],  # lg dict ref key
                              ['layergroup', 'lgnum'])  # dataset key

        return {key: lg_dct[lg_key] for lg_key, key in layergroup_info}


class StructureInfo:
    def __init__(self, atoms: Atoms):
        self.atoms = atoms

    def get_spglib_analysis(self, symprec: float, angle_tolerance: float,
                            hall_number: int = 0, wrap: bool = False) \
            -> SymmetryAnalyzer:
        sa = SymmetryAnalyzer(atoms=self.atoms, wrap=wrap)
        sa.generate_spglib_dataset(symprec=symprec,
                                   angle_tolerance=angle_tolerance,
                                   hall_number=hall_number,
                                   )
        sa.has_inversion()

        return sa

    def formula(self, mode: str = 'metal', empirical: bool = False):
        return self.atoms.get_chemical_formula(mode=mode, empirical=empirical)

    def reduced_formula(self, mode: str = 'metal', empirical: bool = False,
                        anonymous: bool = True):
        formula = self.atoms.get_chemical_formula(mode=mode,
                                                  empirical=empirical)
        return self._reduce_formula_string(formula=formula, anonymous=anonymous)

    @property
    def get_anonymous_reduced_formula(self):
        formula = self.formula(mode='metal', empirical=False)
        return self._reduce_formula_string(formula=formula, anonymous=True)

    # XXX This just seems like a less comprehensive version of what the
    # Formula class is capable of doing. In what way is this meaningfully
    # different from Formula.stoichiometry()[0]. the results look identical
    def _reduce_formula_string(self, formula: str, anonymous: bool = True) \
            -> str:
        """
        Get reduced formula from formula according to the
        ase.atoms.get_chemical_formula().

        :param mode: Chose from all, reduce, hill, and metal.
        :param empirical: Divide the symbol counts by their greatest common
            divisor to yield an empirical formula. Only for mode `metal` and
            `hill`.
        :param anonymous: If True, return the anonymous stoichiometry, for
            example, "AB2" rather than "MoS2."
        :return: The reduced formula
        """
        split = re.findall('[A-Z][^A-Z]*', formula)
        matches = [re.match('([^0-9]*)([0-9]+)', x)
                   for x in split]
        numbers = [int(x.group(2)) if x else 1 for x in matches]
        symbols = [matches[i].group(1) if matches[i] else split[i]
                   for i in range(len(matches))]
        divisor = reduce(gcd, numbers)
        result = ''
        numbers = [x // divisor for x in numbers]
        numbers = [str(x) if x != 1 else '' for x in numbers]
        if anonymous:
            numbers = sorted(numbers)
            symbols = string.ascii_uppercase
        for symbol, number in zip(symbols, numbers):
            result += symbol + number
        return result
