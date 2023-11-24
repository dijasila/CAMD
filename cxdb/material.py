from functools import cached_property

from ase.io import read


class Material:
    header = ['formula',
              'energy [eV]',
              'volume [Å<sup>3</sup>]',
              '']

    def __init__(self, folder, id):
        self.folder = folder
        self.id = id
        self.atoms = read(folder / 'rlx.traj')
        self.energy = self.atoms.get_potential_energy()
        formula = self.atoms.symbols.formula.convert('periodic')
        self.formula_html = formula.format('html')

    @cached_property
    def columns(self) -> list[str]:
        return [
            self.formula_html,
            f'{self.energy:.3f}',
            f'{self.atoms.get_volume():.3f}']
