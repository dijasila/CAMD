from pathlib import Path

import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.db import connect
from cxdb.cmr.projects import ProjectDescription, projects


def create_db_files(path: Path) -> dict[str, ProjectDescription]:
    pds = {}
    for name, cls in projects.items():
        pd = cls()
        pds[name] = pd
        dbfile = path / f'{name}.db'
        with connect(dbfile) as db:
            atoms = Atoms('H2', [(0, 0, 0), (0.7, 0, 0)])
            atoms.center(vacuum=1)
            if pd.pbc is not None:
                atoms.pbc = pd.pbc
            kwargs: dict[str, int | str | float | bool] = {}
            data: dict | None = None
            if pd.uid != 'id':
                kwargs[pd.uid] = 'id-1'
            # int id's should be converted to strings
            if name == 'mp_gllbsc':
                kwargs['mpid'] = 50000
                kwargs['icsd_id'] = 70000
            elif name == 'pv_pec_oqmd':
                kwargs['icsd'] = 70000
            elif name == 'oqmd123':
                kwargs['oqmd_id'] = 'x30000'
                atoms.calc = SinglePointCalculator(
                    atoms=atoms,
                    energy=-27.0,
                    forces=np.zeros((2, 3)),
                    stress=np.zeros(6),
                    magmom=2.3)
            elif name == 'abx2':
                kwargs['KS_gap'] = 0  # int should be converted to float
                kwargs['E_hull'] = np.inf  # should not be shown
            elif name == 'abs3':  # band-structure data:
                data = {'x': [0, 1, 2],
                        'y': [1.0, 1.5, 1.4],
                        'X': [0, 1],
                        'names': ['G', 'X']}
            elif name == 'bidb':
                kwargs['number_of_layers'] = 2
                kwargs['monolayer_uid'] = 'abc-123'
            elif name == 'lowdim':
                kwargs['source'] = 'COD'
                kwargs['doi'] = 'asdf'
                kwargs['dbid'] = 'a1'

            db.write(atoms, abc=27.3, gap=0, **kwargs, data=data)

            if name == 'abs3':
                db.write(atoms)
            elif name == 'bidb':
                db.write(atoms, uid='id-2', number_of_layers=1)
            elif name == 'lowdim':
                db.write(atoms, source='ICSD', dbid='a2')

    return pds


if __name__ == '__main__':
    create_db_files(Path())
