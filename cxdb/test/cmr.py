from pathlib import Path

from ase import Atoms
from ase.db import connect
from cxdb.cmr.projects import _projects, ProjectDescription


def create_db_files(path: Path) -> dict[str, ProjectDescription]:
    pds = {}
    for name, factory in _projects.items():
        pd = factory()
        pds[name] = pd
        dbfile = path / f'{name}.db'
        with connect(dbfile) as db:
            atoms = Atoms('H2', [(0, 0, 0), (0.7, 0, 0)])
            atoms.center(vacuum=1)
            if pd.pbc is not None:
                atoms.pbc = pd.pbc
            kwargs = {}
            if pd.uid != 'id':
                kwargs[pd.uid] = 'id-1'
            if name == 'mp_gllbsc':
                kwargs['mpid'] = 50000
                kwargs['icsd_id'] = 70000
            elif name == 'pv_pec_oqmd':
                kwargs['icsd'] = 70000
            elif name == 'oqmd123':
                kwargs['oqmd_id'] = 'x30000'
            elif name == 'abx2':
                kwargs['KS_gap'] = 0  # int should be converted to float
            db.write(atoms, abc=27.3, gap=0, **kwargs)
    return pds


if __name__ == '__main__':
    create_db_files(Path())
