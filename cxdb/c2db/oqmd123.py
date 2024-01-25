"""Get convex-hull data from OQMD123.

python -m cxdb.c2db.oqmd123 <path-to-oqmd123.db-file>

The OQMD123 database file can be obtained here:

  https://cmr.fysik.dtu.dk/_downloads/oqmd123.db
"""

import gzip
import json
import sys
from pathlib import Path

from ase.db import connect


def read_oqmd123_data() -> tuple[dict[str, float],
                                 dict[str, tuple[dict[str, int], float]]]:
    with gzip.open('oqmd.json.gz', 'rt') as fd:
        data = json.load(fd)
    return (data['atomic_energies'],
            {uid: tuple(x)  # type: ignore
             for uid, x in data['formation_energies'].items()})


def main(oqmd_db_file: Path):
    hform = {}
    atomic_energies = {}
    for row in connect(oqmd_db_file).select():
        count = row.count_atoms()
        if len(count) == 1:
            symb = row.symbols[0]
            assert symb not in atomic_energies
            atomic_energies[symb] = row.energy / row.natoms
        hform[row.uid] = (count, row.hform * row.natoms)
    with gzip.open('oqmd.json.gz', 'wt') as fd:
        json.dump({'atomic_energies': atomic_energies,
                   'formation_energies': hform},
                  fd, indent=0)


if __name__ == '__main__':
    main(Path(sys.argv[1]))
