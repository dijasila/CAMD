"""Get convex-hull data from OQMD123.

python -m cxdb.c2db.oqmd123 <path-to-oqmd123.db-file>

The OQMD123 database file can be obtained here:

  https://cmr.fysik.dtu.dk/_downloads/oqmd123.db
"""

import json
import sys
import gzip

from ase.db import connect


def read_oqmd123_data() -> tuple[dict[str, float],
                                 list[tuple[dict[str, int], float, str]]]:
    with gzip.open('oqmd.json.gz', 'rt') as fd:
        data = json.load(fd)
    return (data['atomic_energies'],
            [tuple(x) for x in data['formation_energies']])


def main():
    hform = []
    atomic_energies = {}
    for row in connect(sys.argv[1]).select():
        count = row.count_atoms()
        if len(count) == 1:
            symb = row.symbols[0]
            assert symb not in atomic_energies
            atomic_energies[symb] = row.energy / row.natoms
        hform.append((count, row.hform * row.natoms))
    with gzip.open('oqmd.json.gz', 'wt') as fd:
        json.dump({'atomic_energies': atomic_energies,
                   'formation_energies': hform},
                  fd, indent='')


if __name__ == '__main__':
    main()
