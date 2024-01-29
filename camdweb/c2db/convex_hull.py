import json
from pathlib import Path

from ase.formula import Formula

from camdweb.c2db.oqmd123 import db2json, read_oqmd123_data
from camdweb.panels.convex_hull import (calculate_ehull_energies,
                                        group_references)


def read_chull_data(oqmd_path: Path) -> tuple[dict[str, float],
                                              dict[str, tuple[dict[str, int],
                                                              float]]]:
    atomic_energies, refs = read_oqmd123_data(oqmd_path)
    return atomic_energies, refs


def update_chull_data(atomic_energies: dict[str, float],
                      refs: dict[str, tuple[dict[str, int], float]]) -> None:
    """Update ehull and hform values.

    Will calculate ehull and hform energies and insert into::

      A*/*/*/data.json

    Also write json files containing all
    reference data for each convex-hull plot::

      convex_hulls/MoS.json
                  /Mo.json
                  /S.json
                  /...
    """
    print('References:', len(refs))
    paths = {}
    c2db_uids = set()
    root = Path()
    for path in root.glob('A*/*/*/'):
        data = json.loads((path / 'data.json').read_text())
        uid = data['uid']
        paths[uid] = path
        f = Formula(path.parent.name)
        count = f.count()
        hform = data['energy'] - sum(n * atomic_energies[symbol]
                                     for symbol, n in count.items())
        refs[uid] = (count, hform)
        c2db_uids.add(uid)
    print('Materials:', len(c2db_uids))

    # Find all convex-hulls:
    # Create map from uid to sorted tuple of symbols:
    tmp: dict[str, tuple[str, ...]] = {}
    for uid, (count, hform) in refs.items():
        tmp[uid] = tuple(sorted(count))
    assert len(tmp) == len(refs)
    groups = group_references(tmp, c2db_uids)
    print('Convex hulls:', len(groups))

    # Write all hull files:
    folder = root / 'convex-hulls'
    folder.mkdir(exist_ok=True)
    for symbols, uids in groups.items():
        data = {}
        for uid in uids:
            source = 'C2DB' if uid in c2db_uids else 'OQMD'
            data[uid] = refs[uid] + (source,)
        (folder / (''.join(symbols) + '.json')).write_text(
            json.dumps(data, indent=2))

    # Calculate ehull:
    ehull_energies = {}
    for symbols, uids in groups.items():
        ehull_energies.update(
            calculate_ehull_energies({uid: refs[uid] for uid in uids},
                                     c2db_uids))

    # Update data.json files:
    for uid, path in paths.items():
        path /= 'data.json'
        data = json.loads(path.read_text())
        count, hform = refs[uid]
        natoms = sum(count.values())
        data['ehull'] = ehull_energies[uid] / natoms
        data['hform'] = hform / natoms
        path.write_text(json.dumps(data, indent=2))
