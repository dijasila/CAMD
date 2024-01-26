import json
from pathlib import Path

from ase.formula import Formula
from ase.phasediagram import PhaseDiagram

from camdweb.c2db.oqmd123 import read_oqmd123_data, db2json
from camdweb.panels.convex_hull import group_references


def update_chull_data(root: Path) -> None:
    oqmd = root / 'oqmd123.json.gz'
    if not oqmd.is_file():
        db = root / 'oqmd123.db'
        if db.is_file():
            db2json(db, oqmd)
        else:
            raise FileNotFoundError(
                'Please download oqmd123.db file:\n\n'
                '   wget https://cmr.fysik.dtu.dk/_downloads/oqmd123.db\n')
    atomic_energies, refs = read_oqmd123_data(oqmd)
    print('References:', len(refs))

    paths = {}
    c2db_uids = set()
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

    tmp: dict[str, tuple[str, ...]] = {}
    for uid, (count, hform) in refs.items():
        tmp[uid] = tuple(sorted(count))
    assert len(tmp) == len(refs)
    groups = group_references(tmp, c2db_uids)
    print('Convex hulls:', len(groups))

    folder = root / 'convex-hulls'
    folder.mkdir(exist_ok=True)
    ehull_energies = {}
    for symbols, uids in groups.items():
        data = {}
        for uid in uids:
            source = 'C2DB' if uid in c2db_uids else 'OQMD'
            data[uid] = refs[uid] + (source,)
        (folder / (''.join(symbols) + '.json')).write_text(
            json.dumps(data, indent=2))

        pd = PhaseDiagram(
            [(count, hform) for (count, hform, source) in data.values()],
            verbose=False)
        for uid in uids:
            if uid in c2db_uids:
                count, hform, source = data[uid]
                if len(count) == len(symbols):  # pragma: no branch
                    ehull = hform - pd.decompose(**count)[0]
                    ehull_energies[uid] = ehull

    for uid, path in paths.items():
        path /= 'data.json'
        data = json.loads(path.read_text())
        count, hform = refs[uid]
        natoms = sum(count.values())
        data['ehull'] = ehull_energies[uid] / natoms
        data['hform'] = hform / natoms
        path.write_text(json.dumps(data, indent=2))


if __name__ == '__main__':
    update_chull_data(Path())
