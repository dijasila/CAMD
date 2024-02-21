from pathlib import Path

from ase.formula import Formula
from camdweb.c2db.copy import PATTERNS, read_result_file

fro = Path('/home/niflheim/joaso/paper_calcs/htp-spiral_paper/n1-magnets/tree')
to = Path('/home/niflheim2/cmr/C2DB-ASR/')


def copy_spinspiral_data_to_tree():  # pragma: no cover
    dirs = [dir
            for pattern in PATTERNS
            for dir in to.glob(pattern)
            if dir.name[0] != '.']
    paths = {}
    for dir in dirs:
        fp = dir / 'results-asr.database.material_fingerprint.json'
        try:
            olduid = read_result_file(fp)['uid']
        except FileNotFoundError:
            print('No fingerprint:', fp)
            continue
        f = Formula(olduid.split('-')[0])
        stoi, reduced, nunits = f.stoichiometry()
        path = fro / f'{stoi}/{f}/{olduid}'
        paths[path] = dir

    print(len(paths))
    for path in fro.glob('A*/*/*/'):
        if path not in paths:
            print(path)
            continue
        dir = paths.pop(path)
        for name in ['collect_spiral',
                     'dmi',
                     'spinorbit',
                     'spinorbit@calculate']:
            f1 = path / f'results-asr.{name}.json'
            if f1.is_file():
                f2 = dir / f1.name
                assert not f2.is_file()
                if 1:
                    f2.write_bytes(f1.read_bytes())
                else:
                    print(f1, f2)
    print(len(paths))


if __name__ == '__main__':
    copy_spinspiral_data_to_tree()
