import argparse
from pathlib import Path

from ase import Atoms
from ase.io import read

from cxdb.material import Material, Materials
from cxdb.panels.atoms import AtomsPanel
from cxdb.web import CXDBApp


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', nargs='+')
    args = parser.parse_args(argv)
    rows: list[Material] = []
    i = 1
    for filename in args.filename:
        path = Path(filename)
        configs = read(path)
        if isinstance(configs, Atoms):
            configs = [configs]
            for atoms in configs:
                print(atoms)
                rows.append(Material(path.parent, str(i), atoms))
                i += 1

    panels = [AtomsPanel()]
    materials = Materials(rows, panels)

    initial_columns = ['uid', 'formula']

    root = Path.cwd()

    CXDBApp(materials, initial_columns, root).app.run(
        host='0.0.0.0', port=8081, debug=True)


if __name__ == '__main__':
    main()
