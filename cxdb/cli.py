import argparse
from pathlib import Path

from ase import Atoms
from ase.io import read

from cxdb.material import Material, Materials
from cxdb.panels.atoms import AtomsPanel
from cxdb.web import CXDBApp


def main(argv: list[str] | None = None,
         run=True) -> CXDBApp:
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', nargs='+',
                        help='Filename of atomic structure file.')
    args = parser.parse_args(argv)
    rows: list[Material] = []
    i = 1
    for filename in args.filename:
        path = Path(filename)
        configs = read(path)
        if isinstance(configs, Atoms):
            configs = [configs]
            for atoms in configs:
                rows.append(Material(path.parent, str(i), atoms))
                i += 1

    panels = [AtomsPanel()]
    materials = Materials(rows, panels)

    initial_columns = ['uid', 'formula']
    for key in ['length', 'area', 'volume']:
        if key in materials.column_names:
            initial_columns.append(key)

    root = Path.cwd()

    app = CXDBApp(materials, initial_columns, root)
    if run:
        app.app.run(host='0.0.0.0', port=8081, debug=True)
    return app


if __name__ == '__main__':
    main()
