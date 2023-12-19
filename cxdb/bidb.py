import json
from collections import defaultdict
from pathlib import Path

from ase.db import connect
from ase.formula import Formula
from cxdb.atoms import AtomsPanel
from cxdb.material import Material, Materials
from cxdb.web import CXDBApp


def expand(db_file: str) -> None:
    monolayers: defaultdict[str, dict[str, int]] = defaultdict(dict)
    for row in connect(db_file).select():
        f = Formula(row.formula)
        ab, xy, n = f.stoichiometry()
        n //= row.number_of_layers
        name = f'{ab}/{n}{xy}'
        ids = monolayers[name]
        if row.monolayer_uid in ids:
            i = ids[row.monolayer_uid]
        else:
            i = len(ids)
            ids[row.monolayer_uid] = i
        folder = Path(name + f'-{i}')
        if row.number_of_layers == 1:
            folder /= 'monolayer'
        else:
            folder /= row.bilayer_uid.split('-', 2)[2]
        folder.mkdir(exist_ok=True, parents=True)
        row.toatoms().write(folder / 'structure.xyz')
        (folder / 'data.json').write_text(json.dumps(row.key_value_pairs))


class BilayerAtomsPanel(AtomsPanel):
    column_names = AtomsPanel.column_names | {
        'binding_energy_zscan': 'Binding energy (zscan)',
        'number_of_layers': 'Number of layers',
        'monolayer_uid': 'Monolayer ID',
        'bilayer_uid': 'Bilayer ID',
        'dynamically_stable': 'Dynamically stable',
        'magnetic': 'Magnetic',
        'interlayer_magnetic_exchange': 'Interlayer Magnetic State',
        'slide_stability': 'Slide Stability',
        'binding_energy_gs': 'Binding Energy (gs) [meV/Å<sup>2</sup>]',
        'ehull': 'Energy above convex hull [eV/atom]',
        'gap_pbe': 'Band gap (PBE)',
        'icsd_id': 'ICSD id of parent bulk structure',
        'cod_id': 'COD id of parent bulk structure',
        'layer_group': 'Layer group',
        'layer_group_number': 'Layer group number',
        'space_group': 'Space group',
        'space_group_number': 'Space group number'}

    columns = list(column_names)

    def get_column_data(self, material):
        dct = json.loads((material.folder / 'data.json').read_text())
        data = {}
        for key, value in dct.items():
            if key not in self.column_names:
                continue
            if key in {'space_group_number', 'space_group_number',
                       'cod_id', 'icsd_id'}:
                value = str(value)
            data[key] = (value,
                         str(value) if not isinstance(value, float)
                         else f'{value:.3f}')
        return data


def main(root: Path) -> CXDBApp:
    mlist: list[Material] = []
    for f in root.glob('*/*/*/'):
        if f.name == 'monolayer':
            uid = f.parent.name
        else:
            uid = f'{f.parent.name}-{f.name}'
        if len(mlist) % 20 == 0:
            print(end='.', flush=True)
        mlist.append(Material(f, uid))
    print()

    panels = [BilayerAtomsPanel()]
    materials = Materials(mlist, panels)

    initial_columns = {'uid', 'energy', 'formula'}

    return CXDBApp(materials, initial_columns, root)


if __name__ == '__main__':
    if 0:
        expand('bidb.db')
    main(Path()).app.run(host='0.0.0.0', port=8081, debug=True)
