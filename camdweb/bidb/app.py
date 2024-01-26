from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

from ase.db import connect
from ase.formula import Formula

from camdweb.material import Material, Materials
from camdweb.panels.atoms import AtomsPanel
from camdweb.panels.panel import Panel
from camdweb.html import table
from camdweb.web import CAMDApp


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
    def __init__(self):
        super().__init__()
        self.column_names.update(
            {'binding_energy_zscan': 'Binding energy (zscan)',
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
             'space_group_number': 'Space group number'})
        self.columns = list(self.column_names)

    def update_data(self, material):
        super().update_data(material)
        dct = json.loads((material.folder / 'data.json').read_text())
        for key, value in dct.items():
            if key not in self.column_names:
                continue
            if key in {'space_group_number', 'cod_id', 'icsd_id'}:
                value = str(value)
            material.add_column(key, value)


class StackingsPanel(Panel):
    title = 'Stackings'
    column_names = {'nstackings': 'Number of Stackings'}

    def get_html(self,
                 material: Material,
                 materials: Materials) -> tuple[str, str]:
        if material.number_of_layers == 2:
            return ('', '')
        rows = []
        for f in self.bilayer_folders(material):
            uid = f'{material.uid}-{f.name}'
            bilayer = materials[uid]
            rows.append([f'<a href="{uid}">{uid}</a>',
                         bilayer['binding_energy_zscan']])
        tbl = table(['Stacking', 'Binding energy'], rows)
        return tbl, ''

    def bilayer_folders(self, material) -> list[Path]:
        return [f for f in material.folder.glob('../*/')
                if f.name != 'monolayer']

    def update_data(self,
                    material: Material) -> None:
        if material.number_of_layers == 1:
            n = len(self.bilayer_folders(material))
            material.add_column('nstackings', n)


def main(root: Path) -> CAMDApp:
    mlist: list[Material] = []
    for f in root.glob('*/*/*/'):
        if f.name == 'monolayer':
            uid = f.parent.name
        else:
            uid = f'{f.parent.name}-{f.name}'
        if len(mlist) % 20 == 0:
            print(end='.', flush=True)
        mlist.append(Material.from_file(f / 'structure.xyz', uid))
    print()

    panels = [BilayerAtomsPanel(),
              StackingsPanel()]

    materials = Materials(mlist, panels)

    initial_columns = ['uid', 'area', 'formula']

    return CAMDApp(materials, initial_columns, root)


if __name__ == '__main__':
    main(Path()).app.run(host='0.0.0.0', port=8083, debug=True)
