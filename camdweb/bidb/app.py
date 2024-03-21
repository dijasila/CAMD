from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import rich.progress as progress
from ase.db import connect
from ase.formula import Formula

from camdweb.html import table
from camdweb.materials import Material, Materials
from camdweb.panels.atoms import AtomsPanel
from camdweb.panels.panel import Panel, PanelData, SkipPanel
from camdweb.web import CAMDApp

COLUMN_DESCRIPTIONS = {
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


class BiDBAtomsPanel(AtomsPanel):
    column_descriptions = COLUMN_DESCRIPTIONS

    def update_material(self, material):
        data = json.loads((material.folder / 'data.json').read_text())
        material.columns.update(data)

    def create_column_one(self,
                          material: Material) -> str:
        rows = self.table_rows(material, COLUMN_DESCRIPTIONS)
        return table(None, rows)


class StackingsPanel(Panel):
    def get_data(self,
                 material: Material) -> PanelData:
        bilayers = material.data.get('bilayers')
        if bilayers is None:
            raise SkipPanel
        rows = []
        for uid, bilayer in bilayers.items():
            e = bilayer.binding_energy_zscan
            rows.append([f'<a href="{uid}">{uid}</a>', f'{e:.3f}'])
        tbl = table(['Stacking', 'Binding energy [meV/Å<sup>2</sup>]'], rows)
        return PanelData(tbl, title='Stackings')


def main(root: Path) -> CAMDApp:
    mlist: list[Material] = []
    paths = list(root.glob('*/*/monolayer/'))
    with progress.Progress() as pb:
        pid = pb.add_task('Reading matrerials:', total=len(paths))
        for f1 in paths:
            uid1 = f1.parent.name
            bilayers = {}
            for f2 in f1.parent.glob('*/'):
                if f2.name != 'monolayer':
                    uid2 = f'{f2.parent.name}-{f2.name}'
                    bilayer = Material.from_file(f2 / 'structure.xyz', uid2)
                    bilayers[uid2] = bilayer
                    mlist.append(bilayer)
            monolayer = Material.from_file(f1 / 'structure.xyz', uid1)
            monolayer.data['bilayers'] = bilayers
            mlist.append(monolayer)
        pb.advance(pid)

    panels = [BiDBAtomsPanel(),
              StackingsPanel()]

    materials = Materials(mlist, panels)

    initial_columns = ['uid', 'area', 'formula']

    return CAMDApp(materials, initial_columns, root=root)


if __name__ == '__main__':
    main(Path()).app.run(host='0.0.0.0', port=8083, debug=True)
