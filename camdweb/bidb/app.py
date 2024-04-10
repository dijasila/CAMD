from __future__ import annotations

import json
from pathlib import Path

import rich.progress as progress

from camdweb.html import table
from camdweb.materials import Material, Materials
from camdweb.panels.atoms import AtomsPanel
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


class BiDBAtomsPanel(AtomsPanel):
    column_descriptions = COLUMN_DESCRIPTIONS

    def update_material(self, material):
        data = json.loads((material.folder / 'data.json').read_text())
        material.columns.update(data)

    def create_column_one(self,
                          material: Material) -> str:
        rows = self.table_rows(material, COLUMN_DESCRIPTIONS)
        return table(None, rows)


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
