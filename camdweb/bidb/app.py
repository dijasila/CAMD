from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Generator

import rich.progress as progress
from ase.atoms import Atoms
from ase.db import connect
from ase.formula import Formula
from ase.io import read

from camdweb.html import table
from camdweb.materials import Material, Materials
from camdweb.panels.atoms import AtomsPanel
from camdweb.panels.panel import Panel
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


class BiDBMaterial(Material):
    binding_energy_zscan: float

    def __init__(self,
                 folder: Path,
                 uid: str,
                 bilayers: dict[str, BiDBMaterial] | None = None):
        atoms = read(folder / 'structure.xyz')
        assert isinstance(atoms, Atoms)
        super().__init__(folder, uid, atoms)
        data = json.loads((folder / 'data.json').read_text())
        for key in COLUMN_DESCRIPTIONS:
            setattr(self, key, data.get(key))
        self.bilayers = bilayers

    def get_columns(self):
        columns = super().get_columns()
        for key in COLUMN_DESCRIPTIONS:
            value = getattr(self, key)
            if value is not None:
                columns[key] = value
        return columns


class BiDBAtomsPanel(AtomsPanel):
    def create_column_one(self,
                          material: Material) -> str:
        rows = []
        for key, desc in COLUMN_DESCRIPTIONS.items():
            value = getattr(material, key)
            rows.append([desc, material.html_format_column(key, value)])

        return table(None, rows)


class StackingsPanel(Panel):
    title = 'Stackings'

    def get_html(self,
                 material: Material) -> Generator[str, None, None]:
        assert isinstance(material, BiDBMaterial)
        bilayers = material.bilayers
        if bilayers is None:
            return
        rows = []
        for uid, bilayer in bilayers.items():
            e = bilayer.binding_energy_zscan
            rows.append([f'<a href="{uid}">{uid}</a>', f'{e:.3f}'])
        tbl = table(['Stacking', 'Binding energy [meV/Å<sup>2</sup>]'], rows)
        yield tbl


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
                    bilayer = BiDBMaterial(f2, uid2)
                    bilayers[uid2] = bilayer
                    mlist.append(bilayer)
            mlist.append(BiDBMaterial(f1, uid1, bilayers))
        pb.advance(pid)

    panels = [BiDBAtomsPanel(),
              StackingsPanel()]

    materials = Materials(mlist, panels, COLUMN_DESCRIPTIONS)

    initial_columns = ['uid', 'area', 'formula']

    return CAMDApp(materials, initial_columns, root)


if __name__ == '__main__':
    main(Path()).app.run(host='0.0.0.0', port=8083, debug=True)
