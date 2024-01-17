"""CRYSP web-app.

This module has code to convert ~crysp*/tree/ folders and friends
(see PATTERNS variable below) to canonical tree layout.

Also contains simple web-app that can run off the tree of folders.

Goal is to have the code decoupled from ASE, GPAW and ASR.
Right now ASR webpanel() functions are still used (see cxdb.asr_panel module).
"""
from __future__ import annotations

import json
import shutil
import numpy as np
from monty.json import jsanitize
from collections import defaultdict
from pathlib import Path

import rich.progress as progress
from ase.io import read

from cxdb.panels.atoms import AtomsPanel
from cxdb.material import Material, Materials
from cxdb.panels.panel import Panel
from cxdb.web import CXDBApp
from taskblaster.repository import Repository
from taskblaster.state import State

RESULT_FILES = [
    'output.json',
    'input.json']

# the directory name you want to gather data from
PATTERNS = [
    'material',
    'groundstate', 'gsresults',
    'relax3d_magnetic', 'postprocess',
    ]


def copy_all_crysp_materials():  # pragma: no cover
    """Copy Crysp files to uniform tree structure.

    Tree structure::

       <stoichiometry>/<formula-units><fomula>/<oqmd_id>/

    Example::

       AB2/1MoS2/1/
       AB2/1MoS2/2/
       ...

    Build tree like this::

        $ cd /tmp
        $ mkdir tree
        $ cd tree
        $ python -c "from cxdb.crysp import *; copy_all_crysp_materials()"

    """
    root = Path('/home/tara/webpage/tree-crysp')  # tree location
    copy_materials(root)


def copy_materials(root: Path) -> None:
    # glob over root to find patterns
    names: defaultdict[str, int] = defaultdict(int)
    with Repository.find(root.as_posix()) as repo:
        for node in repo.tree([root.as_posix()]).nodes():
            # collect only done calculations
            if node.state is not State.done:
                continue
            future = repo.cache[node.name]
            record = future.value()
            if Path(record.node.name).name in PATTERNS:
                copy_material(record, names)


def copy_material(record: Repository, names: defaultdict[str, int]) -> None:
    """
    This is the function we use to convert TaskBlaster trees into an easy
    webpage displayable layout. From dir path, read the output/input.json.
    """
    def get_name(atoms):
        f = atoms.symbols.formula
        ab, xy, n = f.stoichiometry()
        return f'{ab}/{n}{xy}'

    # XXX because this is a work in progress, I haven't added a postprocess
    # for each step I would like to store data for various calculations. So
    # here I am going to make something ugly that can be later changed to
    # only pull postprocess directories which will contain the relevant
    # history one would like to use to create webpages.
    input, output, node = record.inputs, record.output, record.node
    write_structure = False
    shutil_files = []

    # this could be moved to actions on taskblaster since actions are supported
    # on tasks. It would be a programmatic way of getting data from tasks to go
    # into the database.
    if Path(node.name).name == 'material':
        oqmd, atoms = Path(node.name).parts[1], output.atoms
        data = {'oqmd_entry_id': output.oqmd_entry_id}
    if Path(node.name).name == 'relax3d_magnetic':
        atoms, traj = output['atoms'], output['traj_path']
        oqmd, write_structure = Path(node.name).parts[1], True
        shutil_files.append([traj.name, traj])
        if traj.name == 'relax.traj':
            data = {'magstate': 'Mag'}
        else:
            data = {'magstate': 'NM'}
    if Path(node.name).name == 'postprocess':
        oqmd = Path(node.name).parts[1]
        data = {
            'fmax': np.linalg.norm(output['forces']).max(),
            'smax': abs(output['stresses']).max()}
        atoms = output['atoms']
    if Path(node.name).name == 'groundstate':
        oqmd, atoms = Path(node.name).parts[1], input['atoms']
        shutil_files.append([output.name, output.as_posix()])
        data = {}
    if Path(node.name).name == 'gsresults':
        oqmd, atoms = Path(node.name).parts[1], read(input['groundstate'])
        data = {
            'etot': output['etot'],
            'gap': output['gap'],
            'gap_dir': output['gap_dir']}
    # XXX action to webpage: ref_id, atoms, data consisting of kvp (flat
    # namespace?), copy files

    name = get_name(atoms=atoms)
    folder = Path(name) / str(oqmd)
    names[name] = oqmd

    folder.mkdir(exist_ok=True, parents=True)

    # copy relevant data to folder
    for file in shutil_files:
        filename, filepath = file
        shutil.copyfile(filepath, folder / filename)
    if write_structure:
        atoms.write(folder / f'structure.xyz')

    data = jsanitize(data)
    try:
        with open((folder / 'data.json'), 'r') as file:
            existing_data = json.load(file)
        data.update(existing_data)
    except:
        ...

    (folder / 'data.json').write_text(json.dumps(data, indent=0))


class CRYSPAtomsPanel(AtomsPanel):
    def __init__(self):
        super().__init__()
        self.column_names.update(
            etot='Total Energy [eV]',
            magstate='Magnetic state',
            fmax='Maximum Force [eV/A²]',
            smax='Maximum Stress',
            gap='Band gap (PBE) [eV]',
            gap_dir='Direct Band Gap (PBE) [eV]',
            oqmd_entry_id='OQMD Ref. ID')
        self.columns = list(self.column_names)

    def update_data(self, material: Material):
        super().update_data(material)
        data = json.loads((material.folder / 'data.json').read_text())
        for key, value in data.items():
            material.add_column(key, value)


def main(root: Path) -> CXDBApp:
    """Create CRYSP app."""
    mlist: list[Material] = []
    files = list(root.glob('A*/*/*/'))
    with progress.Progress() as pb:
        pid = pb.add_task('Reading materials:', total=len(files))
        for f in files:
            uid = f'{f.parent.name}-{f.name}'
            mlist.append(Material.from_file(f / 'structure.xyz', uid))
            pb.advance(pid)

    panels: list[Panel] = [CRYSPAtomsPanel()]
    # for name in ['bandstructure',
    #              'phonons',
    #              'bader']:
    #     panels.append(ASRPanel(name))

    # breakpoint()
    materials = Materials(mlist, panels)
    initial_columns = ['gap', 'formula']

    return CXDBApp(materials, initial_columns, root)


if __name__ == '__main__':
    main(Path()).app.run(host='0.0.0.0', port=8081, debug=True)
