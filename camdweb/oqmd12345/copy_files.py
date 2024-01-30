"""OQMD12345 web-app.

This module has code to convert ~oqmd*/tree/ folders and friends
(see PATTERNS variable below) to canonical tree layout.

Also contains a simple web-app that can run off the tree of folders.

python -m camdweb.oqmd12345.copy_files
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from ase.io import read

if TYPE_CHECKING:
    from taskblaster.repository import Repository

RESULT_FILES = [
    'output.json',
    'input.json']

# the directory name you want to gather data from
PATTERNS = [
    'material', 'gsresults',
    'relax3d_magnetic', 'postprocess']

ROOT = Path('/home/tara/webpage/tree-crysp')


def copy_materials(root: Path = ROOT,
                   patterns: list = PATTERNS) -> None:  # pragma: no cover
    from taskblaster.repository import Repository
    from taskblaster.state import State
    # glob over root to find patterns
    with Repository.find(root.as_posix()) as repo:
        for node in repo.tree([root.as_posix()]).nodes():
            # collect only done calculations
            if node.state is not State.done:
                continue
            future = repo.cache[node.name]
            record = future.value()
            if Path(record.node.name).name in patterns:
                copy_material(record)


def copy_material(record: Repository) -> None:  # pragma: no cover
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
        # shutil_files.append([traj.name, traj])
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

    folder.mkdir(exist_ok=True, parents=True)

    # copy relevant data to folder
    for file in shutil_files:
        filename, filepath = file
        shutil.copyfile(filepath, folder / filename)
    if write_structure:
        atoms.write(folder / 'structure.xyz')

    with open((folder / 'data.json'), 'r') as fd:
        existing_data = json.load(fd)
    data.update(existing_data)
    (folder / 'data.json').write_text(json.dumps(data, indent=0))


if __name__ == '__main__':
    if len(sys.argv) == 1:
        copy_materials(ROOT, PATTERNS)
    else:
        copy_materials(Path(sys.argv[1]), sys.argv[2:])
