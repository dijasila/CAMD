import copy
from pathlib import Path

import ase.io.ulm as ulm
import numpy as np
from ase import Atoms
from ase.build import mx2
from ase.io.trajectory import write_atoms
from ase.spectrum.band_structure import BandStructure
from ase.io.jsonio import encode


def create_data(dir: Path, atoms: Atoms) -> None:
    from asr.bandstructure import Result

    # Create fake gpw-file:
    writer = ulm.Writer(dir / 'gs.gpw', tag='gpaw')
    with writer:
        writer.write(version=4)
        write_atoms(writer.child('atoms'), atoms)
        writer.child('results').write(energy=-5.9)
        wfs = writer.child('wave_functions')
        wfs.write(fermi_levels=np.array([-0.123]))

    (dir / 'results-asr.magstate.json').write_text('{"magstate": "NM"}')
    (dir / 'results-asr.structureinfo.json').write_text(
        '{"kwargs": {"data": {"has_inversion_symmetry": false}}}')
    (dir / 'results-asr.gs.json').write_text(
        '{"kwargs": {"data": {"gap": 1.8, "evac":4.5}}}')
    (dir / 'results-asr.gs@calculate.json').write_text(
        '{}')
    (dir / 'results-asr.convex_hull.json').write_text(
        '{"kwargs": {"data": {"hform": -0.9}}}')
    (dir / 'results-asr.database.material_fingerprint.json').write_text(
        '{"kwargs": {"data": {"uid": "MoS2-b3b4685fb6e1"}}}')

    kpts = atoms.cell.bandpath('GK', npoints=5)
    bs = BandStructure(kpts, np.zeros((1, 5, 2)), reference=-0.5)
    nosoc = copy.deepcopy(bs.todict())
    nosoc['efermi'] = -0.5
    soc = bs.todict()
    soc['efermi'] = -0.5
    soc['energies'] = np.zeros((4, 5))
    soc['sz_mk'] = np.zeros((4, 5))
    dct = Result.fromdata(bs_soc=soc, bs_nosoc=nosoc).todict()
    (dir / 'results-asr.bandstructure.json').write_text(encode(dct))


def create_tree(dir: Path):
    path = dir / 'MoS2'
    path.mkdir()
    atoms = mx2('MoS2')
    atoms.center(vacuum=5, axis=2)
    create_data(path, atoms)


if __name__ == '__main__':
    create_tree(Path())
