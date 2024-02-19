import copy
from pathlib import Path

import ase.io.ulm as ulm
import numpy as np
from ase import Atoms
from ase.build import mx2
from ase.io.jsonio import encode
from ase.io.trajectory import write_atoms
from ase.spectrum.band_structure import BandStructure
from asr.bandstructure import Result


def create_data(dir: Path, atoms: Atoms) -> None:

    # Create fake gpw-file:
    writer = ulm.Writer(dir / 'gs.gpw', tag='gpaw')
    with writer:
        writer.write(version=4)
        write_atoms(writer.child('atoms'), atoms)
        writer.child('results').write(energy=-5.9)
        wfs = writer.child('wave_functions')
        wfs.write(fermi_levels=np.array([-0.123]))

    # Create fake ASR result-files:
    (dir / 'results-asr.magstate.json').write_text('{"magstate": "NM"}')
    (dir / 'results-asr.magnetic_anisotropy.json').write_text(
        '{"spin_axis": "z"}')
    (dir / 'results-asr.structureinfo.json').write_text(
        """{"has_inversion_symmetry": false,
            "layergroup": "p-6m2",
            "lgnum": 78,
            "spglib_dataset":
                {"rotations":{"__ndarray__":[
                                  [1, 3, 3],
                                  "int32",
                                  [1,0,0,0,1,0,0,0,1]]}}}""")
    (dir / 'results-asr.gs.json').write_text(
        """{"kwargs": {"data": {"gap": 1.8,
                                "gap_dir": 1.8,
                                "gap_dir_nosoc": 1.9,
                                "k_cbm_c": [0.0, 0.0, 0.0],
                                "k_vbm_c": [0.1, 0.0, 0.0],
                                "evac": 4.5,
                                "efermi": 1.5,
                                "gaps_nosoc": {"vbm": 0.5, "cbm": 2.5}}}}""")
    (dir / 'results-asr.gs@calculate.json').write_text(
        '{}')
    (dir / 'results-asr.database.material_fingerprint.json').write_text(
        '{"kwargs": {"data": {"uid": "MoS2-b3b4685fb6e1"}}}')

    # Band-structure:
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

    # Phonons:
    dct = {'minhessianeig': 0.001,
           'dynamic_stability_phonons': 'high',
           'path': kpts,
           'omega_kl': np.zeros((2, 7)),
           'q_qc': np.array([[0, 0, 0], [0, 0.5, 0]]),
           'interp_freqs_kl': np.zeros((5, 7))}
    (dir / 'results-asr.phonons.json').write_text(encode(dct))

    # Bader:
    (dir / 'results-asr.bader.json').write_text("""
    {"kwargs":
      {"data":
       {"bader_charges":
        {"__ndarray__": [[3], "float64", [1.24, -0.62, -0.62]]},
        "sym_a": ["Mo", "S", "S"]}}}""")

    # Shift:
    freqs = np.linspace(0, 10, 11)
    sigma = freqs * (1 + 0j)
    dct = {
        'symm': {'zero': 'xxx=xxz=...', 'yyy': 'yyy=-xyx=-yxx=-xxy'},
        'sigma': {'yyy': sigma},
        'freqs': freqs}
    (dir / 'results-asr.shift.json').write_text(encode(dct))


def create_tree(dir: Path):
    path = dir / 'MoS2'
    path.mkdir()
    atoms = mx2('MoS2')
    atoms.center(vacuum=5, axis=2)
    create_data(path, atoms)


if __name__ == '__main__':
    create_tree(Path())
