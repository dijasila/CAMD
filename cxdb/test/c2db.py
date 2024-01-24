import copy
import json
from pathlib import Path

import ase.io.ulm as ulm
import numpy as np
from ase import Atoms
from ase.build import mx2
from ase.io.jsonio import encode
from ase.io.trajectory import write_atoms
from ase.spectrum.band_structure import BandStructure
from asr.bandstructure import Result

# Convex-hull example data:
CHULL = {
    'hform': -0.920927544752896,
    'references': [
        {
            'hform': 0.0,
            'formula': 'S48',
            'uid': 'S48',
            'natoms': 48,
            'title': 'Bulk crystals (from OQMD123)',
            'legend': 'Bulk crystals',
            'name': 'S48',
            'label': 'S48',
            'link': 'https://cmrdb.fysik.dtu.dk/oqmd123/row/S48',
            'method': 'DFT',
        },
        {
            'hform': 0.0,
            'formula': 'Mo',
            'uid': 'Mo',
            'natoms': 1,
            'title': 'Bulk crystals (from OQMD123)',
            'legend': 'Bulk crystals',
            'name': 'Mo',
            'label': 'Mo',
            'link': 'https://cmrdb.fysik.dtu.dk/oqmd123/row/Mo',
            'method': 'DFT',
        },
        {
            'hform': -0.18018421551178587,
            'formula': 'Mo2S2',
            'uid': 'Mo2S2-925d20f42e31',
            'natoms': 4,
            'title': 'Monolayers (from C2DB)',
            'legend': 'Monolayers',
            'name': 'Mo2S2 (AB-187-hi)',
            'label': 'Mo2S2 (AB-187-hi)',
            'link': '/c2db/row/Mo2S2-925d20f42e31',
            'method': 'DFT',
        },
        {
            'hform': -0.920927544752896,
            'formula': 'MoS2',
            'uid': 'MoS2-b3b4685fb6e1',
            'natoms': 3,
            'title': 'Monolayers (from C2DB)',
            'legend': 'Monolayers',
            'name': 'MoS2 (AB2-187-bi)',
            'label': 'MoS2 (AB2-187-bi)',
            'link': '/c2db/row/MoS2-b3b4685fb6e1',
            'method': 'DFT',
        }]}


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
    (dir / 'results-asr.structureinfo.json').write_text(
        '{"kwargs": {"data": {"has_inversion_symmetry": false}}}')
    (dir / 'results-asr.gs.json').write_text(
        '{"kwargs": {"data": {"gap": 1.8, "evac":4.5}}}')
    (dir / 'results-asr.gs@calculate.json').write_text(
        '{}')
    (dir / 'results-asr.convex_hull.json').write_text(
        json.dumps({'kwargs': {'data': CHULL}}))
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
    create_tree(Path())  # pragma: no cover
