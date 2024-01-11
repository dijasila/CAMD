"""Descriptions of CMR projects.

See https://cmr.fysik.dtu.dk/
"""
from __future__ import annotations
from typing import Callable
from cxdb.utils import Select

_projects: dict[str, Callable[[], ProjectDescription]] = {}


def project(func):
    """Decorator."""
    _projects[func.__name__] = func
    return func


def create_project_description(name):
    """Create ProjectDescription for CMR projects."""
    if name in _projects:
        return _projects[name]()
    # Unknown CMR name.  Create a generic one:
    return ProjectDescription(name, {}, ['formula', 'uid'])


class ProjectDescription:
    def __init__(self,
                 title: str,
                 column_names: dict[str, str],
                 initial_columns: list[str],
                 uid: str = 'id',
                 ndims: int | None = None,
                 pbc: list[bool] | None = None,
                 extra: list[str] | None = None,
                 search: list | None = None):
        self.title = title
        self.column_names = column_names
        self.initial_columns = initial_columns
        self.uid = uid
        self.ndims = ndims
        self.pbc = pbc
        self.extra = extra or []
        self.search = search or []


@project
def solar():
    return ProjectDescription(
        'Organic Donor-Acceptor molecules',
        {'rho': 'Packing density',
         'dip': 'Dipole moment [Debye]',
         'KS_gap': 'HOMO-LUMO gap (B3LYP) [eV]',
         'E_homo': 'HOMO (B3LYP) [eV]',
         'E_lumo': 'LUMO (B3LYP) [eV]',
         'E_gap': 'Tight-binding gap extrapolated [eV]',
         'E_opt': 'Optical gap (singlet-triplet gap) [eV]',
         'Unit': 'Unit',
         'E_homo_TB': 'HOMO (Tight-binding) [eV]',
         'E_lumo_TB': 'LUMO (Tight-binding) [eV]',
         'D_homo': 'Dimer-monomer HOMO difference (B3LYP) [eV]',
         'D_lumo': 'Dimer-monomer LUMO difference (B3LYP) [eV]',
         'DeltaU': 'DeltaU [eV]',
         'Energy': 'GS Energy [eV]',
         'V_oc': 'V_oc for PCBM acceptor [V]'},
        ['uid', 'formula', 'Unit', 'KS_gap', 'E_homo', 'E_lumo',
         'E_gap', 'E_opt', 'rho'],
        extra=['CAN_SMILES', 'InChI', 'SMILES', 'Name', 'fold'],
        ndims=0)


@project
def adsorption():
    return ProjectDescription(
        'Adsorption energy database',
        {'surf_mat': 'Surface Material',
         'adsorbate': 'Adsorbate',
         'mol': 'Reference molecule 1',
         'mol2': 'Reference molecule 2',
         'LDA_adsorp': 'Adsorption energy with LDA [eV]',
         'PBE_adsorp': 'Adsorption energy with PBE [eV]',
         'RPBE_adsorp': 'Adsorption energy with RPBE [eV]',
         'BEEFvdW_adsorp': 'Adsorption energy with BEEF-vdW [eV]',
         'vdWDF2_adsorp': 'Adsorption energy with vdW-DF2 [eV]',
         'mBEEF_adsorp': 'Adsorption energy with mBEEF [eV]',
         'mBEEFvdW_adsorp': 'Adsorption energy with mBEEF-vdW [eV]',
         'EXX_adsorp': 'Adsorption energy with EXX [eV]',
         'RPA_adsorp': 'RPA correlation adsorption energy extrapolated [eV]'},
        ['adsorbate', 'surf_mat', 'LDA_adsorp', 'PBE_adsorp',
         'RPBE_adsorp',
         'BEEFvdW_adsorp', 'vdWDF2_adsorp', 'mBEEF_adsorp',
         'mBEEFvdW_adsorp', 'RPA_EXX_adsorp'],
        ndims=2,
        pbc=[True, True, False],  # overwrite pbc=(1,1,1) in db-file
        search=[Select('Surface material', 'surf_mat',
                       ('Sc Ti V Cr Mn Fe Co Ni Cu '
                        'Y Zr Nb Mo Ru Rh Pd Ag '
                        'Hf Ta W Re Os Ir Pt Au').split()),
                Select('Adsorbate', 'adsorbate',
                       'H O N N2 CO NO CH OH'.split())])


if __name__ == '__main__':
    for k, v in _projects.items():
        for n, x in v().column_names.items():
            if isinstance(x, str):
                break
            i, j, u = x
            j = j or i
            if u:
                j += f' [{u}]'
            print(f'        {n!r}: {j!r},')
