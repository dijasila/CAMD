"""Descriptions of CMR projects.

See https://cmr.fysik.dtu.dk/
"""
from __future__ import annotations
from typing import Callable

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
                 extra: list[str] | None = None):
        self.title = title
        self.column_names = column_names
        self.initial_columns = initial_columns
        self.uid = uid
        self.ndims = ndims
        self.extra = extra or []


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
