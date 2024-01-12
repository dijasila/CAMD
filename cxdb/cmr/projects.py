"""Descriptions of CMR projects.

See https://cmr.fysik.dtu.dk/
TODO:

* c1db
* lowdim
* fix me units
* also add energy, fmax, smax, magmom
* imp2d tables
* bidb

"""
from __future__ import annotations
from typing import Callable
from cxdb.utils import Select, Input
from cxdb.material import Material

# Mapping from project name to ProjectDescription factory:
_projects: dict[str, Callable[[], ProjectDescription]] = {}


def project(func):
    """Decorator for filling in the _projects dict."""
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
                 pbc: list[bool] | None = None,
                 extra: list[str] | None = None,
                 search: list | None = None,
                 postprocess: Callable[[Material], None] | None = None):
        self.title = title
        self.column_names = {name: val or name.title()
                             for name, val in column_names.items()}
        self.initial_columns = initial_columns
        self.uid = uid
        self.pbc = pbc
        self.extra = extra or []
        self.search = search or []
        self.postprocess = postprocess or (lambda mat: None)


@project
def solar() -> ProjectDescription:
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
        extra=['CAN_SMILES', 'InChI', 'SMILES', 'Name', 'fold'])


@project
def adsorption() -> ProjectDescription:
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
        pbc=[True, True, False],  # overwrite pbc=(1,1,1) in db-file
        search=[Select('Surface material', 'surf_mat',
                       [''] + ('Sc Ti V Cr Mn Fe Co Ni Cu '
                               'Y Zr Nb Mo Ru Rh Pd Ag '
                               'Hf Ta W Re Os Ir Pt Au').split()),
                Select('Adsorbate', 'adsorbate',
                       [''] + 'H O N N2 CO NO CH OH'.split())])


@project
def absorption_perovskites() -> ProjectDescription:
    return ProjectDescription(
        'Absorption spectra of perovskites',
        {'name': 'Name of the structure',
         'phase': 'Perovskite phase',
         'gllbsc_dir_gap': 'Direct bandgap GLLB-SC [eV]',
         'gllbsc_ind_gap': 'Indirect bandgap GLLB-SC [eV]',
         'gllbsc_disc': 'Derivative discontinuity GLLB-SC [eV]',
         'eff': 'Calculated efficiency (thickness 1 nm) [eV]',
         'eff_max': 'Calculated efficiency (infinite thickness) [eV]',
         'project': ''},
        ['uid', 'age', 'formula', 'energy', 'pbc', 'volume', 'charge',
         'magmom', 'gllbsc_dir_gap', 'gllbsc_ind_gap', 'eff'])


@project
def abse3() -> ProjectDescription:
    return ProjectDescription(
        'Ternary Selenides ABSe3',
        dict(
            prototype='name of the prototype',
            name='name of the system',
            space_group='Space group',
            pbe_tot_en='PBE total energy [eV]',
            pbe_hof='PBE heat of formation [eV/atom]'),
        ['prototype', 'name', 'space_group', 'pbe_tot_en', 'pbe_hof'])


@project
def abs3() -> ProjectDescription:
    return ProjectDescription(
        'Database of ABS3 materials',
        {'E_hull': 'E-hull [eV]',
         'm_e': 'Effective electron mass [`m_e`]',
         'm_h': 'Effective hole mass [`m_e`]',
         'E_relative_per_atom': 'Energy per atom [eV]',
         'E_uncertainty_hull': 'E-hull Uncertainty [eV]',
         'E_uncertainty_per_atom': 'Energy uncertainty [eV]',
         'GLLB_dir': 'Direct band gap (GLLB-SC) [eV]',
         'GLLB_ind': 'Indirect band gap (GLLB-SC) [eV]',
         'PBEsol_gap': 'Band gap (PBEsol) [eV]',
         'lattice': 'Crystal system',
         'prototype': 'prototype name',
         'ABS3_name': 'Short chemical formula',
         'isreference': 'Is reference'},
        ['formula', 'energy',
         'ABS3_name', 'prototype', 'GLLB_ind', 'GLLB_dir',
         'E_relative_per_atom'],
        search=[
            Select('Prototype name', 'prototype',
                   ',NH4CdCl3/Sn2S3,GdFeO3,YScS3,PbPS3,'
                   'BaNiO3,FePS3,cubic,distorted,Pyroxene-CaIrO3,BiInS3,'
                   'CuTaS3,CeTmS3'.split(','))])
    # XXX:
    # GLLB-SC energy relative to VBM [eV]
    # Electronic band structure


@project
def abx2() -> ProjectDescription:
    return ProjectDescription(
        'Database of ABX2 materials',
        {'E_hull': 'E-hull [eV]',
         'KS_gap': 'Kohn Sham band gap [eV]',
         'm_e': 'Effective electron mass [`m_e`]',
         'm_h': 'Effective hole mass [`m_e`]',
         'Dxc': 'Derivative discontinuity (GLLB-SC) [eV]',
         'E_relative_perAtom': 'Energy per atom [eV]',
         'E_uncertanty_hull': 'Uncertanty of the convex hull energy [eV]',
         'E_uncertanty_perAtom': 'Uncertanty of the total energy [eV]',
         'GLLB_dir': 'Direct band gap (GLLB-SC) [eV]',
         'GLLB_ind': 'Indirect band gap (GLLB-SC) [eV]',
         'TB09_dir': 'Direct band gap (TB09) [eV]',
         'TB09_ind': 'Indirect band gap (TB09) [eV]',
         'lattice': 'Crystal system',
         'phase': 'Phase',
         'short_name': 'Short chemical formula',
         'total_en': 'Total energy  (mBEEF) [eV]'},
        ['short_name', 'phase', 'GLLB_ind', 'GLLB_dir',
         'TB09_ind', 'E_relative_perAtom'],
        search=[
            Select(
                'Phase: ZB=zinckblende, WZ=wurtzite, KT=wz kesterite, '
                'ST=wz stannite, CP=zb chaclopyrite', 'phase',
                ',ST,KT,CP,WZ,ZB'.split(','))])


@project
def agau309() -> ProjectDescription:
    return ProjectDescription(
        'AgAu309',
        {},
        ['calculator', 'formula', 'energy', 'mass'])


@project
def a2bcx4() -> ProjectDescription:
    return ProjectDescription(
        'Database of A2BCX4 materials',
        {'E_relative_per_atom': 'Energy per atom [meV]',
         'E_uncertainty_per_atom': 'Uncertainty of the total energy [meV]',
         'GLLB_dir': 'Direct band gap (GLLB-SC) [eV]',
         'GLLB_ind': 'Indirect band gap (GLLB-SC) [eV]',
         'lattice': 'Crystal system',
         'prototype': 'Prototype',
         'space_group': 'Space group: I-4, I4-2m, P31, Ama2, P1n1, Pmn21',
         'name': 'Short chemical formula',
         'mbeef_en': 'Total energy  (mBEEF) [eV]'},
        ['name', 'space_group', 'GLLB_ind', 'GLLB_dir',
         'E_relative_per_atom', 'E_uncertainty_per_atom'],
        search=[
            Select(
                'Template crystal structure or prototype present '
                'in the ICSD database',
                'prototype',
                ',Cu2ZnSiS4,Cu2ZnSnS4,Cu2CdSnS4,Cu2CdGeS4,'
                'Cu2KVS4,Cu2BaGeSe4'.split(','))])


@project
def catapp() -> ProjectDescription:
    return ProjectDescription(
        'CatApp database',
        {'a': 'Reactant A',
         'b': 'Reactant B',
         'ab': 'Product AB',
         'surface': 'Description of surface',
         'facet': 'Surface facet',
         'site': 'Adsorption site',
         'xc': 'XC-functional',
         'reference': 'Reference',
         'url': 'Link to reference',
         'dataset': 'Description of calculation',
         'project': '',
         'ref': '',
         'ea': '',
         'er': ''},
        ['id', 'formula', 'surface', 'a',
         'b', 'ab', 'facet', 'ea',
         'er', 'xc'])


@project
def cubic_perovskites() -> ProjectDescription:
    return ProjectDescription(
        'Perovskite water-splitting',
        {'A_ion': 'A-ion in the cubic perovskite',
         'B_ion': 'B-ion in the cubic perovskite',
         'anion': 'Anion combination in the perovskite',
         'gllbsc_dir_gap': 'Direct bandgap GLLB-SC [eV]',
         'gllbsc_ind_gap': 'Indirect bandgap GLLB-SC [eV]',
         'heat_of_formation_all': 'Heat of formation [eV]',
         'combination': 'General formula',
         'CB_dir': 'Direct position of the conduction band edge',
         'CB_ind': 'Indirect position of the conduction band edge',
         'VB_dir': 'Direct position of the valence band edge',
         'VB_ind': 'Indirect position of the valence band edge',
         'reference': 'Reference',
         'standard_energy': '',
         'project': ''},
        ['uid', 'age', 'formula', 'energy', 'pbc', 'volume', 'charge',
         'gllbsc_dir_gap', 'gllbsc_ind_gap', 'heat_of_formation_all'])


@project
def dssc() -> ProjectDescription:
    return ProjectDescription(
        'Porphyrin based dyes',
        {'M': 'Metal center',
         'A': 'Anchor group',
         'R1': 'First side group',
         'R2': 'Second side group',
         'R3': 'Third side group',
         'KS_HOMO': 'Kohn-Sham HOMO eigenvalue [eV]',
         'KS_LUMO': 'Kohn-Sham LUMO eigenvalue [eV]',
         'KS_gap': 'KS gap (KS_LUMO - KS_HOMO) [eV]',
         'E_HOMO': 'IP [eV]',
         'E_LUMO': 'EA [eV]',
         'E_gap': 'EA - IP [eV]',
         'E_c': 'Cond. band - IP [eV]',
         'E_1': 'Triplet optical gap [eV]',
         'E_opt_LUMO': 'IP + troiplet [eV]',
         'LQual1': 'Level align. 1 [eV]',
         'LQual2': 'Level align. 2 [eV]',
         'UsedE1': '',
         'UsedEc': '',
         'dssc': '',
         'gpaw_setups': '',
         'gpaw_version': '',
         'project': '',
         'xc': '',
         'DBScreen': ''},
        ['uid', 'age', 'formula', 'calculator', 'energy', 'fmax',
         'charge', 'mass', 'magmom', 'E_gap', 'KS_gap'])


@project
def funct_perovskites() -> ProjectDescription:
    return ProjectDescription(
        'Functional Perovskites',
        {'gllbsc_dir_gap': 'Direct bandgap GLLB-SC [eV]',
         'gllbsc_ind_gap': 'Indirect bandgap GLLB-SC [eV]',
         'gllbsc_disc': 'Derivative discontinuity GLLB-SC [eV]',
         'gllbsc_gamma': 'Bandgap at Gamma GLLB-SC [eV]',
         'comb_A': 'Cubic perovskite labeled as A-combination',
         'comb_B': 'Cubic perovskite labeled as B-combination',
         'sequence': 'Sequence of A and B layers',
         'project': '',
         'gllbsc_gamma_gap': ''},
        ['uid', 'age', 'formula', 'energy', 'pbc', 'volume',
         'charge', 'mass', 'comb_A', 'comb_B', 'sequence',
         'gllbsc_gamma_gap'])


@project
def low_symmetry_perovskites() -> ProjectDescription:
    return ProjectDescription(
        'Low symmetry perovskites',
        {'gllbsc_dir_gap': 'Direct bandgap GLLB-SC [eV]',
         'gllbsc_ind_gap': 'Indirect bandgap GLLB-SC [eV]',
         'gllbsc_disc': 'Derivative discontinuity GLLB-SC [eV]',
         'phase': 'double ruddlesden-popper or dion-jacobsen',
         'ruddlesden_popper': 'Ruddlesden popper phase A2BO4 A3B2O7 or A2BO3N',
         'dion_jacobsen': '1: A=Na Mg K Rb or 2: A=Rb Sr Cs Ba',
         'project': '',
         'composition': ''},
        ['uid', 'age', 'formula', 'energy', 'pbc', 'volume', 'charge',
         'gllbsc_dir_gap', 'gllbsc_ind_gap', 'phase'])


@project
def mp_gllbsc() -> ProjectDescription:
    return ProjectDescription(
        'New Light Harvesting Materials',
        {'gllbsc_dir_gap': 'Direct bandgap GLLB-SC. [eV]',
         'gllbsc_ind_gap': 'Indirect bandgap GLLB-SC. [eV]',
         'gllbsc_disc': 'Derivative discontinuity GLLB-SC. [eV]',
         'mpid': 'ID of materials in Materials project',
         'icsd_id': 'ID of materials in ICSD',
         'g0w0_gap': '`G_0W_0` gap at `\\Gamma` [eV]',
         'gw0_gap': '`GW_0` gap at `\\Gamma` [eV]',
         'gw_gap': '`GW` gap at `\\Gamma` [eV]',
         'hse06_gap': 'HSE06 gap at `\\Gamma` [eV]',
         'lda_gap': 'LDA gap at `\\Gamma` [eV]',
         'gllbsc_gap': 'GLLBSC gap at `\\Gamma` [eV]',
         'project': ''},
        ['uid', 'age', 'formula', 'energy', 'pbc',
         'volume', 'charge', 'magmom',
         'gllbsc_dir_gap', 'gllbsc_ind_gap', 'mpid'],
        postprocess=pp_mp_gllbsc)


def pp_mp_gllbsc(material: Material):
    material.add_column('mpid', str(material.mpid), update=True)
    material.add_column('icsd_id', str(material.icsd_id), update=True)


@project
def oqmd123() -> ProjectDescription:
    return ProjectDescription(
        'One, two and three component references from OQMD',
        {'hform': 'Heat of formation [eV/atom]',
         'oqmd_id': 'OQMD ID'},
        ['formula', 'energy', 'hform', 'stoichiometry', 'volume',
         'fmax', 'smax'],
        uid='uid',
        postprocess=pp_oqmd123)


def pp_oqmd123(material: Material):
    id = material.get('oqmd_id')
    if id:
        material.add_column(
            'oqmd_id',
            id,
            f'<a href="http://oqmd.org/materials/entry/{id}">{id}</a>',
            update=True)


@project
def organometal() -> ProjectDescription:
    return ProjectDescription(
        'Organometal Halide Perovskites',
        {'gllbsc_dir_gap': 'Direct GLLB-SC+SOC bandgap [eV]',
         'gllbsc_ind_gap': 'Indirect GLLB-SC+SOC bandgap [eV]',
         'gllbsc_disc': 'GLLB-SC Derivative discontinuity [eV]',
         'name': 'Name given to the crystal structure',
         'symmetry': 'Symmetry of the crystal',
         'space_group': 'Space group of the crystal',
         'project': ''},
        ['uid', 'age', 'formula', 'energy', 'pbc', 'volume',
         'charge', 'gllbsc_dir_gap', 'gllbsc_ind_gap'])


@project
def pv_pec_oqmd() -> ProjectDescription:
    return ProjectDescription(
        'Screening for PV and PEC materials using the OQMD database',
        {'m_e': 'Effective electron mass [`m_e`]',
         'm_h': 'Effective hole mass [`m_e`]',
         'GLLB_dir': 'Direct band gap (GLLB-SC) [eV]',
         'GLLB_ind': 'Indirect band gap (GLLB-SC) [eV]',
         'PBE_gap': 'Band gap (PBE) [eV]',
         'lattice': 'Crystal system',
         'icsd': 'ICSD number',
         'Dxc': 'Derivative discontinuity (GLLB-SC) [eV]',
         'ntypes': 'Number of species',
         'spacegroup': 'Space group',
         'defect_tolerant': 'Defect tolerant',
         'magnetic': 'Magnetic'},
        ['formula', 'lattice', 'spacegroup',
         'GLLB_ind', 'GLLB_dir', 'Dxc', 'PBE_gap',
         'defect_tolerant', 'magnetic', 'icsd', 'm_e', 'm_h'],
        postprocess=pp_pv_pec_oqmd)


def pp_pv_pec_oqmd(material: Material):
    material.add_column('icsd', str(material.icsd), update=True)


@project
def imp2d() -> ProjectDescription:
    return ProjectDescription(
        'Impurities in 2D Materials Database',
        {'host': 'Host chemical formula',
         'dopant': 'Impurity species',
         'defecttype': 'Impurity type',
         'eform': 'Formation energy [eV]',
         'site': 'Defect setup position',
         'depth': 'Depth parameter',
         'extension_factor': 'Extension factor',
         'en2': 'Total energy [eV]',
         'conv2': 'Total energy convergence [eV]',
         'spin': 'Magnetic moment',
         'supercell': 'Calculation supercell',
         'host_spacegroup': 'Host spacegroup symbol',
         'converged': 'Relaxation converged?',
         'name': 'System name',
         'en1': 'Total energy, first WF step [eV]',
         'conv1': 'Total energy convergence, first WF step [eV]',
         'dopant_chemical_potential': 'Dopant chemical potential [eV]',
         'hostenergy': 'Host energy [eV]'},
        ['host', 'dopant', 'defecttype', 'spin'],
        uid='name',
        pbc=[True, True, False],  # overwrite pbc=(1,1,1) in db-file
        search=[
            Input('Host formiula', 'host', 'e.g. MoS2, WS2, TiO2, etc.'),
            Input('Dopant species', 'dopant', 'e.g. H, He, Li, etc.'),
            Select('Defecttype', 'defecttype',
                   ['', 'interstitial', 'adsorbate']),
            Select('Relaxation fully converged?', 'converged',
                   ['', 'True', 'False']),
            Input('Host spacegroup number', 'host_spacegroup', 'e.g. 187')])


tables = [
    (['Host properties', ''],
     ['host',
      'supercell',
      'host_spacegroup',
      'hostenergy']),
    (['Structural properties', ''],
     ['dopant',
      'defecttype',
      'depth',
      'extension_factor',
      'site']),
    (['Electronic properties', ''],
     ['spin',
      'en2',
      'eform',
      'dopant_chemical_potential',
      'conv2'])]


@project
def ads1d() -> ProjectDescription:
    return ProjectDescription(
        'Chemisorption of gas atoms on 1D transition-metal halides',
        {'x': 'X',
         'y': 'Y',
         'n': 'XY3 formula units',
         'adsorbate': 'Adsorbate',
         'workfunction': 'Workfunction [V]',
         'eads': 'Adsorption energy [eV]',
         'coverage': 'Coverage'},
        ['formula',
         'x', 'y', 'n', 'adsorbate', 'workfunction',
         'coverage', 'eads'],
        uid='uid')


if __name__ == '__main__':
    # Convert cmr.<proj>.custum.Template.raw_key_value_descriptions:
    for k, v in _projects.items():
        print(k)
        for n, x in v().column_names.items():
            if isinstance(x, str):
                break
            i, j, u = x
            j = j or i
            if u:
                j += f' [{u}]'
            print(f'        {n!r}: {j!r},')

    print(' '.join(f'{x}.db' for x in _projects))

    downloads = [
        # ('adsorption', ['adsorption.db', 'surfaces.db']),
        ('absorption_perovskites', ['absorption_perovskites.db',
                                    'perfect_abs.txt']),
        ('abse3', ['abse3.db']),
        ('abs3', ['abs3.db']),
        ('abx2', ['abx2.db']),
        ('agau309', ['agau309.db']),
        ('a2bcx4', ['a2bcx4.db']),
        # ('beef', ['molecules.db', 'solids.db']),
        # ('bondmin', ['bondmin.db']),
        ('catapp', ['catapp.db', 'catappdata.csv']),
        # ('compression', ['compression.db']),
        ('cubic_perovskites', ['cubic_perovskites.db']),
        # ('c1db', ['c1db.db']),
        # ('c2db', ['c2db.db', 'workflow.png', 'c2db-db.png']),
        # ('c2dm', ['c2dm.db']),
        # ('dcdft', ['dcdft.db', 'dcdft_gpaw_pw_paw09.db']),
        ('dssc', ['dssc.db']),
        # ('fcc111', ['fcc111.db']),
        ('funct_perovskites', ['funct_perovskites.db']),
        # ('gbrv', ['gbrv.db']),
        # ('g2', ['g2.db']),
        # ('htgw', ['htgw.db']),
        ('lowdim', ['lowdim.db']),
        ('low_symmetry_perovskites', ['low_symmetry_perovskites.db']),
        ('mp_gllbsc', ['mp_gllbsc.db']),
        # ('oqmd12', ['oqmd12.db']),
        ('oqmd123', ['oqmd123.db']),
        ('organometal', ['organometal.db']),
        ('pv_pec_oqmd', ['pv_pec_oqmd.db']),
        # ('qpod', ['qpod.db', 'workflow_qpod.png']),
        # ('solar', ['solar.db']),
        # ('tmfp06d', ['tmfp06d.db']),
        ('imp2d', ['imp2d.db']),
        ('bidb', ['bidb.db', 'bidb-workflow.png']),
        # ('vdwh', ['chi-data.tar.gz', 'graphene-data.tar.gz',
        #           'chi-data-v2.tar.gz'])
        ('ads1d', ['ads1d.db'])]

    # skip = {
    #     'dcdft', 'beef', 'vdwh', 'c2dm', 'gbrv', 'htgw', 'bondmin', 'qpod'}

    from pathlib import Path
    for name, _ in downloads:
        continue
        p = Path(f'/home/jensj/cmr/cmr/{name}/custom.py')
        print(f"""
@project
def {name}() -> ProjectDescription:
    return ProjectDescription(""", p.read_text())
