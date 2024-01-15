"""Descriptions of CMR projects.

See https://cmr.fysik.dtu.dk/
TODO:

* c1db
* lowdim
* fix me units
* also add energy, fmax, smax, magmom
* imp2d tables
* bidb
* abs3: Electronic band structure
* extra???

"""
from __future__ import annotations
from typing import Callable
from cxdb.utils import Select, Input
from cxdb.material import Material
from cxdb.utils import FormPart

# Mapping from project name to ProjectDescription class:
projects: dict[str, Callable[[], ProjectDescription]] = {}


def project(name: str):
    """Decorator for filling in the _projects dict."""
    def decorator(cls):
        projects[name] = cls
        return cls
    return decorator


def create_project_description(name):
    """Create ProjectDescription for CMR projects."""
    if name in projects:
        return projects[name]()
    # Unknown CMR name.  Create a generic one:
    return ProjectDescription()


class ProjectDescription:
    title: str
    column_names: dict[str, str] = {}
    initial_columns: list[str] = ['formula', 'uid']
    uid: str = 'id'
    pbc: list[bool] | None = None
    extra: list[str] = []
    form_parts: list[FormPart] = []

    def postprocess(self, material: Material) -> None:
        pass


@project('solar')
class SolarProjectDescription(ProjectDescription):
    title = 'Organic Donor-Acceptor molecules'
    column_names = {
        'rho': 'Packing density',
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
        'V_oc': 'V_oc for PCBM acceptor [V]'}
    initial_columns = ['uid', 'formula', 'Unit', 'KS_gap', 'E_homo',
                       'E_lumo', 'E_gap', 'E_opt', 'rho']
    extra = ['CAN_SMILES', 'InChI', 'SMILES', 'Name', 'fold']


@project('adsorption')
class AdsorptionProjectDescription(ProjectDescription):
    title = 'Adsorption energy database'
    column_names = {
        'surf_mat': 'Surface Material',
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
        'RPA_adsorp': 'RPA correlation adsorption energy extrapolated [eV]'}
    initial_columns = ['adsorbate', 'surf_mat', 'LDA_adsorp', 'PBE_adsorp',
                       'RPBE_adsorp',
                       'BEEFvdW_adsorp', 'vdWDF2_adsorp', 'mBEEF_adsorp',
                       'mBEEFvdW_adsorp', 'RPA_EXX_adsorp']
    pbc = [True, True, False]  # overwrite pbc=(1,1,1) in db-file
    form_partsh = [
        Select('Surface material', 'surf_mat',
               [''] + ('Sc Ti V Cr Mn Fe Co Ni Cu '
                       'Y Zr Nb Mo Ru Rh Pd Ag '
                       'Hf Ta W Re Os Ir Pt Au').split()),
        Select('Adsorbate', 'adsorbate',
               [''] + 'H O N N2 CO NO CH OH'.split())]


@project('absorption_perovskites')
class AbsorptionPerovskitesProjectDescription(ProjectDescription):
    title = 'Absorption spectra of perovskites'
    column_names = {
        'name': 'Name of the structure',
        'phase': 'Perovskite phase',
        'gllbsc_dir_gap': 'Direct bandgap GLLB-SC [eV]',
        'gllbsc_ind_gap': 'Indirect bandgap GLLB-SC [eV]',
        'gllbsc_disc': 'Derivative discontinuity GLLB-SC [eV]',
        'eff': 'Calculated efficiency (thickness 1 nm) [eV]',
        'eff_max': 'Calculated efficiency (infinite thickness) [eV]'}
    initial_columns = ['uid', 'age', 'formula', 'energy', 'pbc',
                       'volume', 'charge',
                       'magmom', 'gllbsc_dir_gap', 'gllbsc_ind_gap', 'eff']


@project('abse3')
class ABSe3ProjectDescription(ProjectDescription):
    title = 'Ternary Selenides ABSe3'
    column_names = {
        'prototype': 'name of the prototype',
        'name': 'name of the system',
        'space_group': 'Space group',
        'pbe_tot_en': 'PBE total energy [eV]',
        'pbe_hof': 'PBE heat of formation [eV/atom]'}
    initial_columns = ['prototype', 'name', 'space_group', 'pbe_tot_en',
                       'pbe_hof']


@project('abs3')
class ABS3ProjectDescription(ProjectDescription):
    title = 'Database of ABS3 materials'
    column_names = {
        'E_hull': 'E-hull [eV]',
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
        'isreference': 'Is reference'}
    initial_columns = ['formula', 'energy',
                       'ABS3_name', 'prototype', 'GLLB_ind', 'GLLB_dir',
                       'E_relative_per_atom']
    form_parts = [
        Select('Prototype name', 'prototype',
               ',NH4CdCl3/Sn2S3,GdFeO3,YScS3,PbPS3,'
               'BaNiO3,FePS3,cubic,distorted,Pyroxene-CaIrO3,BiInS3,'
               'CuTaS3,CeTmS3'.split(','))]


@project('abx2')
class ABX2ProjectDescription(ProjectDescription):
    title = 'Database of ABX2 materials'
    column_names = {
        'E_hull': 'E-hull [eV]',
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
        'total_en': 'Total energy  (mBEEF) [eV]'}
    initial_columns = ['short_name', 'phase', 'GLLB_ind', 'GLLB_dir',
                       'TB09_ind', 'E_relative_perAtom']
    form_parts = [
        Select(
            'Phase: ZB=zinckblende, WZ=wurtzite, KT=wz kesterite, '
            'ST=wz stannite, CP=zb chaclopyrite', 'phase',
            ',ST,KT,CP,WZ,ZB'.split(','))]


@project('agau309')
class AgAu309ProjectDescription(ProjectDescription):
    title = 'AgAu309'
    initial_columns = ['calculator', 'formula', 'energy', 'mass']


@project('a2bcx4')
class A2BCX4ProjectDescription(ProjectDescription):
    title = 'Database of A2BCX4 materials'
    column_names = {
        'E_relative_per_atom': 'Energy per atom [meV]',
        'E_uncertainty_per_atom': 'Uncertainty of the total energy [meV]',
        'GLLB_dir': 'Direct band gap (GLLB-SC) [eV]',
        'GLLB_ind': 'Indirect band gap (GLLB-SC) [eV]',
        'lattice': 'Crystal system',
        'prototype': 'Prototype',
        'space_group': 'Space group: I-4, I4-2m, P31, Ama2, P1n1, Pmn21',
        'name': 'Short chemical formula',
        'mbeef_en': 'Total energy  (mBEEF) [eV]'}
    initial_columns = ['name', 'space_group', 'GLLB_ind', 'GLLB_dir',
                       'E_relative_per_atom', 'E_uncertainty_per_atom']
    form_parts = [
        Select(
            'Template crystal structure or prototype present '
            'in the ICSD database',
            'prototype',
            ',Cu2ZnSiS4,Cu2ZnSnS4,Cu2CdSnS4,Cu2CdGeS4,'
            'Cu2KVS4,Cu2BaGeSe4'.split(','))]


@project('catapp')
class CatAppProjectDescription(ProjectDescription):
    title = 'CatApp database'
    column_names = {
        'a': 'Reactant A',
        'b': 'Reactant B',
        'ab': 'Product AB',
        'surface': 'Description of surface',
        'facet': 'Surface facet',
        'site': 'Adsorption site',
        'xc': 'XC-functional',
        'reference': 'Reference',
        'url': 'Link to reference',
        'dataset': 'Description of calculation',
        'ref': 'Ref',
        'ea': 'Ea',
        'er': 'Er'}
    initial_columns = ['id', 'formula', 'surface', 'a',
                       'b', 'ab', 'facet', 'ea',
                       'er', 'xc']


@project('cubic_perovskites')
class CubicPerovskitesProjectDescription(ProjectDescription):
    title = 'Perovskite water-splitting'
    column_names = {
        'A_ion': 'A-ion in the cubic perovskite',
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
        'standard_energy': 'Standard energy'}
    initial_columns = [
        'uid', 'age', 'formula', 'energy', 'pbc', 'volume', 'charge',
        'gllbsc_dir_gap', 'gllbsc_ind_gap', 'heat_of_formation_all']


@project('dssc')
class DSSCProjectDescription(ProjectDescription):
    title = 'Porphyrin based dyes'
    column_names = {
        'M': 'Metal center',
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
        'UsedE1': 'UsedE1',
        'UsedEc': 'UsedEc',
        'dssc': 'DSSC',
        'gpaw_setups': 'GPAW-potential',
        'gpaw_version': 'GPAW-version',
        'xc': 'XC',
        'DBScreen': 'DB-screen'}
    initial_columns = [
        'uid', 'age', 'formula', 'calculator', 'energy', 'fmax',
        'charge', 'mass', 'magmom', 'E_gap', 'KS_gap']


@project('funct_perovskites')
class FunctPerovskitesProjectDescription(ProjectDescription):
    title = 'Functional Perovskites'
    column_names = {
        'gllbsc_dir_gap': 'Direct bandgap GLLB-SC [eV]',
        'gllbsc_ind_gap': 'Indirect bandgap GLLB-SC [eV]',
        'gllbsc_disc': 'Derivative discontinuity GLLB-SC [eV]',
        'gllbsc_gamma': 'Bandgap at Gamma GLLB-SC [eV]',
        'comb_A': 'Cubic perovskite labeled as A-combination',
        'comb_B': 'Cubic perovskite labeled as B-combination',
        'sequence': 'Sequence of A and B layers',
        'gllbsc_gamma_gap': 'GLLBSC Gamma-gap'}
    initial_columns = [
        'uid', 'age', 'formula', 'energy', 'pbc', 'volume',
        'charge', 'mass', 'comb_A', 'comb_B', 'sequence',
        'gllbsc_gamma_gap']


@project('low_symmetry_perovskites')
class LowSymmetryPerovskitesProjectDescription(ProjectDescription):
    title = 'Low symmetry perovskites'
    column_names = {
        'gllbsc_dir_gap': 'Direct bandgap GLLB-SC [eV]',
        'gllbsc_ind_gap': 'Indirect bandgap GLLB-SC [eV]',
        'gllbsc_disc': 'Derivative discontinuity GLLB-SC [eV]',
        'phase': 'double ruddlesden-popper or dion-jacobsen',
        'ruddlesden_popper': 'Ruddlesden popper phase A2BO4 A3B2O7 or A2BO3N',
        'dion_jacobsen': '1: A=Na Mg K Rb or 2: A=Rb Sr Cs Ba',
        'composition': 'Composition'}
    initial_columns = [
        'uid', 'age', 'formula', 'energy', 'pbc', 'volume', 'charge',
        'gllbsc_dir_gap', 'gllbsc_ind_gap', 'phase']


@project('mp_gllbsc')
class MPGLLBSCProjectDescription(ProjectDescription):
    title = 'New Light Harvesting Materials'
    column_names = {
        'gllbsc_dir_gap': 'Direct bandgap GLLB-SC. [eV]',
        'gllbsc_ind_gap': 'Indirect bandgap GLLB-SC. [eV]',
        'gllbsc_disc': 'Derivative discontinuity GLLB-SC. [eV]',
        'mpid': 'ID of materials in Materials project',
        'icsd_id': 'ID of materials in ICSD',
        'g0w0_gap': '`G_0W_0` gap at `\\Gamma` [eV]',
        'gw0_gap': '`GW_0` gap at `\\Gamma` [eV]',
        'gw_gap': '`GW` gap at `\\Gamma` [eV]',
        'hse06_gap': 'HSE06 gap at `\\Gamma` [eV]',
        'lda_gap': 'LDA gap at `\\Gamma` [eV]',
        'gllbsc_gap': 'GLLBSC gap at `\\Gamma` [eV]'}
    initial_columns = [
        'uid', 'age', 'formula', 'energy', 'pbc',
        'volume', 'charge', 'magmom',
        'gllbsc_dir_gap', 'gllbsc_ind_gap', 'mpid']

    def postprocess(self, material: Material):
        material.add_column('mpid', str(material.mpid), update=True)
        material.add_column('icsd_id', str(material.icsd_id), update=True)


@project('oqmd123')
class OQMD123ProjectDescription(ProjectDescription):
    title = 'One, two and three component references from OQMD'
    column_names = {
        'hform': 'Heat of formation [eV/atom]',
        'oqmd_id': 'OQMD ID'}
    initial_columns = [
        'formula', 'energy', 'hform', 'stoichiometry', 'volume',
        'fmax', 'smax']
    uid = 'uid'

    def postprocess(self, material: Material):
        id = material.get('oqmd_id')
        material.add_column(
            'oqmd_id',
            id,
            f'<a href="http://oqmd.org/materials/entry/{id}">{id}</a>',
            update=True)


@project('organometal')
class OrganometalProjectDescription(ProjectDescription):
    title = 'Organometal Halide Perovskites'
    column_names = {
        'gllbsc_dir_gap': 'Direct GLLB-SC+SOC bandgap [eV]',
        'gllbsc_ind_gap': 'Indirect GLLB-SC+SOC bandgap [eV]',
        'gllbsc_disc': 'GLLB-SC Derivative discontinuity [eV]',
        'name': 'Name given to the crystal structure',
        'symmetry': 'Symmetry of the crystal',
        'space_group': 'Space group of the crystal'}
    initial_columns = [
        'uid', 'age', 'formula', 'energy', 'pbc', 'volume',
        'charge', 'gllbsc_dir_gap', 'gllbsc_ind_gap']


@project('pv_pec_oqmd')
class PVPECOQMDProjectDescription(ProjectDescription):
    title = 'Screening for PV and PEC materials using the OQMD database'
    column_names = {
        'm_e': 'Effective electron mass [`m_e`]',
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
        'magnetic': 'Magnetic'}
    initial_columns = [
        'formula', 'lattice', 'spacegroup',
        'GLLB_ind', 'GLLB_dir', 'Dxc', 'PBE_gap',
        'defect_tolerant', 'magnetic', 'icsd', 'm_e', 'm_h']

    def postprocess(self, material: Material):
        material.add_column('icsd', str(material.icsd), update=True)


@project('imp2d')
class Imp2DProjectDescription(ProjectDescription):
    title = 'Impurities in 2D Materials Database'
    column_names = {
        'host': 'Host chemical formula',
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
        'hostenergy': 'Host energy [eV]'}
    initial_columns = [
        'host', 'dopant', 'defecttype', 'spin']
    uid = 'name'
    pbc = [True, True, False]  # overwrite pbc=(1,1,1) in db-file
    form_parts = [
        Input('Host formiula', 'host', 'e.g. MoS2, WS2, TiO2, etc.'),
        Input('Dopant species', 'dopant', 'e.g. H, He, Li, etc.'),
        Select('Defecttype', 'defecttype',
               ['', 'interstitial', 'adsorbate']),
        Select('Relaxation fully converged?', 'converged',
               ['', 'True', 'False']),
        Input('Host spacegroup number', 'host_spacegroup', 'e.g. 187')]


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


@project('ads1d')
class Ads1DProjectDescription(ProjectDescription):
    title = 'Chemisorption of gas atoms on 1D transition-metal halides'
    column_names = {
        'x': 'X',
        'y': 'Y',
        'n': 'XY3 formula units',
        'adsorbate': 'Adsorbate',
        'workfunction': 'Workfunction [V]',
        'eads': 'Adsorption energy [eV]',
        'coverage': 'Coverage'}
    initial_columns = ['formula',
                       'x', 'y', 'n', 'adsorbate', 'workfunction',
                       'coverage', 'eads']
    uid = 'uid'


if __name__ == '__main__':
    # Convert cmr.<proj>.custum.Template.raw_key_value_descriptions:
    for k, v in projects.items():
        print(k)
        for n, x in v().column_names.items():
            if isinstance(x, str):
                break
            i, j, u = x
            j = j or i
            if u:
                j += f' [{u}]'
            print(f'        {n!r}: {j!r},')

    print(' '.join(f'{x}.db' for x in projects))
