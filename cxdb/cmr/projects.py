"""Descriptions of CMR projects.

See https://cmr.fysik.dtu.dk/

TODO:

* c1db
* lowdim
* bidb: magnetic: yes or 1. slide_stability: Stable?
* extra???

"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from ase.db import connect
from cxdb.cmr.lowdim import LowDimPanel, keysfortable0
from cxdb.material import Material, Materials
from cxdb.panels.panel import Panel
from cxdb.utils import FormPart, Input, Range, Select, table, RangeX, RangeS

# Mapping from project name to ProjectDescription class:
projects = {}


def project(name: str):
    """Decorator for filling in the projects dict."""
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
    panels: list[Panel] = []

    def postprocess(self, material: Material) -> None:
        pass

    def create_column_one(self, material: Material, materials: Materials):
        return '', ''

    def create_column_two(self, material: Material, materials: Materials):
        return '', ''


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


class ABS3BandStructurePanel(Panel):
    title = 'Electronic band-structure'

    def get_html(self,
                 material: Material,
                 materials: Materials) -> tuple[str, str]:
        uid = material.uid
        path = material.folder / f'abs3/{uid}.png'
        path.parent.mkdir(exist_ok=True)
        if not path.is_file():
            dbpath = material.folder / 'abs3.db'
            dct = connect(dbpath).get(id=uid).data
            ok = abs3_bs(dct, path)
            if not ok:
                return ('', '')
        return (
            f'<img alt="BS for {uid}" src="/abs3/png/{uid}" />', '')


def abs3_bs(d: dict, path: Path) -> bool:
    if 'X' not in d:
        return False
    n1 = len(d['X'])
    n2 = len(d['names'])
    if n1 != n2:
        print('bad data:', d, path)
        return False
    fig, ax = plt.subplots()
    ax.plot(d['x'], d['y'])
    ax.set_xticks(d['X'], d['names'])
    ax.set_xlim(0, d['X'][-1])
    ax.set_ylim(-6, 5)
    ax.set_ylabel('GLLB-SC energy relative to VBM [eV]')
    plt.savefig(path)
    plt.close()
    return True


@project('abs3')
class ABS3ProjectDescription(ProjectDescription):
    title = 'Database of ABS3 materials'
    column_names = {
        'E_hull': 'E-hull [eV]',
        'm_e': 'Effective electron mass [m<sub>e</sub>]',
        'm_h': 'Effective hole mass [m<sub>e</sub>]',
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
    panels = [ABS3BandStructurePanel()]


@project('abx2')
class ABX2ProjectDescription(ProjectDescription):
    title = 'Database of ABX2 materials'
    column_names = {
        'E_hull': 'E-hull [eV]',
        'KS_gap': 'Kohn Sham band gap [eV]',
        'm_e': 'Effective electron mass [m<sub>e</sub>]',
        'm_h': 'Effective hole mass [m<sub>e</sub>]',
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
        'g0w0_gap': 'G0W0 gap at Gamma [eV]',
        'gw0_gap': 'GW0 gap at Gamma [eV]',
        'gw_gap': 'GW gap at Gamma [eV]',
        'hse06_gap': 'HSE06 gap at Gamma [eV]',
        'lda_gap': 'LDA gap at Gamma [eV]',
        'gllbsc_gap': 'GLLBSC gap at Gamma [eV]'}
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
        'm_e': 'Effective electron mass [m<sub>e</sub>]',
        'm_h': 'Effective hole mass [m<sub>e</sub>]',
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

    def create_column_one(self,
                          material: Material,
                          materials: Materials) -> tuple[str, str]:
        return '\n'.join(
            table([header, ''],
                  materials.table(material, names))
            for header, names in [
                ('Host properties',
                 ['host',
                  'supercell',
                  'host_spacegroup',
                  'hostenergy']),
                ('Structural properties',
                 ['dopant',
                  'defecttype',
                  'depth',
                  'extension_factor',
                  'site']),
                ('Electronic properties',
                 ['spin',
                  'en2',
                  'eform',
                  'dopant_chemical_potential',
                  'conv2'])]), ''


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


@project('bidb')
class BiDBProjectDescription(ProjectDescription):
    title = 'Bilayer database'
    column_names = {
        'binding_energy_zscan': 'Binding energy (zscan) [meV/Å<sup>2</sup>]',
        'number_of_layers': 'nlayers',
        'monolayer_uid': 'Monolayer ID',
        'bilayer_uid': 'Bilayer ID',
        'dynamically_stable': 'Dynamically stable',
        'magnetic': 'Magnetic',
        'interlayer_magnetic_exchange': 'Interlayer Magnetic State',
        'slide_stability': 'Slide Stability',
        'thermodynamic_stability': '',
        'binding_energy_gs': 'Binding Energy (gs) [meV/Å<sup>2</sup>]',
        'ehull': 'Energy above convex hull [eV/atom]',
        'gap_pbe': 'Band gap (PBE)',
        'icsd_id': 'ICSD id of parent bulk structure',
        'cod_id': 'COD id of parent bulk structure',
        'layer_group': 'Layer group',
        'layer_group_number': 'Layer group number',
        'space_group': 'Space group',
        'space_group_number': 'Space group number'}
    uid = 'uid'
    initial_columns = [
        'formula',
        'number_of_layers',
        'binding_energy_gs',
        'slide_stability',
        'uid',
        'magnetic']
    form_parts = [
        Select('Number of layers', 'number_of_layers', ['', '1', '2']),
        Range('Binding energy [meV/Å<sup>2</sup>] (bilayers)',
              'binding_energy_gs'),
        Select('Slide stability (bilayers)', 'slide_stability',
               ['', 'Stable']),
        Range('Band gap range [eV]', 'gap_pbe'),
        Select('Magnetic', 'magnetic', ['', '0', '1'])]

    def postprocess(self, material: Material):
        icsd_id = getattr(material, 'icsd_id', 0)
        if icsd_id:
            material.add_column('icsd_id', str(icsd_id), update=True)
        cod_id = getattr(material, 'cod_id', 0)
        if cod_id:
            material.add_column('cod_id', str(cod_id), update=True)

    def create_column_one(self,
                          material: Material,
                          materials: Materials) -> tuple[str, str]:
        def tab(names):
            return materials.table(material, names)

        if material.number_of_layers == 1:
            tables = [
                (['Monolayer structure info', 'Value'],
                 tab(['layer_group',
                      'layer_group_number',
                      'space_group',
                      'space_group_number',
                      'icsd_id',
                      'cod_id'])),
                (['Stability', ''],
                 tab(['dynamically_stable', 'ehull'])),
                (['Basic properties', ''],
                 tab(['magnetic', 'gap_PBE']))]
        else:
            tables = [
                (['Bilayer structure info', 'Value'],
                 tab(['layer_group',
                      'layer_group_number',
                      'space_group',
                      'space_group_number',
                      'icsd_id',
                      'cod_id',
                      'monolayer_uid'])),
                (['Stability', ''],
                 tab(['binding_energy_zscan', 'slide_stability'])),
                (['Basic properties', ''],
                 tab(['magnetic', 'interlayer_magnetic_exchange', 'gap_PBE']))]
            # Add link to monolayer:
            _, mid = tables[0][1][-1]
            tables[0][1][-1] = (
                'Monolayer in BiDB',
                f'<a href={mid}>{mid}</a>')
            tables[0][1].append(
                ('Monolayer in C2DB',
                 f'<a href=https://cmrdb.fysik.dtu.dk/c2db/row/{mid}>'
                 f'{mid}</a>'))
        return '\n'.join(table(header, rows) for header, rows in tables), ''


@project('lowdim')
class LowDimProjectDescription(ProjectDescription):
    title = ('Definition of a scoring parameter to'
             ' identify low-dimensional materials components')
    uid = 'dbid'
    column_names = {
        's_0': '0D score',
        's_1': '1D score',
        's_2': '2D score',
        's_3': '3D score',
        's_01': '0D+1D score',
        's_02': '0D+2D score',
        's_03': '0D+3D score',
        's_12': '1D+2D score',
        's_13': '1D+3D score',
        's_23': '2D+3D score',
        's_012': '0D+1D+2D score',
        's_013': '0D+1D+3D score',
        's_023': '0D+2D+3D score',
        's_123': '1D+2D+3D score',
        's_0123': '0D+1D+2D+3D score',
        'a_0': 'Start of 0D k-interval',
        'a_1': 'Start of 1D k-interval',
        'a_2': 'Start of 2D k-interval',
        'a_3': 'Start of 3D k-interval',
        'a_01': 'Start of 0D+1D k-interval',
        'a_02': 'Start of 0D+2D k-interval',
        'a_03': 'Start of 0D+3D k-interval',
        'a_12': 'Start of 1D+2D k-interval',
        'a_13': 'Start of 1D+3D k-interval',
        'a_23': 'Start of 2D+3D k-interval',
        'a_012': 'Start of 0D+1D+2D k-interval',
        'a_013': 'Start of 0D+1D+3D k-interval',
        'a_023': 'Start of 0D+2D+3D k-interval',
        'a_123': 'Start of 1D+2D+3D k-interval',
        'a_0123': 'Start of 0D+1D+2D+3D k-interval',
        'b_0': 'End of 0D k-interval',
        'b_1': 'End of 1D k-interval',
        'b_2': 'End of 2D k-interval',
        'b_3': 'End of 3D k-interval',
        'b_01': 'End of 0D+1D k-interval',
        'b_02': 'End of 0D+2D k-interval',
        'b_03': 'End of 0D+3D k-interval',
        'b_12': 'End of 1D+2D k-interval',
        'b_13': 'End of 1D+3D k-interval',
        'b_23': 'End of 2D+3D k-interval',
        'b_012': 'End of 0D+1D+2D k-interval',
        'b_013': 'End of 0D+1D+3D k-interval',
        'b_023': 'End of 0D+2D+3D k-interval',
        'b_123': 'End of 1D+2D+3D k-interval',
        'b_0123': 'End of 0D+1D+2D+3D k-interval',
        'numc_0': 'Number of 0D components',
        'numc_1': 'Number of 1D components',
        'numc_2': 'Number of 2D components',
        'numc_3': 'Number of 3D components',
        'source': 'Source',
        'dbid': 'ID #',
        'doi': 'DOI',
        'publication': 'Publication',
        'spacegroup_number': 'Space group #',
        'dimtype': 'Dimensionality',
        'h': 'Component count',
        'warning': 'Warning'}
    initial_columns = ['formula', 'source', 'dbid',
                       's_0', 's_1', 's_2', 's_3']
    form_parts = [
        RangeX(
            'Score range', 's',
            ['0D', '1D', '2D', '3D',
             '0D+1D', '0D+2D', '0D+3D', '1D+2D', '1D+3D', '2D+3D',
             '0D+1D+2D', '0D+1D+3D', '1D+2D+3D',
             '0D+1D+2D+3D'],
            ['s_0', 's_1', 's_2', 's_3',
             's_01', 's_02', 's_03', 's_12', 's_13', 's_23',
             's_012', 's_013', 's_123',
             's_0123']),
        Select('Database source', 'source', ['', 'COD', 'ICSD'])]
    panels = [LowDimPanel()]

    def postprocess(self, material: Material):
        material.add_column('dbid', str(material.dbid), update=True)

    def create_column_one(self,
                          material: Material,
                          materials: Materials) -> tuple[str, str]:
        rows = materials.table(material, keysfortable0)
        doi = material.get('doi')
        if doi:
            href = f'<a href="https://doi.org/{doi}">{doi}</a>'
            rows.append(('doi', href))
        if material.source == 'COD':
            id = material.dbid
            href = ('<a href="http://www.crystallography.net/cod/' +
                    f'{id}.html">{id}</a>')
            rows.insert(0, ('COD Number', href))
        else:
            assert material.source == 'ICSD'
            rows.insert(0, ('ICSD Number', material.dbid))
        return table(['Basic properties', ''], rows), ''

    def create_column_two(self,
                          material: Material,
                          materials: Materials) -> tuple[str, str]:
        if material.source == 'ICSD':
            return ('not allowed to show atoms', '')
        return ('', '')  # use default AtomsPanel behavior


COD = 'https://www.crystallography.net/cod/'
ICSD = 'https://icsd.products.fiz-karlsruhe.de/en/'


@project('c1db')
class C1DBProjectDescription(ProjectDescription):
    title = 'Computational 1D materials database'
    column_names = {
        'PBE_1D': 'uid for 1D material calculated using PBE',
        'PBED3_1D': 'uid for 1D material calculated using PBE-D3',
        'PBED3_3D': 'uid for 3D material calculated using PBE-D3',
        'Source': 'Source',
        'derived_from': 'derived from',
        'xc': 'XC-functional',
        'ndim': 'Dimensionality'}
    initial_columns = [
        'formula', 'hform', 'gap', 'is_magnetic', 'xc', 'ndim']
    uid = 'uid'
    form_parts = [
        Select('Dimensionality', 'pbc', ['', 'FFT', 'TTT'], ['', '1D', '3D']),
        Select('XC-functional', 'calculator', ['', 'dftd3', 'gpaw'],
               ['', 'PBE+D3', 'PBE']),
        Select('Source', 'source',
               ['',
                'COD',
                'ICSD',
                'Derived by element substitution',
                'Machine learning generated']),
        Select('Dynamically stable (phonons)', 'dyn_phonons',
               ['', 'low', 'high'], ['', 'No', 'Yes']),
        RangeS('Thermodynamic stability', 'thermodynamic_stability_level',
               ['1', '2', '3'], ['Low', 'Medium', 'High']),
        Select('Magnetic', 'is_magnetic',
               ['', 'True', 'False'], ['All', 'Yes', 'No']),
        RangeX('Band gap range [eV]', 'xc',
               ['gap', 'gap_hse'], ['PBE', 'HSE06@PBE'])]

    def create_column_one(self, material, materials):
        rows = materials.table(material, self.column_names)
        source = material.Source
        df = material.get('derived_from')
        if source == 'COD':
            rows.append(
                ['Source',
                 f'<a href={COD}/{df}.html>COD {df}</a>'])
        elif source == 'ICSD':
            rows.append(
                ['Source',
                 f'<a href={ICSD}>ICSD {df}</a>'])
        elif source == 'Machine learning generated':
            rows.append(
                ['Source', source])
        else:
            rows += [
                ['Source', source],
                ['Derived form',
                 f'<a href={df}>{df}</a>']]

        for key, text in [('PBED3_1D', '1D (PBE-D3)'),
                          ('PBE_1D', '1D (PBE)'),
                          ('PBED3_3D', '3D (PBE-D3)')]:
            uid = getattr(material, key, '')
            if uid:
                rows.append(
                    [text, f'<a href={uid}>{uid}</a>'])
        return table(['XXX', ''], rows), ''


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
