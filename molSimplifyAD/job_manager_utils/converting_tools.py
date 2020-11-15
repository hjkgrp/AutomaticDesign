import os, sys
import glob
import numpy as np
from molSimplifyAD.ga_tools import get_mulliken, rename_ligands, find_files_by_name
from molSimplifyAD.process_scf import read_molden_file
from molSimplify.Classes.ligand import get_lig_symmetry


def isCSD(job):
    iscsd = True
    for ii in range(6):
        if not job[ii].isupper():
            iscsd = False
            break
    return iscsd

def isMutation(job):
    return ('mutation' in job)

def isSP(outfile):
    issp = False
    with open(outfile, 'r') as fo:
        txt = "".join(fo.readlines()[:400])
    if "SINGLE POINT ENERGY CALCULATIONS" in txt:
        issp = True
    return issp


def bind_complex_info(this_run, name, spin):
    nn = name.split('_')
    try:
        this_run.metal, this_run.ox, this_run.spin = nn[0], int(nn[1]), int(spin)
    except:
        print("Incompatible namig scheme for db complexes (non-CSD). Should be: <metal>_<ox>_<spin>_<lig1>...")
    this_run.liglist = nn[3:]
    this_run.liglist_compact = rename_ligands(this_run.liglist)
    this_run.refcode = False
    this_run.ligcharge = int(this_run.charge) - this_run.ox


def bind_csd_info(this_run, name, spin, init_mol):
    nn = name.split('_')
    this_run.refcode, possible_suffix, this_run.spin = nn[0], nn[1], int(spin)
    if possible_suffix.isdigit():
        this_run.refcode += '-%s' % possible_suffix
    this_run.liglist, this_run.liglist_compact = [this_run.refcode], [this_run.refcode]
    this_run.ligcharge = False
    this_run.metal = init_mol.getAtom(init_mol.findMetal()[0]).symbol().lower()


def collect_spin_info(this_run, spin, ss_act, ss_target):
    if spin == 1:
        this_run.ss_target = 0
        this_run.ss_act = 0
    else:
        try:
            this_run.ss_target = float(ss_target.strip(')'))
            this_run.ss_act = float(ss_act)
        except:
            this_run.ss_target = np.nan
            this_run.ss_act = np.nan
    return this_run


def check_conv(this_run, tot_time, energy, output):
    this_run.converged = False
    new_opt = output.wordgrab(['Optimization Converged'], 'whole_line')[0][0]
    old_opt = output.wordgrab(['Converged!'], 'whole_line')[0][0]
    if (old_opt != None) or (new_opt != None):
        if (tot_time != None) and (energy != None):
            this_run.converged = True


def calculate_mulliken_spins(this_run):
    multiwfnpath = find_files_by_name(this_run.scrpath_real, ".molden")
    if len(multiwfnpath) > 0:
        multiwfnpath = multiwfnpath[0]
        if (this_run.iscsd or this_run.isMutation):
            mulliken_spin_list = get_mulliken(
                multiwfnpath, this_run.spin, external=True)
        else:
            mulliken_spin_list = get_mulliken(
                multiwfnpath, this_run.spin, this_run.liglist[-1], external=True)
        print(("mulliken spins: ", mulliken_spin_list))
        this_run.net_metal_spin = mulliken_spin_list[0]
        if len(mulliken_spin_list) > 1:
            this_run.net_oxygen_spin = mulliken_spin_list[1]
    else:
        print("No molden path found.")


def get_initgeo(filename):
    with open(filename, "r") as fo:
        numatoms = int(fo.readline().split()[0])
    with open(filename, "r") as fo:
        xyz = fo.readlines()[:numatoms + 2]
    return xyz


def get_lastgeo(filename):
    with open(filename, "r") as fo:
        numatoms = int(fo.readline().split()[0])
    with open(filename, "r") as fo:
        xyz = fo.readlines()[-(numatoms + 2):]
    return xyz


def obtain_wavefunction_molden(this_run):
    scrdir = this_run.scrpath_real
    wavefunc_keys = ["c0", "ca0", "cb0"]
    this_run.wavefunction = {}
    this_run.wavefunction_path = {}
    for key in wavefunc_keys:
        wavefunc_file = scrdir + key
        if os.path.isfile(wavefunc_file):
            # print("found %s. Writting into DFTrun..." % key)
            with open(wavefunc_file, "rb") as fo:
                wf = fo.read()
        else:
            wf = False
            wavefunc_file = False
        this_run.wavefunction.update({key: wf})
        this_run.wavefunction_path.update({key: wavefunc_file})
    this_run.molden_path = find_files_by_name(this_run.scrpath_real, ".molden")
    if this_run.molden_path:
        this_run.molden_path = this_run.molden_path[0]
        with open(this_run.molden_path, "r") as fo:
            this_run.molden = fo.read()
    else:
        this_run.molden_path, this_run.molden = False, False


def grab_dipole_moment(outfile):
    diople_vec, diople_moment = np.nan, np.nan
    with open(outfile, 'r') as fo:
        for line in fo:
            if "DIPOLE MOMENT:" in line:
                diople_vec = [float(x) for x in line.split('{')[-1].split('}')[0].split(",")]
                diople_moment = float(line.split('{')[-1].split('}')[1].split()[-2].strip(')'))
    return diople_vec, diople_moment


def get_ligsymmetry_graphdet(optmol):
    try:
        ligsymmetry = get_lig_symmetry(optmol)
    except:
        ligsymmetry = "undef"
    try:
        det = optmol.get_mol_graph_det(oct=True)
    except:
        det = "undef"
    return ligsymmetry, det



def get_dynamic_feature(this_run):
    pass
