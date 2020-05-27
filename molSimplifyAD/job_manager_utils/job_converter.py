import os, sys
import glob
import numpy as np
from molSimplifyAD.post_classes import DFTRun
from molSimplify.Classes.mol3D import mol3D
from molSimplifyAD.ga_tools import get_mulliken, rename_ligands, translate_job_name
from molSimplifyAD.process_scf import read_molden_file
from molSimplify.job_manager.tools import list_active_jobs, extract_optimized_geo
from molSimplify.job_manager.classes import resub_history, textfile
from molSimplifyAD.job_manager_utils.bind_functions import bind_solvent, bind_water, bind_thermo, bind_vertIP, \
    bind_ligdissociate, bind_vertEA, bind_functionals
from molSimplifyAD.job_manager_utils.converting_tools import *

associated_jobs = {'solvent': bind_solvent,
                   'watercont': bind_water,
                   'thermo': bind_thermo,
                   'vertIP': bind_vertIP,
                   'vertEA': bind_vertEA,
                   'functionals': bind_functionals,
                   'dissociation': bind_ligdissociate}


def collect_base_jobs(path=False):
    if not path:
        path = os.getcwd()
    basejobs = list()
    for dirpath, dirs, files in os.walk(path):
        for file in sorted(files):
            if (file.split('.')[-1] == 'out') and (not any(c in file for c in list(associated_jobs.keys()))) and (
                    not any(c in dirpath for c in list(associated_jobs.keys()))) and not 'nohup' in file:
                basejobs.append([dirpath, '.'.join(file.split('.')[:-1])])
    return basejobs


def common_processing(jobname, basedir, output, outfile, spin, 
                      dbname=False, gene=False):
    dbname = jobname if dbname == False else dbname
    this_run = DFTRun(dbname, external=True)
    this_run.name = dbname
    this_run.octahedral = True
    if gene != False:
        # Assigns a gene, which is used during design
        temp_name = "_".join(jobname.split('_')[0:13]) #stops before job manager extra keywords
        translate_dict = translate_job_name(temp_name)
        this_run.gene = translate_dict['gene']
    energy, ss_act, ss_target, tot_time, thermo_grad_error, solvent_cont, tot_step, d3_energy = output.wordgrab(
        ['FINAL', 'S-SQUARED:', 'S-SQUARED:', 'processing', 'Maximum component of gradient is too large',
         'C-PCM contribution to final energy:', 'Optimization Cycle', 'DFTD Dispersion Correction:'],
        [2, 2, 4, 3, 0, 4, 3, -2], last_line=True)
    iscsd = isCSD(dbname)
    this_run.iscsd = iscsd
    this_run.charge = output.wordgrab(['Total charge'], -1)[0][0]
    try:
        this_run.charge = int(this_run.charge)
    except:
        this_run.charge = False
    if not d3_energy == None:
        this_run.d3opt_flag = True
    else:
        this_run.d3_energy = np.nan
        this_run.d3opt_flag = False
    if not this_run.iscsd:
        print(("jobname: ", dbname))
        bind_complex_info(this_run, dbname, spin)
    else:
        init_mol = mol3D()
        init_mol.readfromxyz(basedir + '/%s.xyz' % jobname)
        bind_csd_info(this_run, dbname, spin, init_mol)
    this_run.outpath = outfile
    this_run.tot_time = tot_time
    this_run.alpha = output.wordgrab(['Hartree-Fock exact exchange:'], -1)[0][0] * 100
    this_run.functional = output.wordgrab(['DFT Functional requested:'], -1)[0][0]
    this_run.terachem_version = output.wordgrab(['TeraChem'], 2)[0][0]
    this_run.alpha_level_shift = output.wordgrab(['Alpha level shift'], -1)[0][0]
    this_run.beta_level_shift = output.wordgrab(['Beta level shift'], -1)[0][0]
    this_run.basis = output.wordgrab(['Using basis set:'], -1)[0][0]
    if not energy == None:
        this_run.energy = float(energy)
    else:
        this_run.energy = np.nan
    this_run = collect_spin_info(this_run, spin, ss_act, ss_target)
    scrpath = basedir + '/scr/'
    this_run.scrpath_real = scrpath
    if os.path.isdir(this_run.scrpath_real):
        calculate_mulliken_spins(this_run)
    return this_run


def process_single_points(this_run, basedir, output):
    energy, ss_act, ss_target, tot_time, thermo_grad_error, solvent_cont, tot_step = output.wordgrab(
        ['FINAL', 'S-SQUARED:', 'S-SQUARED:', 'processing', 'Maximum component of gradient is too large',
         'C-PCM contribution to final energy:', 'Optimization Cycle'],
        [2, 2, 4, 3, 0, 4, 3], last_line=True)
    this_run.geo_opt = False
    this_run.converged = True if not energy == None else False
    if this_run.converged:
        this_run.tot_step = tot_step
        this_run.init_energy = energy
        scrpath = basedir + '/scr/'
        optimpath = scrpath + 'xyz.xyz'
        this_run.scrpath = optimpath
        this_run.init_mol = mol3D()
        this_run.init_mol.readfromtxt(get_initgeo(optimpath))
        this_run.mol = mol3D()
        this_run.mol.readfromtxt(get_initgeo(optimpath))
        obtain_wavefunction_molden(this_run)
        read_molden_file(this_run)
    return this_run


def process_geometry_optimizations(this_run, basedir, outfile, output):
    energy, ss_act, ss_target, tot_time, thermo_grad_error, solvent_cont, tot_step = output.wordgrab(
        ['FINAL', 'S-SQUARED:', 'S-SQUARED:', 'processing', 'Maximum component of gradient is too large',
         'C-PCM contribution to final energy:', 'Optimization Cycle'],
        [2, 2, 4, 3, 0, 4, 3], last_line=True)
    check_conv(this_run, tot_time, energy, output)
    this_run.geo_opt = True
    scrpath = basedir + '/scr/'
    optimpath = scrpath + 'optim.xyz'
    this_run.scrpath = optimpath
    if this_run.converged:
        if os.path.getsize(optimpath) > 45:
            this_run.init_energy = float(output.wordgrab(['FINAL ENERGY'], 'whole_line')[0][0][2])
            try:
                extract_optimized_geo(optimpath)
                this_run.geopath = scrpath + 'optimized.xyz'
            except:
                this_run.geopath = optimpath
            if os.path.isfile(outfile):
                diople_vec, diople_moment = grab_dipole_moment(outfile)
                this_run.diople_vec = diople_vec
                this_run.diople_moment = diople_moment
            else:
                raise ValueError("Output file does not exist: ", outfile)
            read_molden_file(this_run)
            obtain_wavefunction_molden(this_run)
            this_run.init_mol = mol3D()
            this_run.init_mol.readfromtxt(get_initgeo(optimpath))
            this_run.init_ligand_symmetry, this_run.init_mol_graph_det = get_ligsymmetry_graphdet(this_run.init_mol)
            this_run.mol = mol3D()
            this_run.mol.readfromtxt(get_lastgeo(optimpath))
            this_run.ligand_symmetry, this_run.mol_graph_det = get_ligsymmetry_graphdet(this_run.mol)
            try:
                this_run.check_oct_needs_init(debug=False, external=True)
            except:
                print("Warning: falied get geo_check!!!")
            this_run.obtain_rsmd()
            this_run.obtain_ML_dists()
            this_run.get_check_flags()
            this_run.get_optimization_time_step(current_path=outfile)
            if this_run.geo_flag:
                this_run.status = 0
            else:
                this_run.status = 1
        else:
            print(("Warning: optim file %s is empty. Skipping this run." % optimpath))
            this_run.converged = False
    else:
        this_run.status = 7
    return this_run


def jobmanager2mAD(job, active_jobs, dbname=False, gene=False):
    this_run = False
    basedir, jobname = job[0], job[1]
    outfile = basedir + '/' + jobname + '.out'
    if not (os.path.split(outfile.rsplit('_', 1)[0])[-1] in active_jobs) or ('nohup' in outfile):
        output = textfile(outfile)
        try:
            spin = int(output.wordgrab(['Spin multiplicity:'], -1)[0][0])
        except:
            print(('Cannot read file: ', outfile))
            return this_run
        this_run = common_processing(jobname, basedir, output, outfile, spin, dbname=dbname, gene=gene)
        issp = isSP(outfile)
        if not issp:
            this_run = process_geometry_optimizations(this_run, basedir, outfile, output)
            for a in list(associated_jobs.keys()):
                associated_jobs[a](this_run, jobname, basedir)
        else:
            this_run = process_single_points(this_run, basedir, output)
    return this_run


def loop_convert_jobs(path=False, gene=False):
    runs = dict()
    basejobs = collect_base_jobs(path=path)
    active_jobs = list_active_jobs()
    for job in basejobs:
        this_run = jobmanager2mAD(job, active_jobs, gene=gene)
        if this_run and this_run.converged:
            runs.update({'/'.join(job): this_run})
    return runs
