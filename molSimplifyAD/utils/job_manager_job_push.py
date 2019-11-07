'''
This script takes in jobs from the job manager and prepares them for the database.

** Important NOTE**
You must be careful with how the ligands are named in order to be consistent with the database structure. 
The job manager does not care how your jobs are named, and thus if you do not check this, you may push inconsistent data to the DB.

You MUST adjust the subtag and tag that will be used in the db.

'''

import os, sys
import glob
import time
import pandas as pd
import numpy as np
from bson.objectid import ObjectId
from molSimplifyAD.post_classes import DFTRun
from molSimplifyAD.dbclass_mongo import tmcMongo
from molSimplify.Classes.mol3D import mol3D
from molSimplifyAD.utils.pymongo_tools import connect2db, insert, query_db
from molSimplifyAD.ga_tools import get_mulliken, rename_ligands
from molSimplifyAD.process_scf import read_molden_file
from molSimplify.job_manager.tools import *

user = sys.argv[1]
pwd = sys.argv[2]
collection = False
tag = False
if (not collection) or (not tag):
    raise ValueError('Need to provide a collection and tag.')
write = False
outfiles = find('*out')

#### by default, the name of the complex will be assumed to be metal_ox_spin_lig1_lig2_lig3_lig4_lig5_lig6. By default HFX will be 20, unless resampling is done.
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

def get_wavefunction(scrpath):
    wfn = {"c0": False, "ca0": False, "cb0": False}
    if os.path.isfile(scrpath + "ca0"):
        with open(scrpath + "ca0", "rb") as fo:
            wfn.update({"ca0": fo.read()})
    if os.path.isfile(scrpath + "cb0"):
        with open(scrpath + "cb0", "rb") as fo:
            wfn.update({"cb0": fo.read()})
    return wfn



db = connect2db(user=user, pwd=pwd,
                host='localhost', port=27017,
                database='tmc', auth=True)
print("# of complexes: ", db.oct.count())
counter = 0

count = 0
merged = 0
for i, outfile in enumerate(outfiles):
    #if 'HFX' not in outfile:
    #    continue
    #else:
    #    counter += 1
    #    print('-----')
    #    print(outfile)
    #if counter > 10:
    #    sard
    if 'nohup' in outfile:
        continue
    output = textfile(outfile)
    spin = int(output.wordgrab(['Spin multiplicity:'],2)[0][0])
    finalenergy,s_squared,s_squared_ideal,time,thermo_grad_error,implicit_solvation_energy,geo_opt_cycles = output.wordgrab(['FINAL',
                                                                                        'S-SQUARED:',
                                                                                        'S-SQUARED:','processing',
                                                                                        'Maximum component of gradient is too large',
                                                                                        'C-PCM contribution to final energy:',
                                                                                        'Optimization Cycle'],
                                                                                        [2,2,4,3,0,4,3],last_line=True)
    name = output.wordgrab(['XYZ coordinates'],2)[0][0].strip('.xyz').split('HFXresampling')[0]
    metal = name.split('_')[0]
    ox = int(name.split('_')[1])
    lig1 = name.split('_')[3]
    lig2 = name.split('_')[4]
    lig3 = name.split('_')[5]
    lig4 = name.split('_')[6]
    lig5 = name.split('_')[7]
    lig6 = name.split('_')[8]
    this_run = DFTRun(name, external=True)
    this_run.name = name
    this_run.metal = metal
    this_run.ox = ox
    this_run.spin = spin
    this_run.lig1 = lig1
    this_run.lig2 = lig2
    this_run.lig3 = lig3
    this_run.lig4 = lig4
    this_run.lig5 = lig5
    this_run.lig6 = lig6
    this_run.liglist = [lig1,lig2,lig3,lig4,lig5,lig6]
    this_run.liglist_compact = rename_ligands(this_run.liglist)
    this_run.alpha = output.wordgrab(['Hartree-Fock exact exchange:'],-1)[0][0]*100
    this_run.functional = output.wordgrab(['DFT Functional requested:'],-1)[0][0]
    this_run.terachem_version = output.wordgrab(['TeraChem'],2)[0][0]
    this_run.alpha_level_shift = output.wordgrab(['Alpha level shift'],-1)[0][0]
    this_run.beta_level_shift = output.wordgrab(['Beta level shift'],-1)[0][0]
    if spin == 1:
        this_run.ss_target = 0
        this_run.ss_act = 0
    else:
        this_run.ss_target = s_squared_ideal
        this_run.ss_act = s_squared
    this_run.charge = int(output.wordgrab(['Total charge'],-1)[0][0])
    this_run.ligcharge = int(this_run.charge)-int(ox)
    new_opt = output.wordgrab(['Optimization Converged'],'whole_line')[0][0]
    old_opt = output.wordgrab(['Converged!'],'whole_line')[0][0]
    if (old_opt != None) or (new_opt != None):
        if (time != None) and (finalenergy != None):
            this_run.converged = True
            this_run.tot_time = time
            this_run.tot_step = geo_opt_cycles
        else:
            this_run.converged = False
    else:
        this_run.converged = False
    this_run.init_energy = float(output.wordgrab(['FINAL ENERGY'],'whole_line')[0][0][2]) #take the first set of this and take the 3rd value
    basepath = os.path.split(outfile)[0]
    scrpath = basepath+'/scr/'
    optimpath = scrpath+'optim.xyz'
    this_run.scrpath = optimpath
    multiwfnpath = glob.glob(scrpath + "*.molden")
    if len(multiwfnpath) > 0:
        multiwfnpath = multiwfnpath[0]
        mulliken_spin_list = get_mulliken(multiwfnpath, this_run.spin, this_run.liglist[-1],external=True)
        print(mulliken_spin_list)
        this_run.net_metal_spin = mulliken_spin_list[0]
        if len(mulliken_spin_list) > 1:
            this_run.net_oxygen_spin = mulliken_spin_list[1]
    else:
        print("No molden path found.")
    if this_run.converged:
        read_molden_file(this_run)
        extract_optimized_geo(optimpath)
    this_run.geopath = scrpath+'optimized.xyz'
    this_run.init_mol = mol3D()
    this_run.init_mol.readfromtxt(get_initgeo(optimpath))
    this_run.mol = mol3D()
    this_run.mol.readfromtxt(get_lastgeo(optimpath))
    flag_oct, flag_list, dict_oct_info = this_run.mol.IsOct()
    this_run.flag_oct = flag_oct
    if int(flag_oct) == 1:
        this_run.geo_flag = 1
    else:
        this_run.geo_flag = np.nan
    this_run.get_check_flags()
    this_run.flag_oct = flag_oct
    this_run.flag_list = flag_list
    this_run.dict_geo_check = dict_oct_info
    this_run.write_geo_dict()
    this_run.obtain_rsmd()
    this_run.obtain_ML_dists()
    if this_run.converged:
        try:
            wfn = get_wavefunction(scrpath)
        except:
            wfn = False
        this_run.wavefunction = wfn
    this_run.get_descriptor_vector()
    if this_run.converged:
        if this_run.geo_flag:
            this_run.status = 0
        else:
            this_run.status = 1
    else:
        this_run.status = 7
    while 'HFXresampling' in basepath:
        basepath = os.path.split(basepath)[0]
    dbinfo = glob.glob(basepath+'/*dbinfo')
    if len(dbinfo) > 0:
        with open(dbinfo[0]) as f:
            dbval = f.readlines()[0]
        reference_struct = query_db(db,'oct',{'_id':ObjectId(str(dbval))})
        print(reference_struct.count())
        print(reference_struct[0]['subtag'])
        subtag = reference_struct[0]['subtag']+'_additional'
    else:
        subtag = 'job_manager_jobs'
            
    print(this_run.geo_flag,this_run.metal_spin_flag,this_run.ss_flag,this_run.status,subtag)
    if write:
        this_tmc = tmcMongo(this_run=this_run, tag=tag, subtag=subtag)
        _s = time.time()
        inserted = insert(db, collection, this_tmc)
        print("elapse: ", time.time() - _s)
        if inserted:
            count += 1
        else:
            merged += 1

#### Do this after looping over all outfiles.
print("add %d entries in the %s['%s']." % (count, database, collection))
print("merge %d entries in the %s['%s']." % (merged, database, collection))
print("creating index...")
if write:
    db[collection].create_index([("metal", pymongo.ASCENDING),
                                 ("ox", pymongo.ASCENDING),
                                 ("spin", pymongo.ASCENDING),
                                 ("converged", pymongo.ASCENDING),
                                 ("alpha", pymongo.ASCENDING),
                                 ("lig1", pymongo.ASCENDING),
                                 ("lig5", pymongo.ASCENDING),
                                 ("lig6", pymongo.ASCENDING),
                                 ("status", pymongo.ASCENDING),
                                 ("geo_flag", pymongo.ASCENDING),
                                 ("ss_flag", pymongo.ASCENDING)
                                 ])
    dump_databse(user=user, pwd=pwd)
