'''

This is a script to push pre-mAD data to the db.

**IMPORTANT NOTE**

This script is just a scaffold. It has to be modified accordingly to fit one's data structure in order to make it work.

'''
import os
from os.path import expanduser
import glob
import time
import pandas as pd
import numpy as np
from molSimplifyAD.post_classes import DFTRun
from molSimplifyAD.dbclass_mongo import tmcMongo
from molSimplify.Classes.mol3D import mol3D
from molSimplifyAD.utils.pymongo_tools import connect2db, insert
from molSimplifyAD.ga_tools import get_mulliken
from molSimplifyAD.process_scf import read_molden_file


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


def get_final_energy(outfile):
    energy = np.nan
    if os.path.isfile(outfile):
        with open(outfile, 'r') as fo:
            for line in fo:
                if "FINAL ENERGY" in line:
                    energy = float(line.split()[-2])
    return energy


def get_wavefunction(scrpath):
    wfn = {"c0": False, "ca0": False, "cb0": False}
    if os.path.isfile(scrpath + "ca0"):
        with open(scrpath + "ca0", "rb") as fo:
            wfn.update({"ca0": fo.read()})
    if os.path.isfile(scrpath + "cb0"):
        with open(scrpath + "cb0", "rb") as fo:
            wfn.update({"cb0": fo.read()})
    return wfn


def rename_ligands(liglist):
    lig_dict = {"c3c(P(c1ccccc1)c2ccccc2)cccc3": "pph3",
                "c1ccncc1": "pyr",
                }
    for idx, lig in enumerate(liglist):
        if lig in list(lig_dict.keys()):
            liglist[idx] = lig_dict[lig]
    return liglist


update_fields = ['alphaHOMO', 'betaHOMO', 'alphaLUMO', 'betaLUMO', 'status']
print(("update_fields: ", update_fields))
df = pd.read_csv("spectro_all_RACs_format.csv")
home = expanduser("~")
dbconfig = json.load(open(home + "/.db_config"))
db = connect2db(user=dbconfig['user'], pwd=dbconfig['pwd'],
                host='localhost', port=27017,
                database='tmc', auth=True)
print(("# of complexes: ", db.oct.count()))
basepath = "/home/crduan/Binary_Classifier/molecule/"
count, merged = 0, 0
for idx, row in df.iterrows():
    if row["mad"] == "oldgen":
        print("=================")
        print((row["job_name"]))
        filepath = basepath + row["metal"] + '/' + "_".join(row["job_name"].split("_")[:-2]) + '/geometry/'
        outpath = filepath + "_".join(row["job_name"].split("_")[:3]) + "_" + "_".join(
            row["job_name"].split("_")[-2:]) + "_geometry__runlast.out"
        scrpath = filepath + "scr_" + "_".join(row["job_name"].split("_")[-2:]) + "_geometry__runlast/"
        print(filepath)
        this_run = DFTRun(row["job_name"])
        this_run.name = row["job_name"]
        this_run.metal = row["metal"]
        this_run.ox = int(row["oxstate"])
        this_run.spin = int(row["spinmult"])
        this_run.lig1 = row["eqlig"]
        this_run.lig2 = row["eqlig"]
        this_run.lig3 = row["eqlig"]
        this_run.lig4 = row["eqlig"]
        this_run.lig5 = row["axlig1"]
        this_run.lig6 = row["axlig2"]
        this_run.liglist = [this_run.lig1, this_run.lig2, this_run.lig3,
                            this_run.lig4, this_run.lig5, this_run.lig6]
        this_run.liglist_compact = rename_ligands(this_run.liglist)
        this_run.alpha = 20.0
        this_run.basis = 'lacvps_ecp'
        this_run.terachem_version = "v1.9-2018.04-dev"
        this_run.functional = "b3lyp"
        this_run.energy = get_final_energy(outpath)
        this_run.converged = True if int(row["flag_conv"]) else False

        this_run.scrpath = scrpath
        current_folder = this_run.scrpath.strip("optim.xyz")
        multiwfnpath = glob.glob(current_folder + "*.molden")
        if len(multiwfnpath) > 0:
            multiwfnpath = multiwfnpath[0]
            mulliken_spin_list = get_mulliken(multiwfnpath, this_run.spin, this_run.liglist[-1])
            print(mulliken_spin_list)
            this_run.net_metal_spin = mulliken_spin_list[0]
            if len(mulliken_spin_list) > 1:
                this_run.net_oxygen_spin = mulliken_spin_list[1]
        else:
            print("No molden path found.")
        this_run.get_check_flags()
        if this_run.converged:
            read_molden_file(this_run)

        this_run.ss_target = row["ss_target"]
        this_run.ss_act = row["ss_actual"]
        this_run.charge = int(row['charge_lig']) + int(row['oxstate'])
        this_run.ligcharge = int(row['charge_lig'])
        if int(row["geo_use"]):
            this_run.geo_flag = int(row['flag_oct'])
        else:
            this_run.geo_flag = np.nan
        if this_run.converged:
            if this_run.geo_flag:
                this_run.status = 0
            else:
                this_run.status = 1
        else:
            this_run.status = 7
        if int(row["ss_use"]):
            this_run.ss_flag = int(row['ss_label'])
        else:
            this_run.ss_flag = np.nan
        this_run.metal_spin_flag = np.nan
        if not row["spinmult"] == 1 and this_run.converged:
            this_run.wavefunction = get_wavefunction(scrpath)
        else:
            this_run.wavefunction = False
        if this_run.converged:
            this_run.flag_oct = True if int(row["flag_oct"]) else False
        else:
            this_run.flag_oct = -1
        filename = "./oldgen_" + this_run.name + "_optim.xyz"
        this_run.init_mol = mol3D()
        this_run.init_mol.readfromtxt(get_initgeo(filename))
        this_run.mol = mol3D()
        this_run.mol.readfromtxt(get_lastgeo(filename))
        this_run.progmol = mol3D()
        this_run.progmol.readfromtxt(get_lastgeo(filename))

        this_tmc = tmcMongo(this_run=this_run, tag="spectro", subtag="oldgen",
                            publication="Duan_JCTC_2019", update_fields=update_fields)
        print((this_tmc.document['status'], this_tmc.document['alphaHOMO']))
        _s = time.time()
        insetred = insert(db, "oct", this_tmc)
        print(("elapse: ", time.time() - _s))
        if insetred:
            count += 1
        else:
            merged += 1
print(("add %d entries in the %s['%s']." % (count, "tmc", "oct")))
print(("merge %d entries in the %s['%s']." % (merged, "tmc", "oct")))
