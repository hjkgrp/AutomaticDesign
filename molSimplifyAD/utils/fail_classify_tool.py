from molSimplify.Classes import globalvars
from molSimplify.Classes.mol3D import *
import sys, os,shutil,datetime
import time
from molSimplify.Informatics.autocorrelation import*
from molSimplify.Informatics.misc_descriptors import*
from molSimplify.Informatics.RACassemble import*
from molSimplify.Informatics.graph_analyze import*
from molSimplifyAD.ga_tools import* 
from molSimplifyAD.post_classes import* 

def detect_hydrogen_transfer(mol,init_mol):
    hv_f = [] 
    for i,ats in enumerate(mol.getAtoms()):
        if not ats.symbol() == "H":
            
            num_h_fin = mol.getBondedAtomsH(i)
            num_h_init = init_mol.getBondedAtomsH(i)
            hv_f.append(abs(len(num_h_init)-len(num_h_fin)))
    if sum(hv_f) != 0:
        return True ## HAT
    else: 
        return False ## no HAT
def maximum_min_pairwise_dist(mol):
    max_min_dist = 0
    for atoms in mol.getAtoms():
        this_atom_min_dist = 10000
        for other_atoms in mol.getAtoms():
            this_dist = distance(atoms.coords(),other_atoms.coords())
            if this_dist > 0:
                if this_dist < this_atom_min_dist:
                    this_atom_min_dist = this_dist
        if this_atom_min_dist > max_min_dist:
            max_min_dist = this_atom_min_dist
    return(max_min_dist)

def find_problems(initial_path,final_path):
    # types  
    #        geo-miss: no final geo found
    #        lig-break-fail: couldn't unpack ligands 
    #        fall-apart : something falls way off 
    #        metal-undercoord : metal is not oct 
    #        dimerize  : form bridge
    #        HAT  : hydrogen transfer
    #        OTHER: who knows?

    this_type= ""
    init = mol3D()
    fin = mol3D()
    valid = True # flag to continue
   
    ## check if geos exist:
    if os.path.isfile(initial_path):
        try:
            init.readfromxyz(initial_path)
        except:
             valid = False
    else:
        valid = False
        
    if os.path.isfile(final_path):
        try:
            fin.readfromxyz(final_path)
        except:
             valid = False
    else:
        valid = False
        
    if not valid or not init.natoms == fin.natoms:
        this_type= "geo-miss"
        return(this_type)
    else:
        # get max distance 
        max_pd = maximum_min_pairwise_dist(fin)
        if max_pd >= 3.00:
            
            this_type= "fall-apart"
            return(this_type)
        else:
            ## check coordination
            connect_to_metal = fin.getBondedAtomsSmart(fin.findMetal()[0],oct=True)
            if len(connect_to_metal) < 6:
                this_type = "metal-undercoord"            
                return(this_type)
            ## check for unwanted HAT
            elif detect_hydrogen_transfer(mol=fin,init_mol=init):
                this_type = "HAT"
                return(this_type)
    
            ## try to get descriptors:
            try:

                init_names, init_desc  = get_descriptor_vector(init)
                fin_names, fin_desc  = get_descriptor_vector(fin)
                ## check dent
                max_fin_dent  =  max(fin_desc[fin_names.index("misc-dent-eq")],fin_desc[fin_names.index("misc-dent-ax")])
                max_init_dent  =  max(init_desc[init_names.index("misc-dent-eq")],init_desc[init_names.index("misc-dent-ax")])


                if max_fin_dent >  max_init_dent:
                    this_type = "dimmerize"
                else:
                    this_type = "OTHER"
                return(this_type)
            
            except:
                this_type = "lig-break-fail"
                return(this_type)


paths = setup_paths()
print('\n checking convergence of jobs\n')
## set up environment:        
path_dictionary = setup_paths()
base_path_dictionary = setup_paths()
## find converged jobs
converged_jobs = find_converged_job_dictionary()

outcomes_dictionary = dict()

print('found ' + str(len(converged_jobs.keys())) + ' converged jobs to check')


counter = 0
for jobs in converged_jobs.keys():
    job_status = converged_jobs[jobs]
    gene, gen, slot, metal, ox, eqlig, axlig1, axlig2, eqlig_ind, axlig1_ind, axlig2_ind, spin, spin_cat, ahf, base_name, base_gene = translate_job_name(jobs)
    if ahf == 20:
        counter +=1
        path_dictionary = setup_paths()
        path_dictionary = advance_paths(path_dictionary, gen)  ## this adds the /gen_x/ to the paths
        if job_status in ["1","8"]:
                print('bad job ' + str(jobs) + ' with status ' + str(job_status))
                this_run = DFTRun(base_name)
                this_run.geopath = (path_dictionary["optimial_geo_path"] + base_name + ".xyz")
                this_run.progpath = (path_dictionary["prog_geo_path"] + base_name + ".xyz")
                this_run.init_geopath = (path_dictionary["initial_geo_path"] + base_name + ".xyz")
                if job_status == '8':
                    check_geo_path = this_run.progpath
                    reason = find_problems(this_run.init_geopath,this_run.progpath)
                    print(reason,this_run.progpath)
                else:
                    check_geo_path = this_run.geopath
                    reason = find_problems(this_run.init_geopath,this_run.geopath)
                    print(reason,this_run.geopath)
                print('*******')

        elif job_status in ["3"]:
            check_geo_path = path_dictionary["initial_geo_path"] + base_name + ".xyz"
            reason = "ERROR-STAT"
        elif job_status in ["2"]:
            check_geo_path = path_dictionary["prog_geo_path"] + base_name + ".xyz"
            reason = "TO-BE-RESTART"
        elif job_status in ["7"]:
            check_geo_path = path_dictionary["prog_geo_path"] + base_name + ".xyz"
            reason = "MAX-RESUB"
        elif job_status in ["0"]:
            check_geo_path = path_dictionary["optimial_geo_path"] + base_name + ".xyz"
            reason = "SUCCESS"
        else:
            print('unkown status ' + str(job_status))
            check_geo_path = path_dictionary["initial_geo_path"] + base_name + ".xyz"
            reason = "UNKNOWN-STATUS-"+ str(job_status)
        if reason in outcomes_dictionary.keys():
            this_list  = outcomes_dictionary[reason]
            this_list.append(check_geo_path)
            outcomes_dictionary.update({reason:this_list})
        else:
            outcomes_dictionary.update({reason:[check_geo_path]})        
print('*******')
print('Found a total of ' + str(counter) + ' HFX-20 jobs to analyze') 
print('*******')
for keys in outcomes_dictionary.keys():
    print('There are '+ str(len(outcomes_dictionary[keys])) +  ' with status  ' + keys )
    with open('runs-'+keys+'.txt','w') as f:
        for element in outcomes_dictionary[keys]:
            f.write(element + '\n')




            

