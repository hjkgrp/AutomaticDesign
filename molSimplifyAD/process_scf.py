import glob
import string
import sys
import os
import numpy as np
import math
import random
import string
import numpy
from molSimplifyAD.ga_tools import *
from molSimplifyAD.post_classes import *
from molSimplifyAD.ga_oct_check import *
########### UNIT CONVERSION
HF_to_Kcal_mol = 627.509###
eV_to_Kcal_mol = 23.06055## 
###########################


def scfextract(a_run,list_of_props):
    ## function to create a list of run properties
    ##  from a dft run
    #  @param a_run the dft run 
    #  @param list_of_props list of keywords for properties
    #  @return list of properties
    extrct_props = []
    for keys in list_of_props:
        extrct_props.append(a_run.__dict__[str(keys)])
    
    return extrct_props


def test_terachem_sp_convergence(job):
    ## function to test single point convergence
    ##  for terachem files
    #  @param job a job name
    #  @return this_run populated run class
    ### get paths
    path_dictionary = setup_paths()
    gene,gen,slot,metal,ox,eqlig,axlig1,axlig2,eq_ind,ax1_ind,ax2_ind,spin,spin_cat,ahf,basename = translate_job_name(job)
    this_GA = get_current_GA()
    exchange = ahf
    alpha=float(exchange)
    ### flag
    converged =  False
    ## set up up
    this_run=DFTRun(basename)
    this_run.status = 1
    ### test if outfile exits
    this_run.outpath = (path_dictionary["out_path" ] + "/gen_" + str(gen) +"/"
                           + basename + ".out")
    ## load details into run
    this_run.configure(metal,ox,eqlig,axlig1,axlig2,spin,alpha,spin_cat)
    this_run.gene =  gene
    if os.path.exists(this_run.outpath):
        ### file is found,d check if converged
        with open(this_run.outpath) as f: 
            data=f.readlines()
            found_conv =False 
            found_data =False
            found_time = False 
            for i,lines in enumerate(data):
                if str(lines).find('Running Mulliken') != -1:
                    found_conv = True
                if str(lines).find('FINAL ENERGY') != -1:
                    this_run.energy =str(lines.split()[2])
                    found_data = True
                if str(lines).find('Total processing time') != -1:
                    this_run.time=str(lines.split()[3])
                    found_time = True
                if str(lines).find('SPIN S-SQUARED') != -1:
                    this_str=(lines.split())
                    this_run.ss_act =float( this_str[2])
                    this_run.ss_target = float(this_str[4].strip('()'))
        if (found_data == True) and (found_time == True) and (found_conv == True):
            this_run.converged = True
            this_run.status = 0
    return this_run

def process_runs_sp(LS_runs,HS_runs):
    ## function to compare HS/LS
    ##  runs and calc splitting energies
    #  @param LS_runs low spin runs
    #  @param HS_runs high spin runs
    #  @return final_results dictionary 
    #          of comparison classes
    ## this is for GA only 
    final_results=dict()
    matched = False
    number_of_matches  = 0
    for genes in LS_runs.keys():
        matched = 0 
        LS_run = LS_runs[genes]
        this_name = LS_run.name
        this_gene = genes
        if this_gene in HS_runs.keys():
            HS_run = HS_runs[this_gene]
            matched = True
            number_of_matches += 1
        if matched:
            print('matched ID: '+ str(this_gene) + ' files ' + str(HS_run.name) + ' and ' + str(LS_run.name))
            final_results[this_gene] = Comp(this_gene)
            final_results[this_gene].gene = this_gene
            final_results[this_gene].set_properties(LS_run)
            final_results[this_gene].LS_energy = str(float(LS_run.energy))
            final_results[this_gene].HS_energy = str(float(HS_run.energy))
            final_results[this_gene].LS_time = str(float(LS_run.time))
            final_results[this_gene].HS_time = str(float(HS_run.time))
            final_results[this_gene].HS_status = HS_run.status
            final_results[this_gene].LS_status = LS_run.status
            final_results[this_gene].HS_time = HS_run.time
            final_results[this_gene].HS_time = HS_run.time
            final_results[this_gene].process()
            final_results[this_gene].HS_ss_act = HS_run.ss_act
            final_results[this_gene].LS_ss_act = LS_run.ss_act
            final_results[this_gene].LS_ss_target = LS_run.ss_target
            final_results[this_gene].HS_ss_target = HS_run.ss_target
            final_results[this_gene].max_spin_error =max(abs( float(HS_run.ss_target) - float(HS_run.ss_act)),abs(float(LS_run.ss_target - LS_run.ss_act)))
        else:
            print('unmatched ID: '+ str(this_gene) + ' files ' + str(LS_run.name)+ ' has no partner' )
    return final_results
def process_runs_geo(all_runs,list_of_prop_names,local_spin_dictionary,local_metal_list=False):
    ## function to find mathcing runs by gene
    ## and extract their properties
    ##  for terachem GO runs
    #  @param all_runs list of runs
    #  @param list_of_prop_names list of properties
    #                            to carry over
    #  @param local_spin_dictionary metals and spin states
    #                               used to define low and
    #                               high spins expected
    #  @return final_results dictionary of comparisons keyed by gene
    final_results=dict()
    if not local_metal_list:
        local_metal_list = get_metals()
    matched = False
    number_of_matches  = 0
    print(local_spin_dictionary)
    print('processing all converged runs')
    for runkeys in all_runs.keys():
        skip = False
        duplication = False
        this_run = all_runs[runkeys]
        if this_run.metal in local_metal_list:
            this_metal = this_run.metal
        else:
            this_metal = local_metal_list[int(this_run.metal)]

        this_name = "_".join([this_metal,'eq',str(this_run.eqlig),'ax1',str(this_run.axlig1),'ax2',str(this_run.axlig2),'ahf',str(this_run.alpha)])
                ### add alpha value to list owned by this_comp:
                
        if this_name not in final_results.keys():
            ## need to create a new holder to store this gene
            this_comp = Comp(this_name)
            this_comp.set_properties(this_run)
            
        else:
            this_comp = final_results[this_name]
        print(runkeys)
        this_comp.attempted += 1 # advance number of attempts
        ## get basic details
        this_ox = int(this_run.ox)
        metal_spins  = local_spin_dictionary[this_metal][this_ox]
        if this_run.spin not in metal_spins:
           print('ERROR! not in metal spins : ' +  str(this_run) + ' not in ' +  str(metal_spins))
        else:
            spin_ind = metal_spins.index(this_run.spin)
            if spin_ind == 0:
                 spin_cat = 'LS'
            else:
                spin_cat = 'HS'
        print('spin ind is found to be ' + str(this_run.spin) + ' interpretted as ' + str(spin_cat))
        ## check if this a duplicate:
        this_attribute = "_".join(['ox',str(this_ox),spin_cat,'converged'])
        if getattr(this_comp,this_attribute):
            duplication = True
            print('run duplication at  ' +str(this_name))
            if this_ox == 2:
                this_ox2RN = this_run.number
                old_ox2RN = this_comp.ox2RN
                if this_ox2RN <= old_ox2RN:
                    skip = True
                    ## use this one, get rid of the old one 
            else:
                this_ox3RN = this_run.number
                old_ox3RN = this_comp.ox3RN
                if this_ox3RN <= old_ox3RN:
                    skip = True
        if not skip:
            if duplication:
                # set back conv counter
                this_comp.convergence -= 1
                if this_run.coord == 6:
                    ## replace the descriptor if set
                    this_comp.set_desc = False
            
            ## find oxidation state:
            if this_ox == 2:
                 this_comp.ox2RN = max(this_run.number,this_comp.ox2RN)
            else:
                 this_comp.ox3RN = max(this_run.number,this_comp.ox3RN)
            if this_comp.gene =="undef":
                this_comp.gene = this_run.gene
            if this_run.converged and this_run.coord == 6:
                this_comp.convergence += 1
            if this_run.coord == 6 and not this_comp.set_desc:
                try:
                    if not os.path.isdir('used_geos/'):
                        os.mkdir('used_geos/')
                    this_run.mol.writexyz('used_geos/'+this_name+'.xyz')
                    this_comp.axlig1 = this_run.axlig1
                    this_comp.axlig2 = this_run.axlig2
                    this_comp.eqlig = this_run.eqlig
                    this_comp.set_rep_mol(this_run)
                    this_comp.get_descriptor_vector(loud=False,name=this_name)
                except:
                    if not os.path.isdir('bad_geos/'):
                        os.mkdir('bad_geos/')
                    this_run.mol.writexyz('bad_geos/'+this_name+'.xyz')
                    this_comp.convergence -= 1
                    this_run.coord = 'error'
#            print(dir(this_run))
            for props in list_of_prop_names:
                     this_attribute = "_".join(['ox',str(this_ox),spin_cat,props])
#                     print(this_attribute)
                     setattr(this_comp,this_attribute,getattr(this_run,props))
#            if this_run.coord == 6 and spin_cat == 'HS' and this_ox == 2:
#                if not os.path.isdir('coulomb_geos/'):
#                        os.mkdir('coulomb_geos/')
#                this_run.mol.writexyz('coulomb_geos/'+this_name+'.xyz')
#                this_comp.get_coulomb_descriptor(size=85)
        ## the hack to get around expecting 
        ## spins
        this_comp.get_some_split()
        ###
        final_results.update({this_name:this_comp})
    return final_results
    
def check_solvent_file(this_run):
    ## function to test solvent single point convergence
    ##  for terachem files
    #  @param this_run a run class
    #  @return this_run populated run class
	solvent_contribution = False
	found_solvent_energy = False
	found_solvent_cont = False
	if os.path.exists(this_run.solvent_outpath):
		### file is found, check if converged
		with open(this_run.solvent_outpath) as f:
			sol_data=f.readlines()
			for i,lines in enumerate(sol_data):
				if str(lines).find('FINAL ENERGY') != -1:
					found_solvent_energy = True
					print('found solvent energy ' + this_run.name)
				if str(lines).find('C-PCM contribution ') != -1:
						solvent_contribution =str(lines.split()[4]).split(':')[1]
						found_solvent_cont = True
						print('found solvent contri ' + this_run.name)
        			if str(lines).find('Total processing time') != -1:
	        			this_run.solvent_time = str(lines.split()[3])

        if (found_solvent_energy == True) and (found_solvent_cont == True):
			this_run.solvent_cont = solvent_contribution
	return(this_run)
    
def check_thermo_file(this_run):
    ## function to test thermodynamic contribution
    ##  for terachem files
    #  @param this_run a run class
    #  @return this_run populated run class
	found_vib_correction = False
	found_grad_error = False
	vib_correction = False
	if os.path.exists(this_run.thermo_outpath):
		with open(this_run.thermo_outpath) as f:
			thermo_data=f.readlines()
			for i,lines in enumerate(thermo_data):
				if str(lines).find('Total Vibrational Energy Correction') != -1:
					vib_correction =str(lines.split()[5])
					found_vib_correction = True
				        print('found vib correction ' + this_run.name)
			        if str(lines).find('imaginary frequencies') != -1:
				        this_run.imag =str(lines.split()[0])
		        		print('found imag ' + this_run.name)
			        	print(lines)
        			if str(lines).find('Maximum component of gradient is too large') != -1:
	        			this_run.thermo_status ="GRAD_TOO_LARGE"
		        		found_grad_error = True
			        	print('found GRAD error ' + this_run.name)
				        print(lines)
        			if str(lines).find('Total processing time') != -1:
	        			this_run.thermo_time = str(lines.split()[3])

		if (found_vib_correction == True) and (found_grad_error == False):
			this_run.thermo_cont = vib_correction
		if (found_vib_correction == True) and (found_grad_error == True):
			this_run.thermo_cont = "grad_error"
			this_run.comment +="grad_error\n"

            

	return(this_run)

def check_init_sp(this_run):
    ## function to test initial single point convergence
    ##  for terachem files
    #  @param this_run a run class
    #  @return this_run populated run class
    found_data = False
    if os.path.exists(this_run.init_outpath):
        with open(this_run.outpath) as f: 
            data=f.readlines()
            found_conv =False 
            found_data =False
            found_time = False 
            for i,lines in enumerate(data):
                if str(lines).find('FINAL ENERGY') != -1:
                    energy =str(lines.split()[2])
                    found_data = True
        if (found_data == True):
			this_run.init_energy = energy
    return this_run
def check_mopac(this_run):
    ## function to test mopac convergence
    #  @param this_run a run class
    #  @return this_run populated run class
    print(this_run.moppath)
    if os.path.exists(this_run.moppath):
        print('mopac file exists')
        with open(this_run.moppath) as f:
            conv_flag =  False
            data=f.readlines()
            found_geo = False
            this_geo = list()
            mop_converged = False
            in_cord = False
            for i,lines in enumerate(data):
                if str(lines).find('TOTAL ENERGY') != -1:
                    this_run.mop_energy =str(float(lines.split()[3])*eV_to_Kcal_mol)
                if str(lines).find('Converged!') != -1:
                    unconv = 0
                if  (str(lines).find('SCF FIELD WAS ACHIEVED ') != -1):
                    conv_flag =  True
                    mop_converged = True
                if (str(lines).find('CARTESIAN COORDINATES') != -1) and  (conv_flag):
                    in_cord = True
                    print('found mopac geo')
                if in_cord:
                    if (str(lines).find('Empirical Formula') != -1):
                        in_cord =  False
                        found_geo = True
                        print('final line of mopac of geo')
                    else:
                        if lines.strip():
                            this_geo.append(list(lines[1:]))
        if mop_converged and found_geo:
            with open(this_run.mop_geopath,'w') as f:
                f.write(str(int(len(this_geo))-1)+'\n')
                f.write('#'+this_run.name+' mopac \n')
                for i,elements in enumerate(this_geo):
                    if not i==0:
                        line_tw = ''.join(elements[5:])
                        line_tw= line_tw.lstrip()
                        f.write(line_tw)
            this_run.obtain_mopac_mol()
            this_run.check_coordination()
        return(this_run)
            



def test_terachem_go_convergence(this_run):
    ## function to test geometry optimization convergence
    ##  for terachem files
    #  @param this_run a run class
    #  @return this_run populated run class    
    print('we have access go to test_terachem_go function' )
    if not this_run.logpath:
        this_run.logpath = get_run_dir()
    print('logging to ' +this_run.logpath)
    if os.path.exists(this_run.geopath):
        this_run.geo_exists = True
        print('geo exists ' +this_run.geopath)
    else:
        this_run.comment += 'no geo found\n'
        print(' no  geo exists ' +this_run.geopath)
        if os.path.exists(this_run.scrpath):
            this_run.extract_geo()
            print('  geo extracted to  ' +this_run.geopath)
        else:
            print(' cannot find scr:   ' +this_run.scrpath)

    if os.path.exists(this_run.outpath):
        read_terachem_go_output(this_run)
    else:
        this_run.comment += ' no outfile found\n'
    print('has run converged : ' + str(this_run.converged))
    if this_run.converged:
        logger(this_run.logpath, str(this_run.name) + ' run converged ' +  ' and now testing geo '+this_run.geopath )
        # check the geo
        if os.path.exists(this_run.geopath):

                print(this_run.geopath + ' found')
                # get mol3D file
                this_run.obtain_mol3d()
                # check if inidicators are good
                # check coordinattion
                this_run.check_coordination()
                logger(this_run.logpath, str(this_run.name) + ' cooridination is ' +str(this_run.coord))
                print(str(this_run.name) + ' cooridination is ' +str(this_run.coord))
                # ML distances, geo
                this_run.obtain_ML_dists()
                
                ## check intial conditions:
                if os.path.exists(this_run.init_geopath):
                    this_run.obtain_init_mol3d()
                    this_run.check_oct_needs_init()
                    this_run.obtain_rsmd() # copmare to initial
                else:
                    this_run.check_oct_needs_final_only()
        print('this flag oct is '+ str(this_run.flag_oct))
        if this_run.coord == 6 and this_run.converged and this_run.flag_oct == 1:
            this_run.status = 0
            if not this_run.tspin == this_run.spin:
                print(this_run.tspin)
                print(this_run.spin)
                # sardines
                
        else:
            this_run.status = 1
            this_run.comment += 'coord not good ' +str(this_run.coord) +'\n '
            this_run.comment += 'flag_oct_list: %s\n'%(this_run.flag_oct_list)
            

def read_terachem_go_output(this_run):
    ## function to parse geometry optimization outfile
    ##  for terachem files
    #  @param this_run a run class
    #  @return this_run populated run class  
    found_conv =False 
    found_data =False
    converged = False
    found_time = False 
    found_init = False
    if os.path.exists(this_run.outpath):
        ### file is found, check if converged
        with open(this_run.outpath) as f:
            this_run.attempted  = True
            data=f.readlines()
            for i,lines in enumerate(data):
                if str(lines).find('TeraChem v') != -1:
                    this_run.terachem_version = lines.split()[2]
                    #print('TeraChem Version: ' + this_run.terachem_version)
                if str(lines).find('Hg Version') != -1:
                    this_run.terachem_detailed_version = lines.split()[3]
                    #print('TeraChem Hg build: ' + this_run.terachem_detailed_version)
                if str(lines).find('Using basis set') != -1:
                    this_run.basis = lines.split()[3]
                    #print('TeraChem basis: ' + this_run.basis)
                if str(lines).find('Spin multiplicity') != -1:
                    this_run.tspin = int(lines.split()[2])
                    #print('TeraChem spin: ' + str(this_run.tspin))
                if str(lines).find('Total charge:') != -1:
                    this_run.charge = int(lines.split()[2])
                    #print('TeraChem charge: ' + str(this_run.charge))
                if str(lines).find('Alpha level shift') != -1:
                    this_run.alpha_level_shift = float(lines.split()[3])
                    #print('Alpha level: ' + str(this_run.alpha_level_shift)  )
                if str(lines).find('Beta level shift') != -1:
                    this_run.beta_level_shift = float(lines.split()[3])
                if str(lines).find('DFT Functional requested:') != -1:
                    this_run.functional = lines.split()[3]                    
                    #print('TC functional: ' + this_run.functional  )
                if (str(lines).find('Optimization Converged.') != -1) or (str(lines).find('Converged!') != -1):
                   found_conv = True
                if str(lines).find('FINAL ENERGY') != -1:
                    this_run.energy =str(lines.split()[2])
                    found_data = True
                    if not found_init:
                        this_run.init_energy = str(lines.split()[2])
                        found_init = True
                if str(lines).find('Total processing time') != -1:
                    this_run.time=str(lines.split()[3])
                    found_time = True
                if str(lines).find('SPIN S-SQUARED') != -1:
                    this_str=(lines.split())
                    this_run.ss_act =float( this_str[2])
                    this_run.ss_target = float(this_str[4].strip('()'))
    else:
        this_run.attempted  = False
    if (found_data == True) and (found_time == True) and (found_conv == True):
        this_run.converged = True
        print('run outfile converged')
        
        
