import glob
import glob
import math
import numpy
import subprocess
import argparse
import os
import random
import shutil
import sys

from ga_tools import *

class octahedral_complex:
    def __init__(self,ligands_list):
        self.free_sites = [1,2,3,4,5,6]
        ### mark 3x bidentate
        self.three_bidentate = False
        self.ready_for_assembly = False
        self.has_vacant = False
        self.ligands_list=get_ligands()
        self.metals_list= get_metals()
        self.ox_list = get_ox_states()
        self.core  = 'U'
        self.ox = 'U'
        self.ax_dent =  False
        self.eq_dent = False
        self.eq_ligands = list()
        self.ax_ligands= list()
        self.ax_inds = list()
        self.eq_inds = list()
        GA_run = get_current_GA()
        self.ahf = int(GA_run.config["exchange"]) # % HFX, B3LYP = 20 
    def random_gen(self):
        self._get_random_metal()
        self._get_random_ox()
        self._get_random_equitorial()
        self._get_random_axial()
        self._name_self()
    def _name_self(self):
        ## this creates the gene for a given complex
        GA_run = get_current_GA()
        self.name = "_".join([str(self.core),str(self.ox),str(self.eq_inds[0]),str(self.ax_inds[0]),str(self.ax_inds[1]),str(self.ahf)])

    def copy(self,partner):
         self.core = partner.core
         self.ox = partner.ox
         self.ax_dent = partner.ax_dent
         self.ax_inds = partner.ax_inds
         self.ax_ligands = partner.ax_ligands
         self.eq_dent = partner.eq_dent
         self.eq_inds = partner.eq_inds
         self.eq_oc = partner.eq_oc
         self.eq_ligands = partner.eq_ligands
         self.three_bidentate = partner.three_bidentate
         self.ahf = partner.ahf
         self._name_self()
    def _get_random_metal(self):
        n = len(self.metals_list)
        metal_ind = numpy.random.randint(low = 0,high = n)
        self.core = metal_ind
    def _get_random_ox(self):
        possible_ox_states = get_ox_states()
        self.ox = numpy.random.choice(possible_ox_states)
    def _get_random_equitorial(self):
        ### choose the equitorial ligand
        n = len(self.ligands_list)
        eq_ind = numpy.random.randint(low = 0,high = n)
        #print('choosing '+str(eq_ind) + ' from '+str(n))
        # now we test if it is a SMILEs or molsimplify ligand    
        #print('name  is ' +str(self.ligands_list[eq_ind][0]))
        eq_ligand_properties  = self.ligands_list[eq_ind][1]
        self.eq_dent = eq_ligand_properties[0]
        self.eq_oc  = int(4/self.eq_dent)
        self.eq_ligands = [self.ligands_list[eq_ind][0] for i in range(0,self.eq_oc)]
        self.eq_inds = [eq_ind]
        
    def _get_random_axial(self):
        ### choose axial ligands:
        n = len(self.ligands_list)
        self.ready_for_assembly =  False
        self.ax_ligands = list()
        self.ax_inds = list()
        while not self.ready_for_assembly:
            ax_ind = numpy.random.randint(low = 0,high = n)
            ax_ligand_properties  = self.ligands_list[ax_ind][1]
            ax_dent = ax_ligand_properties[0]
            if ax_dent > 1:
                if ((self.eq_dent == 2) and (ax_dent == 2) and (len(self.ax_ligands) == 0)):
                    three_bidentate  = True
                    self.ax_ligands = [self.ligands_list[ax_ind][0],self.ligands_list[ax_ind][0]]
                    self.ax_inds = [ax_ind, ax_ind]
                    self.ready_for_assembly = True
                    self.ax_dent = 2
                    self.ax_oc = [2]
                else:
                    self.ready_for_assembly = False
            elif ax_dent == 1:
                GA_run = get_current_GA()
                if GA_run.config["symclass"] == "strong":
                    self.ax_ligands = [self.ligands_list[ax_ind][0],self.ligands_list[ax_ind][0]]
                    self.ax_inds = [ax_ind, ax_ind]
                    if (len(self.ax_ligands) ==2):
                        self.ax_dent = 1
                        self.ax_oc = [1,1]
                        self.ready_for_assembly = True
                else:
#                    This section is intended to allow for vacant axial sites
#                    It is not currently implemented               
#                    if ((ax_ind == 0) and (not has_zero)):
#                        has_zero = True
#                        self.ax_ligands.append(self.ligands_dict[ax_ind][0])
#                        self.ax_inds.append(ax_ind)
#                        self.geo = "spy"
#         
#                    elif ax_ind !=  0:
                    self.ax_ligands.append(self.ligands_list[ax_ind][0])
                    self.ax_inds.append(ax_ind)
                    if (len(self.ax_ligands) ==2):
                        self.ax_dent = 1
                        self.ax_oc = [1,1]
                        self.ready_for_assembly = True
                    else:
                         self.ready_for_assembly = False

    def examine(self):
        print("name is " + self.name)
        print("eq", self.eq_ligands, self.eq_inds)
        print("axial", self.ax_ligands, self.ax_inds)

    def encode(self,gene):
        self.random_gen()
        ll = gene.split("_")
        ll = [int(item) for item in ll]
        #print('ll is '+str(ll))
        self.replace_metal(ll[0])
        self.replace_ox(ll[1])
        self.replace_equitorial([ll[2]])
        self.replace_axial(ll[3:5])
        self.ahf=ll[5]
        self._name_self()
    def replace_metal(self,new_metal_ind):
        self.core = new_metal_ind
    def replace_ox(self,new_ox):
        self.ox = new_ox
    def replace_equitorial(self,new_eq_ind):
        #print('in repcoding, setting eq to ' + str(new_eq_ind))
        
        eq_ligand_properties  = self.ligands_list[new_eq_ind[0]][1]
        self.eq_dent = eq_ligand_properties[0]
        self.eq_oc  = int(4/self.eq_dent)
        self.eq_ligands = [self.ligands_list[new_eq_ind[0]][0] for i in range(0,self.eq_oc)]
        self.eq_inds = new_eq_ind
        if (self.ax_dent == 1) or ((self.ax_dent == 2) and (self.eq_dent ==2)):
                ## everything is ok!
                if (self.ax_dent == 2):
                    self.three_bidentate = True
        else: ## this complex cannot exist. keeping  equitorial,
              ## regenerating axial ligands
           #print("complex with" + str(self.eq_ligands) + " and " + str(self.ax_ligands) + " cannot exist, randomizing axial")
           self._get_random_axial()
        self._name_self()

    def replace_axial(self,new_ax_ind):
        self.ax_ligands = list()
        self.ax_inds = list()
        n = len(self.ligands_list)
        #print(new_ax_ind)
        self.three_bidentate =  False
        for i,indices in enumerate(new_ax_ind):
            #print('trying to add ',indices )
            ax_ligand_properties = self.ligands_list[indices][1]
            ax_dent = ax_ligand_properties[0]
            #print('ax dent ' + str(ax_dent))
            if (ax_dent > 1):
                if (ax_dent == 2) and (i == 0):
                    three_bidentate  = True
                    self.ax_ligands = [self.ligands_list[indices][0],self.ligands_list[indices][0]]
                    self.ax_inds = [indices, indices]
                    self.ax_dent = 2
                    self.three_bidentate = True
                    break
                else:
                    print('impossible, ax_dent  = ' + str(ax_dent))
            else:
                self.ax_ligands.append(self.ligands_list[indices][0])
                self.ax_inds.append(indices)
                self.ax_dent = 1

        if (self.three_bidentate) and not (self.eq_dent == 2):
            ## this complex cannot exist
            ## regenerating equitorial ligands
            #print("complex with" + str(self.eq_ligands) + " and " + str(self.ax_ligands) +
            #      " cannot exist, regenerating equitorial ligands")
            self.ready_for_assembly  = False
            while not self.ready_for_assembly:
                eq_ind = numpy.random.randint(low = 0,high = n)
                eq_ligand_properties  = self.ligands_list[eq_ind][1]
                eq_dent = eq_ligand_properties[0]
                if (eq_dent == 2):
                    self.eq_dent = 2
                    self.eq_inds = [eq_ind]
                    self.eq_oc = int(4/eq_dent)
                    self.eq_ligands = [self.ligands_list[eq_ind][0] for i in range(0,self.eq_oc)]
                    self.ready_for_assembly =  True
        self._name_self()

    def exchange_ligands(self,partner,eq_swap):
        child = octahedral_complex(self.ligands_list)
        child.copy(self) # copies this parent
        
        print("swapping from",partner.name," to ",self.name)
        self.examine()
        if eq_swap:
            print("swapping equitorial " + str(child.eq_inds) + ' -> ' + str(partner.eq_inds))
            child.replace_equitorial(partner.eq_inds)
        else:
            print("swapping axial"+ str(child.ax_inds) + ' -> ' + str(partner.ax_inds))
            child.replace_axial(partner.ax_inds)
        child.examine()
        child._name_self()
        return child
        
    def exchange_metal(self,partner):
        child = octahedral_complex(self.ligands_list)
        child.copy(self) # copies this parent
        print("swapping from",partner.name," to ",self.name)
        self.examine()
        child.replace_metal(partner.core)
        child._name_self()
        return child
        
    def exchange_ox(self,partner):
        child = octahedral_complex(self.ligands_list)
        child.copy(self) # copies this parent
        print("swapping from",partner.name," to ",self.name)
        self.examine()
        child.replace_ox(partner.ox)
        child._name_self()
        return child

    def mutate(self):
        ## mutates either the axial
        ## or equitorial ligand a random
        GA_run = get_current_GA()
        lig_to_mutate = random.randint(0,3)
        child = octahedral_complex(self.ligands_list)
        child.copy(self) # copies this parent
        n = len(self.ligands_list)
        self.examine()
        print('I think this is 3x bidentate: ',self.three_bidentate,self.ax_dent)
        if (lig_to_mutate == 0):
            print("mutating equitorial ligand")
            rand_ind = numpy.random.randint(low = 0,high = n)
            child.replace_equitorial([rand_ind])
        elif (lig_to_mutate == 1):
            print('mutating axial ligand')
            ready_flag = False
            while not ready_flag:
                new_ax_list = list()
                rand_ind = numpy.random.randint(low = 0,high = n)
                ax_ligand_properties  = self.ligands_list[rand_ind][1]
                ax_dent = ax_ligand_properties[0]
                if (ax_dent == self.ax_dent):
                    if (lig_to_mutate == 1):
                        print("mutating axial 1 ")
                        if GA_run.config['symclass'] =="strong":
                            new_ax_list = [rand_ind,rand_ind]
                        else:
                            new_ax_list = [rand_ind,self.ax_inds[1]]
                    elif (lig_to_mutate == 2):
                        print("mutating axial 2 ")
                        if GA_run.config['symclass'] =="strong":
                            new_ax_list = [rand_ind,rand_ind]
                        else:
                            new_ax_list = [self.ax_inds[0],rand_ind]
                        
                    child.ax_dent = 1
                    child.three_bidentate = False
                    ready_flag = True
                elif (ax_dent  == 2) and (self.ax_dent == 1):
                    ## here, we want to add a bidentate but
                    ## need to swap the second ligand too
                    print("swapping both axial ")
                    new_ax_list = [rand_ind,rand_ind]
                    child.ax_dent = 2
                    child.three_bidentate = True
                    ready_flag =  True
            print("trying to add " + str(new_ax_list))
            child.replace_axial(new_ax_list)
        elif (lig_to_mutate == 3): ## metal mutation
            print('mutating metal')
            child._get_random_metal()
            print('mutating ox')
            child._get_random_ox()
        child._name_self()
        child.examine()
        return child

    def generate_geometery(self,prefix,spin,path_dictionary,rundirpath,gen):
        # get path info
        ligloc_cont = True
        # set metal properties:
        this_metal = self.metals_list[self.core]
        # set oxidation state
        if self.ox == 2:
            ox_string = "II"
        elif self.ox == 3:
            ox_string = "III"
        mol_name = prefix + self.name + "_" + str(spin)
        
        smicat = False # holder for connection atoms calls for SMILES ligands
        if self.ax_dent == 1:
            # assemble SMILEs ligands
            liglist = "" #empty string
            for eq_lig in self.eq_ligands:
                #print('processing ligand ' + str(eq_lig))
                if not hasattr(eq_lig,'__iter__'): # test if SMILES:
                    #print('molSimplify ligand in eq position  '+ str(eq_lig))
                    liglist += " " + str(eq_lig).strip("'[]'")
                elif  hasattr(eq_lig,'__iter__'): # this is the mark of SMILES strings:
                    #print('SMILEs ligand in eq position  '+ str(eq_lig))
                    liglist += " " +  "'"+str(eq_lig[0])+ "'"
                    if not smicat: # false on first hit 
                        smicat = " ["    + str(eq_lig[1]).replace("'","") # cat list 
                    else:
                        smicat += ",  " + str(eq_lig[1]).replace("'","") # cat list 
            #print(' after eq-shell, liglist is ' + liglist)
            #print(' after eq-shell, smicat is ' + str(smicat))
            for ax_lig in self.ax_ligands:
                #print('processing ligand ' + str(ax_lig))
                if not hasattr(ax_lig,'__iter__'): # test if SMILES:
                    #print('molSimplify ligand in ax position  '+ str(ax_lig))
                    liglist += " " + str(ax_lig).strip("'[]'")
                elif  hasattr(ax_lig,'__iter__'): # this is the mark of SMILES strings:
                    #print('SMILEs ligand in ax position  '+ str(ax_lig))
                    liglist += " " +  "'"+ str(ax_lig[0])+  "'"
                    if not smicat: # false on first hit 
                        smicat = " [ "    + str(ax_lig[1]).replace("'","") # cat list 
                    else:
                        smicat += ",  " + str(ax_lig[1]).replace("'","")  # cat list 
            #liglist.replace("'","")
            if smicat:
                smicat += "]"
            #print(' after ax-shell, liglist is ' + liglist)
            #print(' after ax-shell, smicat is ' + str(smicat))            
            ligloc = 1
            ligalign = 0
        elif self.ax_dent == 2:
            half = len(self.eq_ligands)/2
            liglist = "" #empty string
            # first eq ligand 
            eq_lig = self.eq_ligands[0]
            #print('processing ligand ' + str(eq_lig))
            if not hasattr(eq_lig,'__iter__'): # test if SMILES:
                #print('molSimplify ligand in eq position  '+ str(eq_lig))
                liglist += " " + str(eq_lig).strip("'[]'")
            elif  hasattr(eq_lig,'__iter__'): # this is the mark of SMILES strings:
                #print('SMILEs ligand in eq position  '+ str(eq_lig))
                liglist += " " + str(eq_lig[0])
                if not smicat: # false on first hit 
                    smicat = " ["    + str(eq_lig[1]).replace("'","") # cat list 
                else:
                    smicat += ",  " + str(eq_lig[1]).replace("'","") # cat list 
            #print(' after eq-shell 1, liglist is ' + liglist)
            #print(' after eq-shell 1, smicat is ' + str(smicat))
            
            # only one axial ligand allowed
            ax_lig = self.ax_ligands[0]
            #print('processing ligand ' + str(ax_lig))
            if not hasattr(ax_lig,'__iter__'): # test if SMILES:
                #print('molSimplify ligand in ax position  '+ str(ax_lig))
                liglist += " " + str(ax_lig).strip("'[]'")
            elif  hasattr(ax_lig,'__iter__'): # this is the mark of SMILES strings:
                #print('SMILEs ligand in ax position  '+ str(ax_lig))
                liglist += " " + str(ax_lig[0])
                if not smicat: # false on first hit 
                    smicat = " ["    + str(ax_lig[1]).replace("'","") # cat list 
                else:
                    smicat += ",  " + str(ax_lig[1]).replace("'","")  # cat list 
            

            #print(' after ax-shell, liglist is ' + liglist)
            #print(' after ax-shell, smicat is ' + str(smicat))            
            
            # second eq ligand 
            eq_lig = self.eq_ligands[1]
            #print('processing ligand ' + str(eq_lig))
            if not hasattr(eq_lig,'__iter__'): # test if SMILES:
                #print('molSimplify ligand in eq position  '+ str(eq_lig))
                liglist += " " + str(eq_lig).strip("'[]'")
            elif  hasattr(eq_lig,'__iter__'): # this is the mark of SMILES strings:
                #print('SMILEs ligand in eq position  '+ str(eq_lig))
                liglist += " " + str(eq_lig[0])
                if not smicat: # false on first hit 
                    smicat = " ["    + str(eq_lig[1]).replace("'","") # cat list 
                else:
                    smicat += ",  " + str(eq_lig[1]).replace("'","") # cat list 
            #print(' after eq-shell 2, liglist is ' + liglist)
            #print(' after eq-shell 2, smicat is ' + str(smicat))
            liglist.replace("'","")
            if smicat:
                smicat += "]"            
#            liglist = (str([str(element).strip("'[]'") for element in (self.eq_ligands[:half])] +
#             [str(self.ax_ligands[0]).strip("'[]'")] +
#              [str(element).strip("'[]'") for element in (self.eq_ligands[half:])]).strip("[]")).replace("'", "")
            ligloc = 'false'
            ligalign = 0
        ## disable force field opt
        ## if not DFT
        if isDFT():
            ff_opt = 'A'
        else: # optional denticity based FF control
            if max([self.ax_dent,self.eq_dent]) ==1:
                ff_opt = 'A'
            else:
                ff_opt = 'A'

        ## get custom exchange fraction
        this_GA = get_current_GA()
        exchange = this_GA.config['exchange']
        optimize = this_GA.config['optimize']
        
        if optimize:
            print(' setting up GEO optimization ')
            rty = 'minimize'
            scrpath =  "scr/geo/gen_"+str(gen) + '/'  + mol_name
        else:
            rty = 'energy'
            scrpath =  "scr/sp/gen_"+str(gen) + '/'  + mol_name
            
        geometry = "oct"
        
        ## set paths for generation
        ms_dump_path = path_dictionary["molsimplify_inps"] +  'ms_output.txt'
        ms_error_path = path_dictionary["molsimplify_inps"] +  'ms_errors.txt'
        jobpath = path_dictionary["job_path"]  + mol_name + '.in'
        inpath = path_dictionary["infiles"]  + mol_name + '.in'
        
        

        geometry_path = path_dictionary["initial_geo_path"] +'/'+ mol_name + '.xyz'
        
        ## check if already exists:
        geo_exists = os.path.isfile(path_dictionary["initial_geo_path"] + mol_name + '.xyz')

        if not (geo_exists):
                print('generating '+ str(mol_name) + ' with ligands ' + str(self.eq_ligands) + ' and'  + str(self.ax_ligands))
                try:
                #if True:
                    with open(ms_dump_path,'a') as ms_pipe:
                        with open(ms_error_path,'a') as ms_error_pipe:
                            call = " ".join(["molsimplify " ,'-core ' + this_metal,'-lig ' +liglist,'-ligocc 1,1,1,1,1,1',
                                     '-rundir ' +"'"+ rundirpath.rstrip("/")+"'",'-keepHs yes,yes,yes,yes,yes,yes','-jobdir','temp',
                                     '-coord 6','-ligalign '+str(ligalign),'-ligloc ' + str(ligloc),'-calccharge yes','-name '+"'"+mol_name+"'",
                                     '-geometry ' + geometry,'-spin ' + str(spin),'-oxstate '+ ox_string, '-exchange '+str(exchange),
                                     '-qccode TeraChem','-runtyp '+rty,'-method UDFT',"-ffoption "+ff_opt,' -ff UFF'])
                            if smicat:
                                call += ' -smicat ' + smicat
                            print(call)
                            p2 = subprocess.call(call,stdout = ms_pipe,stderr=ms_error_pipe, shell=True)
                    assert(os.path.isfile(rundirpath + 'temp'+'/' + mol_name + '.molinp'))
                    shutil.move(rundirpath + 'temp'+'/' + mol_name + '.molinp', path_dictionary["molsimplify_inps"]+'/' + mol_name + '.molinp')
                    shutil.move(rundirpath + 'temp'+'/' + mol_name + '.xyz', geometry_path)
                except:
                        print('Error: molSimplify failure when calling ')
                        print(call)
                        sys.exit()
                    
                if this_GA.config['symclass']=="strong":                        
                    with open(rundirpath + 'temp' +'/' + mol_name + '.report') as report_f:
                        for line in report_f:
                                if ("pred_split" in line):
                                    #print('****')
                                    #print(line)
                                    ANN_split = float(line.split(",")[1])
                                    print('ANN_split is ' +"{0:.2f}".format(ANN_split))
                                if("ANN_dist_to_train" in line):
                                    #print('****')
                                    #print(line)
                                    ll = line.split(',')[1]
                                    ANN_distance = float(ll)
                                    print('ANN_distance is ' +"{0:.2f}".format(ANN_distance))
                elif this_GA.config['symclass']=="weak": ## ANN not currently supported!
                    ANN_split = False
                    ANN_distance = False
                       
                
                shutil.move(rundirpath + 'temp' +'/' + mol_name + '.report', path_dictionary["ms_reps"] +'/'+ mol_name + '.report')

                ## write the job file
                with open(jobpath,'w') as newf:
                    with open(rundirpath + 'temp/' + mol_name + '.in','r') as oldf:
                        for line in oldf:
                            if not ("coordinates" in line) and (not "end" in line) and not ("scrdir" in line):
                                newf.writelines(line)
                    newf.writelines("scrdir " +scrpath + "\n")
                os.remove(rundirpath + 'temp/' + mol_name + '.in')
                
                ### make an infile!
                create_generic_infile(jobpath,restart=False)
                    
                
        else:
            
            ANN_split = False
            ANN_distance = False
        if not 'ANN_split' in dir():
            ANN_split = False
            ANN_distance = False
        print('!!!ANN_split:', ANN_split)
        return jobpath,mol_name,ANN_split,ANN_distance


