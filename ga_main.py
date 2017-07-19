import glob
import operator
import datetime
import math
import numpy
import subprocess
import argparse
import os
import random
import shutil
from ga_tools import *
from tree_classes import *
from ga_check_jobs import *

########################
def update_current_gf_dictionary(gene,fitness):
     ## set up environment:        
     path_dictionary = setup_paths()
     new_tree = tree_generation('temp tree')
     ## read in info
     new_tree.read_state()
     new_tree.gene_fitness_dictionary.update({gene:fitness})
     logger(path_dictionary['state_path'],str(datetime.datetime.now())
                            + " Gen "+ str(new_tree.status_dictionary['gen']) + " :  updating gene-fitness dictionary")
     ## save
     new_tree.write_state()
########################

class tree_generation:
        def __init__(self,name):
                path_dictionary = setup_paths()
                ligands_list = get_ligands()
                self.base_path_dictionary = path_dictionary
                self.name = name
                self.genes =  dict()
                self.gene_fitness_dictionary = dict()
                self.ligands_list = ligands_list
                self.status_dictionary = dict()
                self.gene_compound_dictionary = dict()

        def configure_gen(self,gen_num,npool,ncross,pmut,genmax,scoring_function="split",split_parameter = 15.0,distance_parameter = 1,DFT =True, RTA = False,mean_fitness =  0):
                self.current_path_dictionary = advance_paths(self.base_path_dictionary,gen_num)
                self.status_dictionary.update({'gen':gen_num})
                self.status_dictionary.update({'scoring_function': scoring_function})
                self.status_dictionary.update({'split_parameter': split_parameter})
                self.status_dictionary.update({'distance_parameter': distance_parameter})
                self.status_dictionary.update({'npool':npool,'genmax':genmax})
                self.status_dictionary.update({'ncross': ncross})
                self.status_dictionary.update({'pmut': pmut})
                self.status_dictionary.update({'ready_to_advance':RTA})
                self.status_dictionary.update({'mean_fitness': mean_fitness})
                self.status_dictionary.update({'DFT': DFT})

                #FITNESS DEBUGGING  print "______________ scoring function: " + scoring_function

        def populate_random(self):
                ## clear the pool
                self.gene_compound_dictionary = dict()
                self.genes = dict()

                ### fill the pool with random structures 
                for i in range(0,self.status_dictionary['npool']):
                        this_complex = octahedral_complex(self.ligands_list)
                        this_complex.random_gen()
                        this_gene = this_complex.name
                        self.genes[i] = this_gene
                        self.gene_compound_dictionary[i] = this_complex

        def write_state(self):
                ## first write genes to path
                state_path = self.current_path_dictionary["state_path"] +"current_genes.csv"
                if not os.path.isfile(state_path):
                        open(state_path,'a').close()
                else:   ## backup state data
                        shutil.copyfile(state_path,self.current_path_dictionary["state_path"] +"current_genes.csv.bcp")
                emsg = write_dictionary(self.genes,state_path)
                ## second write live info to base directory
                state_path = self.base_path_dictionary["state_path"] +"/current_status.csv"
                if not os.path.isfile(state_path):
                        open(state_path,'a').close()
                emsg = write_dictionary(self.status_dictionary,state_path)
                if emsg:
                        print(emsg)

                ## third,  write gene-fitness info to path
                state_path = self.current_path_dictionary["state_path"] +"/gene_fitness.csv"
                if not os.path.isfile(state_path):
                        open(state_path,'a').close()
                emsg = write_dictionary(self.gene_fitness_dictionary,state_path)
                if emsg:
                        print(emsg)



        def read_state(self):
                ## first read live info from base directory
                state_path = self.base_path_dictionary["state_path"] +"/current_status.csv"
                emsg,read_dict = read_dictionary(state_path)
                if emsg:
                        print(emsg)
                self.configure_gen(gen_num = int(read_dict["gen"]),
                                   npool = int(read_dict["npool"]),
                                   ncross = int(read_dict["ncross"]),
                                   pmut = float(read_dict["pmut"]),
                                   genmax = int(read_dict["genmax"]),
                                   scoring_function =read_dict["scoring_function"],
                                   split_parameter = float(read_dict["split_parameter"]),
                                   distance_parameter = float(read_dict["distance_parameter"]),
                                   RTA = bool((read_dict["ready_to_advance"] == 'True')),
                                   mean_fitness = float(read_dict["mean_fitness"]),
                                   DFT =  bool((read_dict["DFT"] == 'True')))
                ## next read  genes from path
                state_path = self.current_path_dictionary["state_path"] +"current_genes.csv"
                emsg,gene_dict = read_dictionary(state_path)
                if emsg:
                        print(emsg)
                for keys in gene_dict.keys():
                    self.genes[int(keys)] = gene_dict[keys]
                for keys in self.genes.keys():
                        genes = self.genes[keys]
                        this_complex = octahedral_complex(self.ligands_list)
                        this_complex.encode(genes)
                        self.gene_compound_dictionary[keys] = this_complex

                ## third,  read gene-fitness info to path
                state_path = self.current_path_dictionary["state_path"] +"/gene_fitness.csv"
                emsg,fit_dict = read_dictionary(state_path)
                if emsg:
                        print(emsg)
                self.gene_fitness_dictionary = fit_dict

        def check_results(self):
                ## load gene fitness dict
                fitkeys  = self.gene_fitness_dictionary.keys()

                ## if doing a DFT run, we need to check the filestytem for updates
                if self.status_dictionary["DFT"]:
                        final_results = check_all_current_convergence()
                        for genes in final_results.keys():
                                if genes in fitkeys:
                                        print('gene ' + str(genes) + ' already in dict, no action')
                                else:
                                        this_split_energy = float(final_results[genes].split)
                                        if self.status_dictionary['scoring_function'] == "split+dist":
                                                print('error, cannot using aplit+dist fitness with ANN only. Switching to split only.')
                                                logger(self.base_path_dictionary['state_path'],str(datetime.datetime.now()) + ": Gen " +
                                                       str(self.status_dictionary['gen'] ) +
                                                      ' error, cannot using aplit+dist fitness with ANN only. Switching to split only')
                                        fitness =  find_split_fitness(this_split_energy,self.status_dictionary['split_parameter'])
                                        logger(self.base_path_dictionary['state_path'],str(datetime.datetime.now()) + ": Gen " +
                                               str(self.status_dictionary['gen'] ) +
                                               ' setting fitness to ' + str(fitness) + ' for new genes ' + str(genes))
                                        self.gene_fitness_dictionary.update({genes:fitness})
        def assess_fitness(self):
            ## loop all over genes in the pool and the selected set
            fitkeys  = self.gene_fitness_dictionary.keys()
            fitness_values  = dict()
            print('is code ready to advance?: '+str(self.status_dictionary["ready_to_advance"]))
            logger(self.base_path_dictionary['state_path'],str(datetime.datetime.now()) + ":  Gen " 
                       + str(self.status_dictionary['gen'])
                       + " is code ready to advance? " +str(self.status_dictionary["ready_to_advance"]))
            self.ready_to_advance = False
            self.outstanding_jobs = dict()
            for genekeys in self.genes.keys():
#                print('genekey is ' + str(genekeys))
                genes = self.genes[genekeys]
                ## see if this gene is in the fitness dictionary
                if genes in fitkeys:
                    fitness_values[genes] = self.gene_fitness_dictionary[genes]
                else:
                    self.outstanding_jobs.update({genekeys:self.gene_compound_dictionary[genekeys]})
            logger(self.base_path_dictionary['state_path'],str(datetime.datetime.now()) + ":  Gen " 
                       + str(self.status_dictionary['gen'])
                       + " with " + str(len(self.outstanding_jobs.keys())) + " calculations to be completed")
            print('length of outstanding jobskeys',len(self.outstanding_jobs.keys()))
            if (len(self.outstanding_jobs.keys()) ==0):
                logger(self.base_path_dictionary['state_path'],str(datetime.datetime.now()) 
                               + ": Gen " + str(self.status_dictionary['gen'])
                               + " all jobs completed, ranking ")
                self.status_dictionary["ready_to_advance"] = True
            else:
                self.job_dispatcher()
            ## if we are using the ANN only, populate the gene-fitnes dictionary
            if self.status_dictionary["DFT"] == False:
                    self.ANN_fitness()

        def random_fitness(self):
                ## test function for validating GA = white noise fitness 
                for keys in self.genes:
                        gene = self.genes[keys]
                        random_fitness = random.uniform(0,1)
                        logger(self.base_path_dictionary['state_path'],str(datetime.datetime.now()) 
                               + ": Gen " + str(self.status_dictionary['gen'])
                               + " assign random fitness  " + str(random_fitness) + ' to  gene ' + str(gene))
                        self.gene_fitness_dictionary.update({gene:random_fitness})
        def ANN_fitness(self):
                msg, ANN_dict = read_dictionary(self.current_path_dictionary["ANN_output"] +'ANN_results.csv')
                
                for keys in ANN_dict.keys():
                        gene,gen,slot,metal,ox,eq,ax1,ax2,spin,spin_cat,basename = translate_job_name(keys)
                        this_split_energy = float(ANN_dict[keys].split(',')[0])
                        this_ann_dist = float(ANN_dict[keys].split(',')[1].strip('\n'))

                        if self.status_dictionary['scoring_function'] == "split+dist":
                            fitness =  find_split_dist_fitness(this_split_energy,self.status_dictionary['split_parameter'],this_ann_dist,self.status_dictionary['distance_parameter'])
                        else:
                            fitness =  find_split_fitness(this_split_energy,self.status_dictionary['split_parameter'])

                        #FITNESS DEBUGGING
                        ##print gene + " fitness: " + str(fitness)

                        logger(self.base_path_dictionary['state_path'],str(datetime.datetime.now()) 
                               + ": Gen " + str(self.status_dictionary['gen'])
                               + " fitness from ANN  " + str(fitness) + ' assigned to  gene ' + str(gene))
                        self.gene_fitness_dictionary.update({gene:fitness})
        def job_dispatcher(self):
                jobpaths = list()
                emsg,ANN_results_dict = read_dictionary(self.current_path_dictionary["ANN_output"] +'/ANN_results.csv')
                current_outstanding = get_outstanding_jobs()
                converged_jobs = find_converged_job_dictionary()
                for keys in self.outstanding_jobs.keys():

                        jobs = self.outstanding_jobs[keys]
                        spins_dict = spin_dictionary()
                        metal = jobs.metals_list[jobs.core]
                        spin_list = spins_dict[metal][jobs.ox]
                        for spins in spin_list:
                                job_prefix = "gen_" + str(self.status_dictionary["gen"]) + "_slot_" + str(keys) + "_"
                                ## generate HS/LS
                               ## convert the gene into a job file and geometery
                                jobpath,mol_name,ANN_split,ANN_distance = jobs.generate_geometery(prefix = job_prefix, spin = spins,path_dictionary = self.current_path_dictionary,
                                                                      rundirpath = get_run_dir(), molsimpath = self.base_path_dictionary["molsimp_path"])
                                if (jobpath not in current_outstanding) and (jobpath not in converged_jobs.keys()):
                                        ## save result
                                        ANN_results_dict.update({mol_name:",".join([str(ANN_split),str(ANN_distance)])})
                                        jobpaths.append(jobpath)
                                        logger(self.base_path_dictionary['state_path'],str(datetime.datetime.now()) + ":  Gen " 
                                        + str(self.status_dictionary['gen'])
                                        + " missing information for gene number  " + str(keys) + ' with  name ' + str(jobs.name) )
                write_dictionary(ANN_results_dict,self.current_path_dictionary["ANN_output"] +'ANN_results.csv')
                set_outstanding_jobs(current_outstanding+jobpaths)


        def select_best_genes(self):
                ## first write genes to path
                summary_path = self.current_path_dictionary["state_path"] +"all_genes.csv"
                outcome_list = list()
                npool  =  self.status_dictionary["npool"]
                mean_fitness = 0
                for keys in self.genes.keys():
                    outcome_list.append((keys,self.genes[keys],float(self.gene_fitness_dictionary[self.genes[keys]])))
                    logger(self.base_path_dictionary['state_path'],str(datetime.datetime.now()) +
                               ": Gen " + str(self.status_dictionary['gen']) +  '  gene is   ' + str(keys) 
                             + " fitness is = " +  str(float(self.gene_fitness_dictionary[self.genes[keys]])))
 
                outcome_list.sort(key=lambda tup: tup[2], reverse = True)

                full_size = len(outcome_list)

                if not os.path.isfile(summary_path):
                       open(summary_path,'a').close()
                emsg = write_summary_list(outcome_list,summary_path)
                self.genes = dict()
                self.gene_compound_dictionary = dict()
                for i in range(0,npool):
                        self.genes[i] = outcome_list[i][1]
                        this_complex = octahedral_complex(self.ligands_list)
                        this_complex.encode(self.genes[i])
                        self.gene_compound_dictionary[i] = this_complex
                        mean_fitness += float(outcome_list[i][2])
                mean_fitness = mean_fitness/npool # average fitness
                self.status_dictionary.update({'mean_fitness':mean_fitness})
                logger(self.base_path_dictionary['state_path'],str(datetime.datetime.now()) +
                       ": Gen " + str(self.status_dictionary['gen']) 
                     + " complete, mean_fitness = " +  str(mean_fitness))
        def advance_generation(self):
                ## advance counter
                self.status_dictionary['gen'] +=1
                logger(self.base_path_dictionary['state_path'],str(datetime.datetime.now()) +
                       ": Gen " + str(self.status_dictionary['gen']-1) 
                     + " advancing to Gen " +  str(self.status_dictionary['gen']))
                self.status_dictionary['ready_to_advance'] = False
                self.current_path_dictionary = advance_paths(self.base_path_dictionary,self.status_dictionary['gen'])

                npool  =  self.status_dictionary["npool"]
                ncross =  self.status_dictionary["ncross"]
                pmut   =  self.status_dictionary["pmut"]

                ## generation selected set
                selected_genes = dict()
                selected_compound_dictionary = dict()
                number_selected = 0
                ## populate selected pool
                while number_selected < npool:
                        this_int = random.randint(0,npool -1)
                        this_barrier = random.uniform(0,1)
                        this_gene = self.genes[this_int]
                        if self.gene_fitness_dictionary[this_gene] > this_barrier:
                                selected_genes[number_selected + npool] = this_gene
                                number_selected += 1
                ## populate compound list
                for keys in selected_genes.keys():
                        genes = selected_genes[keys]
                        this_complex = octahedral_complex(self.ligands_list)
                        this_complex.encode(genes)
                        selected_compound_dictionary[keys] = this_complex
                ## now perfrom ncross exchanges
                number_of_crosses = 0
                while number_of_crosses < ncross:
                        these_partners = random.sample(range(npool,(2*npool - 1)),2)
                        keep_axial = selected_compound_dictionary[these_partners[0]]
                        keep_equitorial = selected_compound_dictionary[these_partners[1]]
                        old_genes = [selected_genes[key] for key in these_partners]
                        new_complex_1 = keep_axial.exchange_ligands(keep_equitorial,True)
                        new_complex_1 = new_complex_1.exchange_metal(keep_equitorial)
                        new_complex_1 = new_complex_1.exchange_ox(keep_equitorial)
                        new_complex_2 = keep_equitorial.exchange_ligands(keep_axial,True)
                        new_complex_2 = new_complex_2.exchange_metal(keep_axial)
                        new_complex_2 = new_complex_2.exchange_ox(keep_axial)
                        new_gene_1 = new_complex_1.name
                        new_gene_2 = new_complex_2.name
                        selected_genes[these_partners[0]] = new_gene_1
                        selected_compound_dictionary[these_partners[0]] = new_complex_1
                        selected_genes[these_partners[1]] = new_gene_2
                        selected_compound_dictionary[these_partners[1]] = new_complex_2
                        new_genes = [selected_genes[key] for key in these_partners]

                        number_of_crosses +=1
                        logger(self.base_path_dictionary['state_path'],str(datetime.datetime.now()) +
                               ":  Gen " + str(self.status_dictionary['gen'])
                               + " crossing " + str(these_partners) + " " +
                              str(old_genes) + " -> " + str(new_genes)  )

                ## mutate
                for keys in selected_genes.keys():
                        does_mutate = random.uniform(0,1)
                        if does_mutate < pmut:
                                old_gene = selected_genes[keys]
                                mutant = selected_compound_dictionary[keys].mutate()
                                selected_compound_dictionary[keys] = mutant
                                selected_genes[keys] = selected_compound_dictionary[keys].name
                                logger(self.base_path_dictionary['state_path'],str(datetime.datetime.now()) +
                                       ":  Gen " + str(self.status_dictionary['gen'])
                                       + " mutating " + str(keys) + ": "  + old_gene +  " -> " + mutant.name) 

                ## merge the lists 
                self.genes.update(selected_genes)
                self.gene_compound_dictionary.update(selected_compound_dictionary)








