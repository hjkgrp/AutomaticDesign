import glob
import datetime
import math
import numpy
import subprocess
import argparse
import os
import random
import shutil
from molSimplifyAD.ga_complex import *
from molSimplifyAD.ga_main import *
from molSimplifyAD.ga_check_jobs import *
from molSimplifyAD.ga_io_control import *
def initialize_GA_calc(path = ''):
        
        ## load in run info
        GA_run = GA_run_defintion()
        

        GA_run.deserialize(path + '.madconfig')
        
        path_dictionary = setup_paths()
        
        new_tree = GA_generation('current_GA')
        new_tree.configure_gen(0,**GA_run.config)

        new_tree.populate_random()
        new_tree.write_state()
        logger(new_tree.base_path_dictionary['state_path'],str(datetime.datetime.now())
               + ": <new tree>  Gen : " + str(new_tree.status_dictionary['gen']) + ' commencing')
        return new_tree

def wake_up_routine():
        ## set up environment:
        path_dictionary = setup_paths()
        ## initialize class
        new_tree = GA_generation('current_gen')
        ## read in info
        new_tree.read_state()
        current_gen = new_tree.status_dictionary["gen"]
        maxgen = new_tree.status_dictionary["maxgen"]

        logger(new_tree.base_path_dictionary['state_path'],str(datetime.datetime.now())
               + ": <resuming>  Gen : " + str(new_tree.status_dictionary['gen']))
        ## assess current fitness
        new_tree.check_results()
        new_tree.assess_fitness()
        print(new_tree.status_dictionary["ready_to_advance"])
        if current_gen <= maxgen and GA_run.config["type"]=="split" :
               ## check if there is still
                ## work to be done
                if (new_tree.status_dictionary["ready_to_advance"]):
                        print('ready to advance')
                        ## if all jobs are complete:
                        ## choose genes for next gen
                        new_tree.select_best_genes()
                        new_tree.write_state()
                        ## advance to next gen:
                        new_tree.advance_generation()
                        new_tree.assess_fitness()
                        new_tree.write_state()
                else:
                        logger(new_tree.base_path_dictionary['state_path'],str(datetime.datetime.now())
                               + ":  Gen  " + str(new_tree.status_dictionary['gen']) + ' has calcs oustanding,waiting')
                        new_tree.write_state()
        else:
            logger(new_tree.base_path_dictionary['state_path'],str(datetime.datetime.now())
                               + ":  Gen  " + str(new_tree.status_dictionary['gen']) + ' has reached maxgen')

        return(new_tree)
