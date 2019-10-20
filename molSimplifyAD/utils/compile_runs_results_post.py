#Hand in a list of mad directories to this program, and it will concatenate the runs results post from each one
#Writes all of the runs results into a single file called compiled_runs_results_post.csv, in the current working directory
#Use like:
# python compile_runs_results_post.py <name of mad directory 1> <name of mad directory 2> etc.


import os
import pandas as pd
from sys import argv


args = argv[1:]
home = os.getcwd()

results = []
for i in args:
    os.chdir(i)
    results.append(pd.read_csv('runs_results_post.csv'))
    os.chdir(home)
    
results = pd.concat(results,axis=0)
results.to_csv('compiled_runs_results_post.csv')
