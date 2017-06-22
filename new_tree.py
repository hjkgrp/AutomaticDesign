from ga_init import *
maxgen = 200
npool = 10
pmut  = 1
ncross = 3
scoring_function = "split"
scoring_function = "split+dist"
t1   = initialize_GA_calc(npool,ncross,pmut,maxgen,scoring_function)
