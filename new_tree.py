from ga_init import *
maxgen = 200
npool = 20
pmut  = 0.25
ncross = 5
scoring_function = "split"
#scoring_function = "split+dist"
split_parameter = 10
distance_parameter = 0.5
t1   = initialize_GA_calc(npool,ncross,pmut,maxgen,scoring_function,split_parameter,distance_parameter)
