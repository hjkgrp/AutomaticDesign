from molSimplifyAD.ga_init import *
from molSimplifyAD.ga_io_control import *
from molSimplifyAD.ga_get_general import *
from molSimplifyAD.get_distances import *
import os, datetime,time
def resume_run(args):
  print(args)
  
  ## change to  run directory
  with switch_to_rundir(args.resume):
    print(os.getcwd())
    if args.reps:
      reps = args.reps
    else:
      reps = 1
    its = 0
    while its < reps:
      wake_up_routine()
      format_freqeuncies()
      format_distances()
      its +=1
      if args.sleep:
        print('sleeping for '+ str(args.sleep))
        time.sleep(args.sleep) 
        
      

  

