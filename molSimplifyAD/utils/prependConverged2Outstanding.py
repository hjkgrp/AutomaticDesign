import os, sys, shutil


slotsToResub = [72,73,74] # here we can put certain slots to fix

pathToConvergedJobs = "converged_job_dictionary.csv"
newPathToConvergedJobs = "new_converged_job_dictionary.csv"
PathToOutstandingJobs = "outstanding_job_list.txt"
newPathToOutstandingJobs = "new_outstanding_job_list.txt"

# holder of lines to add to outstanding
OutstandingLinesToAdd = []

# backup jobs:
shutil.copy(pathToConvergedJobs,pathToConvergedJobs + '.backup')
shutil.copy(PathToOutstandingJobs,PathToOutstandingJobs + '.backup')

with open(pathToConvergedJobs,'r') as f: # this is the file with the converged jobs
    with open(newPathToConvergedJobs,'w') as newf: # this is the new file kept converged jobs
        for line in f.readlines():
            thisLineStatus = int(line.split(",")[1].strip()) # get status of this line
            thisSlot = int(os.path.basename(line.split(",")[0]).split("_")[3]) # this is the slot number (as int)
            if thisSlot in slotsToResub:
                # check that it is status 6:
                if thisLineStatus == 6:
                    print(('removing '+line.strip()+ ' from converged'))
                    OutstandingLinesToAdd.append(line.split(",")[0]) # store these for later
                else:
                    print(('error '+line+ ' would be removed but has stats ' + str(thisLineStatus)))
                    sys.exit() # check it is six
            else: # this is not the lines you were looking for
                newf.write(line)

# now we **could** copy and overwrite

#shutil.copy(newPathToConvergedJobs,pathToConvergedJobs)

# now they are removed from outstanding, need to be added to the outstanding jobs
with open(newPathToOutstandingJobs, 'w') as newf: # start a new outstanding jobs list
    for line in OutstandingLinesToAdd: # add our new jobs to the top
        newf.write(line + '\n')
    # add the other outstanding jobs
    with open(PathToOutstandingJobs, 'r') as f:
        for line in f.readlines():
            newf.write(line)

# now we **could** copy and overwrite

#shutil.copy(newPathToOutstandingJobs,PathToOutstandingJobs)



