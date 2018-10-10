from molSimplifyAD.ga_tools import *
import shutil

rundir = str(isKeyword('rundir'))
pathToConvergedJobs = rundir+'jobs/converged_job_dictionary.csv'
pathToOutstandingJobs = rundir+'jobs/outstanding_job_list.txt'
pathToSubmittedJobs = rundir+'jobs/submitted_jobs.csv'
shutil.copy(pathToConvergedJobs,pathToConvergedJobs + '.backup')
shutil.copy(pathToOutstandingJobs,pathToOutstandingJobs + '.backup')
shutil.copy(pathToSubmittedJobs,pathToSubmittedJobs + '.backup')
with open(pathToConvergedJobs+'.backup', 'rb') as infile:
	outfile = open(pathToConvergedJobs,'wb')
	for line in infile:
		if rundir in line:
			continue
		else:
			index = line.find("jobs/")
			if index == -1:
				print('could not change path in line:' +line)
			else:
				line = line.replace(line[:index],rundir)
				outfile.write(line)
	outfile.close()
with open(pathToOutstandingJobs+'.backup', 'rb') as infile:
	outfile = open(pathToOutstandingJobs,'wb')
	for line in infile:
		if rundir in line:
			continue
		else:
			index = line.find("jobs/")
			if index == -1:
				print('could not change path in line:' +line)
			else:
				line = line.replace(line[:index],rundir)
				outfile.write(line)
	outfile.close()
with open(pathToSubmittedJobs+'.backup', 'rb') as infile:
	outfile = open(pathToSubmittedJobs,'wb')
	for line in infile:
		if rundir in line:
			continue
		else:
			index = line.find("jobs/")
			if index == -1:
				print('could not change path in line:' +line)
			else:
				line = line.replace(line[:index],rundir)
				outfile.write(line)
	outfile.close()
	
