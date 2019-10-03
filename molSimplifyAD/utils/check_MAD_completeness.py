#Call this function from the directory aboce the mad directory. You can hand it the name of multiple mad runs at once to check if they are all complete
#Use like:
# python check_mad_completeness.py <name of mad directory 1> <name of mad directory 2> etc.

import os
import pandas as pd
from sys import argv
def load(path):
	fil = open(path,'r')
	lines = fil.readlines()
	fil.close()
	new_lines = []
	for i in lines:
		if i[-1] =='\n':
			new_lines.append(i[:-1])
		else:
			new_lines.append(i)
	return new_lines

def check_finished(debug=False):
	#Checks if a MAD directory located at the current working directory is finished running or not
	if not os.path.exists(os.path.join('jobs','outstanding_job_list.txt')):
		raise ValueError('Outstanding job list does not exist!')
	if not os.path.exists(os.path.join('jobs','live_jobs.csv')):
		raise ValueError('Live jobs csv does not exist!')
	if not os.path.exists(os.path.join('jobs','converged_job_dictionary.csv')):
		raise ValueError('Converged job list does not exist')
	if not os.path.exists('.madconfig'):
		raise ValueError('No madconfig!')
	if not os.path.exists(os.path.join('jobs','submitted_jobs.csv')):
		raise ValueError('No submitted job list!')

	#Find the max_resub limit
	madconfig = load('.madconfig')[0]
	for i in madconfig.split(','):
		if 'max_resubmit' in i:
			max_resub = int(i[-1])
	#Find a dictionary of the number of times each job has been submitted
	submitted = load(os.path.join('jobs','submitted_jobs.csv'))
	submitted = [i.split(',') for i in submitted]
	submitted_dict = dict()
	for i in submitted:
		submitted_dict[i[0]]=int(i[1])

	#identify outstanding jobs which are not status 7
	outstanding_jobs = load(os.path.join('jobs','outstanding_job_list.txt'))
	live_jobs = load(os.path.join('jobs','live_jobs.csv'))
	converged_jobs = load(os.path.join('jobs','converged_job_dictionary.csv'))
	submission_limit_jobs = [i.rsplit(',',1)[0] for i in converged_jobs if int(i.rsplit(',',1)[1])==7]

	to_submit = list(set(outstanding_jobs)-set(submission_limit_jobs))

	#Eliminate jobs which are above the sumbission limit
	passed_limit = []
	for i in to_submit:
		try:
			if submitted_dict[i] >= (max_resub+1):
				passed_limit.append(i) 
		except:
			#if the job hasn't been submitted a single time yet, then it won't be in the submit dict
			return 'In Progress'

	to_submit = list(set(to_submit)-set(passed_limit))

	if debug:
		for i in to_submit:
			for counter,ii in enumerate([iii.rsplit(',',1)[0] for iii in converged_jobs]):
				if i == ii:
					print('job: '+i+' is reportng status: '+converged_jobs[counter].rsplit(',',1)[1])
	if len(to_submit) == 0:
		return 'Finished'
	else:
		return 'In Progress'

if __name__ == '__main__':
	args = argv[1:]
	home = os.getcwd()
	for i in args:
		os.chdir(i)
		finished = check_finished()
		os.chdir(home)
		print('Status of '+i+': '+finished)
		
