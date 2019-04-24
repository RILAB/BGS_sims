import sys
import numpy as np
import random
import math
import commands
from scipy import stats

#This script will take in a data file of X number of simulation replicates containing a particular statistic of interest (pi, xi, or tajD [Tajima's D]) for a given demographic and selection (neutral or bgs) model and either:
#1) bootstrap (sample with replacement) those replicates and generate a mean of the statistic from the newly sampled replicates
# - OR - 
#2) not resample the simulation replicates, but calculate the mean of the statistic directly from all replicates

#arguments
directory = sys.argv[1] #directory simulation data files reside in
outdirectory = sys.argv[2] #directory to write bootstrapped files to
model = sys.argv[3] #population demographic model (e.g., 'model1', 'model2', 'model3')
typeof = sys.argv[4] #statistic of interest, either 'pi', 'xi', or 'tajD'
bootstrap_itrs = sys.argv[5] #index of bootstrap iteration performed. If == '0', then this script will only return the mean of all the simulation replicates

neutfile = open(directory + '/' + model + '_neutral_' + typeof + '.csv', 'r') #neutral simulations data file
bgsfile = open(directory + '/' + model + '_bgs_' + typeof + '.csv', 'r') #background selection (bgs) simulations data file
myfile = directory + '/' + model + '_neutral_' + typeof + '.csv' #variable representing a filename for running UNIX commands on (see below)

neutfile.readline()
header = bgsfile.readline()

#make dictionaries for storing data file records in
myneutdic = {}
mybgsdic = {}

#get total number of independent sims (simulation replicates) that were used to generate the data files
tot_sims = int(commands.getoutput("awk 'BEGIN {FS=\",\"} {print $NF}' " + myfile + " | sort -n | tail -1").strip().split()[0])+1

#each dictionary will have keys representing the simulation replicates of the data file, which ranges from 0 to tot_sims-1 (tot_sims being the total number of simulation replicates)
#fill dictionaries with keys, make corresponding value a list
for i in xrange(0, tot_sims):
	myneutdic[str(i)] = []
	mybgsdic[str(i)] = []

#read each line of the neutral sims file, assign a key, and store in the dictionary
for line in neutfile:
	line = line.strip().split(',')
	rep = line[-1] #sim replicate number (will be our key)
	gen = line[0] #first column of line (this column is also the generation number of the particular demographic model)
	windows = [gen] #start a list for all columns in the line
	windows.extend([float(x) for x in line[1:-1]]) #store all other columns of the line in the list (these other columns contain values about the statistic of interest)
	myneutdic[rep].append(windows) #add simulation replicate to dictionary

#read each line of the bgs sims file, assign a key, and store in the dictionary
for line in bgsfile:
	line = line.strip().split(',')
	rep = line[-1] #sim replicate number (will be our key)
	gen = line[0] #first column of line (this column is also the generation number of the particular demographic model)
	windows = [gen] #start a list for all columns in the line
	windows.extend([float(x) for x in line[1:-1]]) #store all other columns of the line in the list (these other columns contain values about the statistic of interest for every genomic window in the simulation)
	mybgsdic[rep].append(windows) #add simulation replicate to dictionary

#make a new dictionary for putting the means of our re-sampled (or directly sampled if bootstrap_itrs==0) replicates
#each key will be a specific generation number from the demographic model
mean_neut_gendic = {}
mean_bgs_gendic = {}

#get generation numbers of the particular demographic model
genlist = commands.getoutput("awk 'BEGIN {FS=\",\"} {print $1}' " + myfile + " | sort -n | uniq | grep -v 'gen'").strip().split()

#add keys to dictionary representing generation numbers, make corresponding value a list
for gen in genlist:
	mean_neut_gendic[gen] = []
	mean_bgs_gendic[gen] = []

#make empty dictionaries to put list of sampled (or directly sampled if bootstrap_itrs==0) replicates
rep_neut_dic = {} 
rep_bgs_dic = {}

#add keys to the dictionaries, make corresponding value a list
for gen in genlist:
	rep_neut_dic[gen] = []
	rep_bgs_dic[gen] = []

#if bootstrap_itrs!=0, sample with replacement the simulation replicates. We will calcluate a mean from these new replicate samples for the statistic of interest. This mean is one full bootstrapped iteration!
#if bootstrap_itrs==0, we will not sample with replacement the simulation replicates but rather just sample each replicate directly. From this we will calculate the mean for the statistic of interest. This is the observed mean from all simulation replicates!
for rep in xrange(0, tot_sims):
	if bootstrap_itrs == '0':
		r = str(rep) #index of replicate, not sampled with replacement (i.e., not bootstrapped)
	else:
		r = str(random.randint(0, tot_sims-1)) #index of random replicate (sampled with replacement for bootstrapping)
	for neutgen, bgsgen in zip(myneutdic[r], mybgsdic[r]): #go through generations of selected replicate for both neutral and bgs simulations
		rep_neut_dic[neutgen[0]].append(neutgen[1:]) #add replicate to dictionary of re-sampled replicates
		rep_bgs_dic[bgsgen[0]].append(bgsgen[1:]) #add replicate to dictionary of re-sampled replicates

#now generate the means of the statistic of interest from the sampled replicates
for gen in genlist:
	rep_neut_dic[gen] = np.asarray(rep_neut_dic[gen])
	rep_bgs_dic[gen] = np.asarray(rep_bgs_dic[gen])
	meanwindow_neut_list = [] #make list to put mean of tot_sims windows for each particular genomic window
	meanwindow_bgs_list = [] #make list to put mean of tot_sims windows for each particular genomic window
	for window in xrange(0, len(rep_neut_dic[gen][0])):
		mean_neut_window = sum(rep_neut_dic[gen][:,window])/(float(tot_sims)) #mean of genomic window across all sims
		mean_bgs_window = sum(rep_bgs_dic[gen][:,window])/(float(tot_sims)) #mean of genomic window across all sims
		meanwindow_neut_list.append(mean_neut_window)
		meanwindow_bgs_list.append(mean_bgs_window)
	mean_neut_gendic[gen].append(meanwindow_neut_list)
	mean_bgs_gendic[gen].append(meanwindow_bgs_list)

'''structure of mean_neut_gendic and mean_bgs_gendic:
Both are dictionaries where the key is a generation (i.e., '-0.01', '0.0', '0.01' ... '1.0')
and the value is a list of a single sublist, where the sublist is a list of mean values corresponding to every
ith window of the simulated region for a particular set of tot_sims simulations.
If bootstrap_itrs == '0', then mean_neut_gendic and mean_bgs_gendic is list containing the mean values for each
window for each generation for all unique tot_sims (i.e., not randomly sampled with replacement).'''

#if boostrap_itrs == '0', then calculate ratios of bgs to neutral (or difference if typeof == 'tajD') and write out files for means of bgs and neutral generations and their ratios
#if boostrap_itrs != '0', then write out files with the means of the bgs and neutral randomly sampled replicate generations created from the above sampling procedure
if bootstrap_itrs == '0':
	ratio_dic = {}
	
	for gen in genlist:
		ratio_dic[gen] = []
	
	for gen in genlist:
		newlist = []
		for j in xrange(0, len(mean_bgs_gendic[gen][0])):
			if typeof == 'tajD':
				newlist.append(mean_bgs_gendic[gen][0][j]-mean_neut_gendic[gen][0][j])
			else:
				newlist.append(mean_bgs_gendic[gen][0][j]/mean_neut_gendic[gen][0][j])
		ratio_dic[gen].append(newlist)
			
	mean_neut_gendic['outfile'] = open(outdirectory + '/neut_mean_' + model + '_' + typeof + '.txt', 'w')
	mean_bgs_gendic['outfile'] = open(outdirectory + '/bgs_mean_' + model + '_' + typeof + '.txt', 'w')
	ratio_dic['outfile'] = open(outdirectory + '/ratio_mean_' + model + '_' + typeof + '.txt', 'w')
	
	for gen in genlist:
		mean_neut_gendic[gen] = [str(x) for x in mean_neut_gendic[gen][0]]
		mean_bgs_gendic[gen] = [str(x) for x in mean_bgs_gendic[gen][0]]
		ratio_dic[gen] = [str(x) for x in ratio_dic[gen][0]]
	
	for i in [mean_neut_gendic, mean_bgs_gendic, ratio_dic]:
		for gen in genlist:
			i['outfile'].write(gen + ' ' + ' '.join(i[gen]) + '\n')
		i['outfile'].close()
else:
	outfile_neut = open(outdirectory + '/neut_mean_' + model + '_' + typeof + '_BSitr' + str(bootstrap_itrs) + '.txt', 'w')
	outfile_bgs = open(outdirectory + '/bgs_mean_' + model + '_' + typeof + '_BSitr' + str(bootstrap_itrs) + '.txt', 'w')
	
	for gen in genlist:
		mean_neut_gendic[gen] = [str(x) for x in mean_neut_gendic[gen][0]]
		mean_bgs_gendic[gen] = [str(x) for x in mean_bgs_gendic[gen][0]]
		outfile_neut.write(gen + ' ' + ' '.join(mean_neut_gendic[gen]) + '\n') 
		outfile_bgs.write(gen + ' ' + ' '.join(mean_bgs_gendic[gen]) + '\n')
		
	outfile_neut.close()
	outfile_bgs.close()
