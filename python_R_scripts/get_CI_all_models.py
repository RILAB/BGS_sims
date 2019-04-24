import sys
import os
import numpy as np
import random
import math
import commands
from scipy import stats

#This script will take in X number of bootstrapped simulation iterations (the bootstrap distribution)
#containing a particular statistic of interest (pi, xi, or tajD [Tajima's D])
#for a given demographic and selection (neutral or bgs) model and in turn
#generate 95% confidence intervals and standard errors of the mean

#arguments
directory = sys.argv[1] #directory where bootstrap data files reside in
outdirectory = sys.argv[2] #directory to write output files to
model = sys.argv[3] #population demographic model (e.g., 'model1', 'model2', 'model3')
type = sys.argv[4] #statistic of interest, either 'pi', 'xi', or 'tajD'

os.chdir(directory)

#get list of bootstrap data files (list of all bootstrap iterations)
bgsfilelist = commands.getoutput("ls bgs_mean_" + model + "_" + type + "_BSitr*.txt").strip().split()
neutfilelist = commands.getoutput("ls neut_mean_" + model + "_" + type + "_BSitr*.txt").strip().split()

bootstrap_itrs = len(neutfilelist)

myfile = "neut_mean_" + model + "_" + type + "_BSitr1.txt"

genlist = commands.getoutput("awk '{print $1}' " + myfile + " | sort -n | uniq").strip().split()

#make dictionaries to store list of sampled means for each generation from all bootstrap iterations
mean_neut_gendic = {}
mean_bgs_gendic = {}

for gen in genlist:
	mean_neut_gendic[gen] = []
	mean_bgs_gendic[gen] = []

for neutfile in neutfilelist:
	infile = open(neutfile, 'r')
	for line in infile:
		line = line.strip().split()
		mean_neut_gendic[line[0]].append([float(x) for x in line[1:]])
	infile.close()

for bgsfile in bgsfilelist:
	infile = open(bgsfile, 'r')
	for line in infile:
		line = line.strip().split()
		mean_bgs_gendic[line[0]].append([float(x) for x in line[1:]])
	infile.close()

'''structure of mean_neut_gendic and mean_bgs_gendic:
Both are dictionaries where the key is a generation (i.e., '-0.01', '0.0', '0.01' ... '1.0')
and the value is a list of lists (bootstrap_itrs sublists in each list), where each sublist is list of
mean values corresponding to the ith window of the simulated region for a particular set of simulations.
The final value in each sublist of mean_bgs_gendic is the central region experiencing selection.'''

windowlength = len(mean_neut_gendic[genlist[0]][0])

ratio_dic = {}

for gen in genlist:
	ratio_dic[gen] = []

for gen in genlist:
	for i in xrange(0, int(bootstrap_itrs)):
		newlist = []
		for j in xrange(0, windowlength):
			if type == 'tajD':
				newlist.append(mean_bgs_gendic[gen][i][j]-mean_neut_gendic[gen][i][j])
			else:
				newlist.append(mean_bgs_gendic[gen][i][j]/mean_neut_gendic[gen][i][j])
		ratio_dic[gen].append(newlist)

sems_neut_dic = {}
sems_neut_dic['outfile'] = open(outdirectory + '/neut_sems_' + model + '_' + type + '.txt', 'w')
sems_bgs_dic = {}
sems_bgs_dic['outfile'] = open(outdirectory + '/bgs_sems_' + model + '_' + type + '.txt', 'w')
sems_ratio_dic = {}
sems_ratio_dic['outfile'] = open(outdirectory + '/ratio_sems_' + model + '_' + type + '.txt', 'w')

neut_95CI_lower_dic = {}
neut_95CI_lower_dic['outfile'] = open(outdirectory + '/neut_95CI_lower_' + model + '_' + type + '.txt', 'w')
neut_95CI_upper_dic = {}
neut_95CI_upper_dic['outfile'] = open(outdirectory + '/neut_95CI_upper_' + model + '_' + type + '.txt', 'w')

bgs_95CI_lower_dic = {}
bgs_95CI_lower_dic['outfile'] = open(outdirectory + '/bgs_95CI_lower_' + model + '_' + type + '.txt', 'w')
bgs_95CI_upper_dic = {}
bgs_95CI_upper_dic['outfile'] = open(outdirectory + '/bgs_95CI_upper_' + model + '_' + type + '.txt', 'w')

ratio_95CI_lower_dic = {}
ratio_95CI_lower_dic['outfile'] = open(outdirectory + '/ratio_95CI_lower_' + model + '_' + type + '.txt', 'w')
ratio_95CI_upper_dic = {}
ratio_95CI_upper_dic['outfile'] = open(outdirectory + '/ratio_95CI_upper_' + model + '_' + type + '.txt', 'w')

lower = int(math.floor(bootstrap_itrs*0.025)-1)
upper = int(math.floor(bootstrap_itrs*0.975)-1)

for gen in genlist:
	mean_neut_gendic[gen] = np.asarray(mean_neut_gendic[gen])
	mean_bgs_gendic[gen] = np.asarray(mean_bgs_gendic[gen])
	ratio_dic[gen] = np.asarray(ratio_dic[gen])

	sems_neut_dic[gen] = [str(x) for x in list(np.std(mean_neut_gendic[gen], dtype=np.float64, axis=0))]
	sems_bgs_dic[gen] = [str(x) for x in list(np.std(mean_bgs_gendic[gen], dtype=np.float64, axis=0))]
	sems_ratio_dic[gen] = [str(x) for x in list(np.std(ratio_dic[gen], dtype=np.float64, axis=0))]

	neut_95CI_lower_dic[gen] = [str(x) for x in list(np.sort(mean_neut_gendic[gen], axis=0)[lower,:])]
	neut_95CI_upper_dic[gen] = [str(x) for x in list(np.sort(mean_neut_gendic[gen], axis=0)[upper,:])]

	bgs_95CI_lower_dic[gen] = [str(x) for x in list(np.sort(mean_bgs_gendic[gen], axis=0)[lower,:])]
	bgs_95CI_upper_dic[gen] = [str(x) for x in list(np.sort(mean_bgs_gendic[gen], axis=0)[upper,:])]

	ratio_95CI_lower_dic[gen] = [str(x) for x in list(np.sort(ratio_dic[gen], axis=0)[lower,:])]
	ratio_95CI_upper_dic[gen] = [str(x) for x in list(np.sort(ratio_dic[gen], axis=0)[upper,:])]

for i in [sems_neut_dic, sems_bgs_dic, sems_ratio_dic, neut_95CI_lower_dic, neut_95CI_upper_dic, bgs_95CI_lower_dic, bgs_95CI_upper_dic, ratio_95CI_lower_dic, ratio_95CI_upper_dic]:
	for gen in genlist:
		i['outfile'].write(gen + ' ' + ' '.join(i[gen]) + '\n')
	i['outfile'].close()
