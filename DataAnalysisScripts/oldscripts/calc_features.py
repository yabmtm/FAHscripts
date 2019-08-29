
import os
import sys
import mdtraj as md
import numpy as np
from analysis_featurizer import *

base_dir = '/home/server/server2/data/SVR166219/'


def calc_distances(proj_num):

    # creates a distances.npy file for the project if there isn't one and loads it in
    create_npy_file('distances.npy', proj_num)
    npy_file = os.path.join(base_dir, 'PROJ' + proj_num, 'features', 'distances.npy')
    np_arr = np.load(npy_file)

    # loops through all the runs, clones, gens of the project calculating its distances
    for run_num, run in enumerate(np_arr):
	for clone_num, clone in enumerate(run):
	    for gen_num, gen in enumerate(clone):
		# if this index of the array has no data in it and if there is data for it, calculate its distances
		if gen is None and os.path.isdir(os.path.join(base_dir, 'PROJ' + str(proj_num), 'RUN' + str(run_num), 'CLONE' + str(clone_num), 'results' + str(gen_num))):
		    traj = Featurizer(proj_num, run_num, clone_num, gen_num)
		    if traj.xtc_file is None:
			print('Could not find xtc_file for run %d, clone %d, results %d' % (str(run_num), str(clone_num), str(gen_num)))
			continue
	   	    if traj.gro_file is None:
			print('Could not find gro_file for run %d, clone %d, results %d' % (str(run_num), str(clone_num), str(gen_num)))
			continue
		    print('Calculating distances for run %d, clone %d, results %d' % (run_num, clone_num, gen_num))
		    np_arr[run_num][clone_num][gen_num] = traj.calculate_distances(False)
		    # saves distances for every 10 gens that are calculated
		    if gen_num % 10 == 0 and gen_num is not 0:
			print('Saving distances to distances.npy...')
			np.save(npy_file, np_arr)
			print('Distances saved')
		

						

try:	
    calc_distances(sys.argv[1])
except Exception as e:
    print('Invalid input: ' + repr(e))






























