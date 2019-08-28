import os
import mdtraj as md
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from itertools import combinations
import sys
from bs4 import BeautifulSoup
from analysis_featurizer import *

# base directories
base_proj_dir = '/home/server/server2/data/SVR166219/'
base_gro_dir = '/home/server/server2/projects/Gromacs/'


def RMSD_analysis(proj_num):
    # loads in the distances.npy file for this project
    dist_file = os.path.join(base_proj_dir, 'PROJ' + proj_num, 'features', 'distances.npy')
    runs, clones, gens = 0, 0, 0
    try:	
        dist = np.load(dist_file)
        runs, clones, gens = dist.shape[0], dist.shape[1], dist.shape[2]
	# if there is already an RMSD_data file it will load it in, otherwise it will create one and load it
        RMSD_file = os.path.join(base_proj_dir, 'PROJ' + proj_num, 'features', 'RMSD_data.npy')
        RMSD_data = [[[]]]
        if os.path.isfile(RMSD_file):
	    RMSD_data = np.load(RMSD_file)
        else:
	    create_npy_file('RMSD_data.npy', proj_num)
	    RMSD_data = np.load(RMSD_file)
	# loops through every run, clone, gen of the projects calculating RMSD values for all data found
        for i in range(runs):
            for j in range(clones):
                for k in range(gens):
		    # if that index in RMSD_data has no data and there is data for it in the distances.npy file, calculate RMSD values
		    if RMSD_data[i][j][k] is None and dist[i][j][k] is not None:
                        try:
                            RMSD = Analysis(proj_num, i, j, k)
                            RMSD.calculate_RMSD()
                        except Exception as e:
                            print(repr(e) + '\nError with calculating values for run %d, clone %d, gen %d. Skipping to next clone' % (i, j, k))
                            break

    except Exception as e:
	print(repr(e))
	exit()
    print('All RMSD values have been calculated')


def RMSD_plot(proj_num):
    # uses analysis_featurizer.py analysis class to plot the RMSD values for this project
    analysis_object = Analysis(proj_num, 0, 0, 0)
    analysis_object.plot_RMSD()




if str(sys.argv[1]) == 'plot':
    RMSD_plot(sys.argv[2])
elif str(sys.argv[1]) == 'update':
    RMSD_analysis(sys.argv[2])
else:
    print 'input invalid'











