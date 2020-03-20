#!/home/server/anaconda2/bin/python 
# Written by Samir Singh
# Reads what data is coming in from the fah-work.log and retrieves that data to calculate
# distances between specified atom indices and saves them to a file.
import time
import subprocess
import select
import os
import fnmatch
import mdtraj as md
from bs4 import BeautifulSoup
import numpy as np
from itertools import combinations
from datetime import datetime
from analysis_featurizer import DistanceCalculator, check_file_size

base_dir1 = '/home/server/server2/data/SVR166219/'
base_dir2 = '/array0/data/'
gro_dir = '/home/server/server2/projects/Gromacs/'
MAX_FILE_SIZE = 4000000000
proj_blacklist = []

# reads the log as data comes in
filename = '/home/server/server2/fah-work.log'
f = subprocess.Popen(['tail','-F',filename],\
        stdout=subprocess.PIPE,stderr=subprocess.PIPE)
p = select.poll()
p.register(f.stdout)

while True:
    if p.poll(1):
        line = f.stdout.readline()

	# grabs a line from the fah-work.log that contains the numbers for run, clone, gen
        if 'remaining' in line.split():
	    log_file = open('/home/server/git/fah-scripts/DataAnalysisScripts/readlog.log', 'a+', 0)

	    # isolates the run, clone, and gen number from fah-work.log
            xtc_dir = line.split(':')[4:8]
	    proj_n = xtc_dir[0].lstrip('P')
	    run_n = xtc_dir[1].lstrip('R')
	    clone_n = xtc_dir[2].lstrip('C')
	    gen_n = str(int(xtc_dir[3].lstrip('G'))-1)
	    
	    if proj_n in proj_blacklist:
		log_file.write(datetime.now().strftime('%Y-%m-%d %H:%M:%S: ') + 'The distances file for project %s has reached its size limit\n' % proj_n)
		continue
	    
            # creates a featurizer object from the featurizer class of analysis_featurizer.py
	    try:
	        traj = DistanceCalculator(int(proj_n), int(run_n), int(clone_n), int(gen_n))
	    except Exception as e:
		log_file.write(datetime.now().strftime('%Y-%m-%d %H:%M:%S: ') + str(e) + ' There is most likely no ndx.gro and index file for this project (%s, %s, %s, %s)' % (proj_n, run_n, clone_n, gen_n) + '\n')
		continue
	   # checks if the max file size allowed has been reached and adds the project to a list of projects that will not be calculated if the max size is reached
	    try:
	        print_string, size_lim_reached = check_file_size(os.path.join('/home/server/server2/data/SVR166219/PROJ%s/features/distances.npy' % proj_n), MAX_FILE_SIZE)
	        log_file.write(datetime.now().strftime('%Y-%m-%d %H:%M:%S: ') + print_string + '\n')
	        if size_lim_reached:
		    log_file.write(datetime.now().strftime('%Y-%m-%d %H:%M:%S: ') + 'Adding project %s to project blacklist\n' % proj_n)
		    proj_blacklist.append(proj_n)
	    except Exception as e:
		print(repr(e))
		pass

       	    log_file.write(datetime.now().strftime('%Y-%m-%d %H:%M:%S: ') + 'Calculating features and saving data (%s, %s, %s, %s)...\n' % (proj_n, run_n, clone_n, gen_n))

	    # use built in fucntion from featurizer class that calculates the distances and saves them to a features/distances.npy file in the data directory
	    try:
	    	traj.calculate_distances()
		del traj
	    	log_file.write(datetime.now().strftime('%Y-%m-%d %H:%M:%S: ') + 'Saved data for project %s to index [%s][%s][%s] of distances.npy\n' % (proj_n, run_n, clone_n, gen_n))
	    except Exception as e:
		log_file.write(datetime.now().strftime('%Y-%m-%d %H:%M:%S: ') + 'Unable to save data successfully for project %s: ' % proj_n + repr(e) + '\n')
	    log_file.close()		
	






