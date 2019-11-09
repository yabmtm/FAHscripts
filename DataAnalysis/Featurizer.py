import os, sys, re, glob, subprocess
import time, tqdm, socket
import pandas as pd


import mdtraj as md
import numpy as np
import fnmatch
import itertools
from itertools import groupby, count, combinations
import msmbuilder.utils
from bs4 import BeautifulSoup
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.lines import Line2D
from msmbuilder.cluster import KCenters, KMedoids, KMeans
from msmbuilder.decomposition import tICA
from msmbuilder.featurizer import AtomPairsFeaturizer
from msmbuilder.msm import ContinuousTimeMSM, implied_timescales, MarkovStateModel
from operator import itemgetter
from sklearn.externals import joblib
from sklearn.pipeline import Pipeline

base_dir = '/home/server/server2/data/SVR166219/'
gro_dir = '/home/server/server2/projects/Gromacs/'


class Featurizer:
    """
    This is the base class for calculating features for a frame of data.   
	
    A Featurizer() is an object for featurizing trajectory data for a given FAH project.
    It knows the the 

    Attributes:

        proj_num: The project number
	filename: Name of file where data is stored
	xtc_file: The traj_comp.xtc file for the frame
	gro_file: The .gro file for the run
	traj: The trajectory created by loading in the xtc_file and gro_file together
	data_path: The path to the features file of the project. This path has
		   filename appended to it which will get the full path of the data
		   (default None)
    """

    def __init__(self, proj_num,
                       base_data_dir = '/home/server/server2/data/SVR166219/',
                       base_setup_dir = '/home/server/server2/projects/Gromacs/'):
        """ Constructor for Featurizer class """

	self.proj_num = proj_num
	self.base_data_dir = base_data_dir
	self.dataframe_path = os.path.join(base_dir, 'PROJ' + str(self.proj_num), 'p'+str(self.proj_num)+'.h5')

	# If the *.h5 DataFrame doesn't yet exist, make a new one
	if not os.path.exists(self.dataframe_path):
	    self.create_dataframe()
	else:
	    print('Using existing pandas DataFrame:', self.dataframe_path) 

	### create a file handle to store the data in (the dict-like) HDF5 format
        self.store = pd.HDFStore(self.dataframe_path)

        # store paths for all the setup files
	self.base_setup_dir = base_setup_dir
	assert os.path.exists(self.base_setup_dir)

	self.setup_dir = os.path.join(self.base_setup_dir, 'p'+str(self.proj_num))
	assert os.path.exists(self.setup_dir)

	self.feature_dir = os.path.join(self.setup_dir, 'features')
	assert os.path.exists(self.feature_dir)

	self.config_filename = os.path.join(self.feature_dir, 'config.dat')
	assert os.path.exists(self.feature_dir)

        # The Featurizer class needs to know how many frames per gen, and ns per frame, etc.
	self.nruns = None
	self.nclones = None
	self.ngens = None
	self.ns_per_gen = None
	self.ns_per_frame = None
	self.gmx_index_groups = []     # needed to process the pbc conversion for each run
	self.gmx_index_filenames = []  # needed to process the pbc conversion for each run

	# fill the these attributes using information from the config.dat file
	self.read_config(self.config_filename)

        # create some temp storage for the gro file, and run,
        # so we don't do extra work finding these it and loading it in each time featurize() is called
	self.gro_file = None
        self.run = None

        """
	self.xtc_file = rm_periodic_boundary_cond(self)
	self.traj = md.load(self.xtc_file, top=self.gro_file)
	os.system('rm /home/server/git/fah-scripts/DataAnalysisScripts/temp_xtc/*.xtc')
        """

    def read_config(self, config_filename):
        """Read in the configuration file, and return a dictionary.

        Example features/config.dat file:

        # Configuration file for featurization

        nruns    8
        nclones  1500
        ngens  400

        ns_per_gen     2.50
        ns_per_frame   0.1

        # For pbc removal
        ### Note the index group and filename must be supplied for each run
        pbc_index_info
        28 /home/server/server2/projects/p14188/RUN0/index.ndx
        28 /home/server/server2/projects/p14188/RUN1/index.ndx
        28 /home/server/server2/projects/p14188/RUN2/index.ndx
        28 /home/server/server2/projects/p14188/RUN3/index.ndx
        28 /home/server/server2/projects/p14188/RUN4/index.ndx
        28 /home/server/server2/projects/p14188/RUN5/index.ndx
        28 /home/server/server2/projects/p14188/RUN6/index.ndx
        28 /home/server/server2/projects/p14188/RUN7/index.ndx
        """

	self.nruns = None
	self.nclones = None
	self.ngens = None
	self.ns_per_gen = None
	self.ns_per_frame = None
	self.gmx_index_filenames = []  # needed to process the pbc conversion for each run
	self.gmx_index_groups = []     # needed t

	fin = open(config_filename, 'r')
	lines = fin.readlines()
	fin.close()

	while len(lines) > 0:
		fields = lines[0].strip().split()
		if len(fields) == 0:
			lines.pop(0)
		else:
			if fields[0][0] == '#':
				lines.pop(0)
			elif fields[0] == 'nruns':
				self.nruns = int(fields[1])
				lines.pop(0)
			elif fields[0] == 'nclones':
				self.nclones = int(fields[1])
				lines.pop(0)
			elif fields[0] == 'ngens':
				self.ngens = int(fields[1])
				lines.pop(0)
			elif fields[0] == 'ns_per_gen':
				self.ns_per_gen = float(fields[1])
				lines.pop(0)
			elif fields[0] == 'ns_per_frame':
				self.ns_per_frame = float(fields[1])
				lines.pop(0)
			elif fields[0] == 'pbc_index_info':
				lines.pop(0)
				for i in range(self.nruns):
				    s = lines.pop(0).split()
				    self.gmx_index_groups.append(int(s[0]))
				    self.gmx_index_filenames.append(s[1])

    def create_dataframe(self):
        """Create an empty DataFrame for this project."""

        df = pd.DataFrame({'date': [],
                           'RUN': [],
                           'CLONE': [],
                           'GEN': pd.Series(0, index=[], dtype='int'),
                           'frame': pd.Series([], index=[], dtype='int'),
                           'time (ns)': [] })  # set the index
        df.set_index('date')
        print(df)

        # Save the DataFrame to disk

        ### create a file handle to store the data in (a dict-like) HDF5 format
        store = pd.HDFStore(self.dataframe_path)
        print(store)
        store.put('df', df)
        return store


    def featurize(self, run, clone, gen):
        """
        Featurizes the trajectory data for the given run, clone, gen
            
        INPUT
        
        run: The run number
        clone: The clone number
        gen: The gen number).
        """

        # Find the gro file structure template, only if it's a new RUN.
        # Otherwise we can just use the old one
        if run != self.run:
            self.gro_file = self.find_gro(self.setupdir, run)
            self.run = run

        # Get the xtc file trajectory for this RUN, CLONE, GEN
        xtc_filename = self.find_traj()

        # Remove the periodic boundary conditions in a new temporary file
        xtc = self.rm_periodic_boundary_cond(xtc_filename)

        ### DO THE FEATURIZATION HERE ###
        pass

        # Close the HDF5 handle to the DataFrame
        self.store.close()

        # Close the file handle to the pbc-processed xtc
        xtc.close()

    def rm_periodic_boundary_cond(self, filename):
        """ 
        Creates a new trajectory file with the periodic boundary conditions 
        removed. It saves the file in a directory where it should be deleted
        after this function is called and data is taken from the file

        INPUT:
        filename: the filename of the xtc to process

        Returns:
        xtc: a tempfile.NamedTemporaryFile() file handle for the processed
            trajectory.  xtc.name is the filename.
        """

        import tempfile
        tmp = tempfile.NamedTemporaryFile(suffix='.xtc')
        print(tmp.name)


        f = open('/home/server/server2/projects/Gromacs/p' + str(feat.proj_num) + '/features/index')

        list_index = f.read().split('\n')
        indx1, indx2 = int(list_index[0]), int(list_index[1])
        os.system('echo "%d\n%d\n" | gmx trjconv -f /home/server/server2/data/SVR166219/PROJ%d/RUN%d/CLONE%d/results%s/traj_comp.xtc -s \
				  /home/server/server2/data/SVR166219/PROJ%d/RUN%d/CLONE%d/frame0.tpr -n /home/server/server2/projects/Gromacs/p%d/features/xtc.ndx \
				  -pbc mol -center -o /home/server/git/fah-scripts/DataAnalysisScripts/temp_xtc/traj_%d_%d_%d_%d.xtc &>/dev/null'
                  % (
                  indx1, indx2, feat.proj_num, feat.run_num, feat.clone_num, feat.gen_num, feat.proj_num, feat.run_num,
                  feat.clone_num, feat.proj_num, feat.proj_num,
                  feat.run_num, feat.clone_num, feat.gen_num))
        new_xtc = '/home/server/git/fah-scripts/DataAnalysisScripts/temp_xtc/traj_%d_%d_%d_%d.xtc' % (
        feat.proj_num, feat.run_num, feat.clone_num, feat.gen_num)
        return new_xtc


    def find_traj(self, run, clone, gen):
        """Searches for the trajectory xtc file of a certain run, clone, gen.	
           
        INPUTS
        run, clone, gen

        Returns:
            xtc_file: The directory of the xtc file.
            None: If the file with the given parameters does not exist
        """

        xtc_dir = 'PROJ%d/RUN%d/CLONE%d/results%d/traj_comp.xtc' % (self.proj_num, run, clone, gen)
        xtc_file = os.path.join(base_dir, xtc_dir)
        if os.path.isfile(xtc_file):
            return xtc_file
        else:
            return None #return None if there is no traj_comp.xtc


    def find_gro(self, gro=None):
        """
	Searches for a compatible .gro file for a certain frame.
	  
	Parameters:
		gro: Name of .gro file (default None).
	
	Returns:
		gro_file: Directory of the .gro file.
		None: If no compatible .gro file was found
        """

	if gro == None:
            gro_file_dir = os.path.join(gro_dir, 'p' + str(self.proj_num))
	    for filename in os.listdir(gro_file_dir):
		if filename.endswith('.gro'):
		    try:
			gro_file = os.path.join(gro_file_dir, filename)
			t = md.load(self.xtc_file, top=gro_file)
			return gro_file
		    except:
			pass
		else:
		    pass

		#if not found yet, searches in the given RUN file
		gro_file_dir = os.path.join(gro_file_dir, 'RUN' + str(self.run))
		files = []
		for root, dirnames, filenames in os.walk(gro_file_dir):
		    for filename in fnmatch.filter(filenames, '*.gro'):
			files.append(os.path.join(root, filename))
		for filename in files:
		    file_n = filename.split('/')[-1]
		    try:
			filename = os.path.join(gro_file_dir, filename)
			t = md.load(self.xtc_file, top=filename)
			return filename
		    except:
			pass
	else:
	    f1 = os.path.join(gro_dir, 'p' + str(self.proj_num), gro)
	    f2 = os.path.join(gro_dir, 'p' + str(self.proj_num), 'RUN' + str(self.run_num), gro)
	    try:
		t = md.load(self.xtc_file, top=f1)
	        return f1
	    except:
		try:
		    t = md.load(self.xtc_file, top=f2)
		    return f2
		except:
		    pass
    	    # return None if no compatible gro file found
	    return None


    def load_npy_file(self):
	"""
	Loads the .npy file from the features directory.
	
	Returns:
	    npy_file: The data that is stored in the .npy file
	    None: If there was no .npy file of filename found
	"""

	if os.path.isfile(self.data_path):
	    npy_file = np.load(self.data_path)
	    return npy_file
	else:
	    print('Could not find specified file')
	    return None


