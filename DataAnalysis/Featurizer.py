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
			pass

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



		self.xtc_file = rm_periodic_boundary_cond(self)

		self.traj = md.load(self.xtc_file, top=self.gro_file)


		os.system('rm /home/server/git/fah-scripts/DataAnalysisScripts/temp_xtc/*.xtc')

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
					self.ns_per_frame = int(fields[1])
					lines.pop(0)
				elif fields[0] == 'pbc_index_info':
					for i in range(self.nruns):
						s = lines.pop(0).split()
						self.gmx_index_groups.append(int(s[0]))
						self.gmx_index_filenames.append(s[1])
					lines.pop(0)

	def create_dataframe(self):
		"""Create an empty DataFrame for this project."""

        df = pd.DataFrame({'date': [],
                           'RUN': [], dtype='int'),
                           'CLONE': [], dtype='int'),
                           'GEN': pd.Series(0, index=l[], dtype='int'),
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


    def find_traj(self):
            """ 
            Searches for the trajectory xtc file of a certain frame.	
            
            Parameters:
            feat: Featurizer object
        
            Returns:
            xtc_file: The directory of the xtc file.
            None: If the file with the given parameters does not exist
            """

            xtc_dir = 'PROJ%s/RUN%s/CLONE%s/results%s/traj_comp.xtc' % (str(feat.proj_num), str(feat.run_num), str(feat.clone_num), str(feat.gen_num))
            xtc_file = os.path.join(base_dir, xtc_dir)
            if os.path.isfile(xtc_file):
                return xtc_file
            else:
                return None #return None if there is no traj_comp.xtc


	def find_gro(self, gro=None):
		"""
		Searches for a compatible .gro file for a certain frame.
	  
		Parameters:
		feat: Featurizer object
		gro: Name of .gro file (default None).
	
		Returns:
		gro_file: Directory of the .gro file.
		None: If no compatible .gro file was found
		"""

		if gro == None:
		gro_file_dir = os.path.join(gro_dir, 'p' + str(feat.proj_num))
		for filename in os.listdir(gro_file_dir):
			if filename.endswith('.gro'):
			try:
				gro_file = os.path.join(gro_file_dir, filename)
				t = md.load(feat.xtc_file, top=gro_file)
				return gro_file
			except:
				pass
			else:
			pass

			#if not found yet, searches in the given RUN file
		gro_file_dir = os.path.join(gro_file_dir, 'RUN' + str(feat.run_num))
		files = []
		for root, dirnames, filenames in os.walk(gro_file_dir):
			for filename in fnmatch.filter(filenames, '*.gro'):
			files.append(os.path.join(root, filename))
		for filename in files:
			file_n = filename.split('/')[-1]
			try:
			filename = os.path.join(gro_file_dir, filename)
			t = md.load(feat.xtc_file, top=filename)
			return filename
			except:
			pass
		else:
		f1 = os.path.join(gro_dir, 'p' + str(feat.proj_num), gro)
		f2 = os.path.join(gro_dir, 'p' + str(feat.proj_num), 'RUN' + str(feat.run_num), gro)
		try:
			t = md.load(feat.xtc_file, top=f1)
			return f1
		except:
			try:
			t = md.load(feat.xtc_file, top=f2)
			return f2
			except:
			pass
		# return None if no compatible gro file found
		return None


	def load_npy_file(self):
		"""
		Loads the .npy file from the features directory.
	
		Parameters:
		feat: Featurizer object
	
		Returns:
		npy_file: The data that is stored in the .npy file
		None: If there was no .npy file of filename found
		"""

		if os.path.isfile(feat.data_path):
			npy_file = np.load(feat.data_path)
		return npy_file
		else:
		print('Could not find specified file')
		return None


	def create_feat_dir(self):
		""" Creates a features directory in the data directory of the given Featurizer object """

		os.mkdir(os.path.join(base_dir, 'PROJ' + str(feat.proj_num), 'features'))


	def create_npy_file(self, default=True):
		"""
		Creates a .npy file in the features directory of the project of the shape
		of the number of runs, clones, and gens in the project.
	
		Parameters:
		feat: Featurizer object
			default: If False will create the file initialized with zeros 
			 (default True: Initializes it as an empty array)
		"""

		feat_dir = os.path.join(base_dir, 'PROJ' + str(feat.proj_num), 'features')
		if os.path.isfile(feat.data_path):
		print('File already exists')
		# parses project.xml file of the project to find the number of runs, clones, and gens of the projects
		elif os.path.isdir(feat_dir):
		project_xml = open(os.path.join(gro_dir, 'p' + str(feat.proj_num), 'project.xml'), 'r')
		xml = project_xml.read()
		soup = BeautifulSoup(xml, 'lxml')
		if default == True:
			arr = np.empty((int(soup.find('runs')['v']), int(soup.find('clones')['v']), int(soup.find('gens')['v'])), dtype=np.ndarray)
		else:
			arr = np.zeros((int(soup.find('runs')['v']), int(soup.find('clones')['v']), int(soup.find('gens')['v'])))
		np.save(feat.data_path, arr)
		# if a features file does not exist, creates one and runs this function again
		else:
			print('features file does not exist for this project. Creating one now...')
		create_feat_dir(feat)
		create_npy_file(feat, default)

	def load_and_save(feat, data):
		"""
		Loads a .npy file, writes data to its specified index, and saves it.
	
		Parameters:
		feat: Featurizer object
		data: The data to be saved
		"""

		try:
		npy_file = np.load(feat.data_path)
		npy_file[int(feat.run_num), int(feat.clone_num), int(feat.gen_num)] = data
		np.save(feat.data_path, npy_file)
		except Exception as e:
		print(repr(e))

	def check_file_size(filename, max_size=sys.maxint):
		""" Checks the size of the given file and will stop a script if a file reaches the max_size parameter. """

		file_size = os.path.getsize(filename) # strips last character which is the letter 'l'
		if file_size > max_size:
		string = 'File %s size has reached its user specified limit' % filename.split('/')[-1]
		return string, True
		else:
		file_size_mb = file_size/1000000
		max_size_mb = max_size/1000000
		pct_used = float(file_size)/max_size * 100
		size_info = 'The file %s has used %dmb/%dmb (%.2f%%) of the max user specified file size' % (filename.split('/')[-1], file_size_mb, max_size_mb, pct_used)
		return size_info, False


class DistanceCalculator(Featurizer):

	"""
	A child of the Featurizer class for calculating Distances between atom indices.

	Attributes:
		proj_num: The project number
		run_num: The run number
		clone_num: The clone number
		gen_num: The gen number
	filename: Name of the file that the data will be saved to (default 'distances.npy')
		xtc_file: The traj_comp.xtc file for the frame
		gro_file: The .gro file for the run
		traj: The trajectory created by loading in the xtc_file and gro_file together
	features_file: the file path to features.npy which holds atom indices for calculation
	data_path: The path to filename where all the distances for that project
		   is stored
	"""

	def __init__(self, proj_num, run_num, clone_num, gen_num, filename='distances.npy'):
		""" Constructor for the DistanceCalculator class """

		Featurizer.__init__(self, proj_num, run_num, clone_num, gen_num, filename)
		self.features_file = os.path.join(gro_dir, 'p' + str(proj_num), 'features', 'features.npy')
		if not os.path.isfile(self.data_path):
		   create_npy_file(self)

	def get_distance_indices(self, atom_desc=None):
		""" 
	Returns atom indices for calculating distances
	
	parameters:
		atom_desc: Will retrieve atom indices with description atom_desc
			   (default None: will check for features.npy and if that
			does not exit it will return alpha and beta carbons)

	returns:
		indice_pairs: An itertools list of all combinations of atom indices
	 """
	print('Loading in atom indices...')
	if atom_desc is not None:
		indices = self.traj.topology.select(atom_desc)
			indice_pairs = combinations(indices, 2)
			return indice_pairs
	if os.path.isfile(self.features_file):
		indice_pairs = np.load(self.features_file)
			return indice_pairs
	else:
		indices = self.traj.topology.select('name CA or name CB')
			indice_pairs = combinations(indices, 2)
		return indice_pairs

	def calculate_distances(self, indice_pairs=None, save=True):

	""" 
	Using atom indices from get_distance_indices(), computes the distances between
	those indices and saves them to distances.npy.
	
	Parameters:
		save: Saves the calculated distances if save is True (default True)

	Returns:
		distances: A numpy array of the computed distances
	"""
	if indice_pairs == None:
		indice_pairs = self.get_distance_indices()
	print('Calculating distances...')
	distances = md.compute_distances(self.traj, indice_pairs)
	if save:
			load_and_save(self, distances)
		print('Distances saved for index [%d][%d][%d]' % (self.run_num, self.clone_num, self.gen_num))
	return distances

	def load_distances(self):
	""" Loads and returns the distances calculated for this gen. """

	return np.load(self.data_path)[self.run_num][self.clone_num][self.gen_num]


class DihedralCalculator(Featurizer):

	"""
	A child class of Featurizer used for calculating dihedral angles.

	Attributes:
		proj_num: The project number
		run_num: The run number
		clone_num: The clone number
		gen_num: The gen number
	angle: The angle you want calculated
	filename: Name of file where data is stored in
		xtc_file: The traj_comp.xtc file for the frame
		gro_file: The .gro file for the run
		traj: The trajectory created by loading in the xtc_file and gro_file together
		data_path: The path to filename which holds all the dihedral data for that 
		   angle of the whole project (default None: if None makes dihedrals_(angle).npy)
	"""
	
	def __init__(self, proj_num, run_num, clone_num, gen_num, angle, filename=None):
	""" Constructor for DihedralCalculator class """
	if filename is None:
		filename = 'dihedrals_%s.npy' % angle
	Featurizer.__init__(self, proj_num, run_num, clone_num, gen_num, filename)
	self.angle = angle
	if not os.path.isfile(self.data_path):
		create_npy_file(self)

	def calculate_dihedrals(self, save=True):
	""" 
	Calculates dihedral angles using mdtraj functions given omega, psi, or phi

	Parameters:
		save: Saves to dihedrals_(angle).npy if True (default True)
	
	Returns:
		indices: The indices that dihedrals are being calculated for
		dihedrals: The calculated dihedral angles 	
	"""
	print('Calculating dihedrals for angle %s...' % self.angle) 
	if self.angle == 'phi':
		indices, dihedrals = md.compute_phi(self.traj)
	elif self.angle == 'psi':
		indices, dihedrals = md.compute_psi(self.traj)
	else:
		indices, dihedrals = md.compute_omega(self.traj)
	if save:
		load_and_save(self, dihedrals)
		print('Dihedrals saved for index [%d][%d][%d]' % (self.run_num, self.clone_num, self.gen_num))
	return indices, dihedrals

	def load_dihedrals(self):
	""" Loads in the dihedral data for that gen and returns it """
	return np.load(self.data_path)[self.run_num][self.clone_num][self.gen_num]


class COMCalculator(Featurizer):

	"""
	A child class of Featurizer used for calculating Center of Mass distances.

	Attributes:
		proj_num: The project number
		run_num: The run number
		clone_num: The clone number
		gen_num: The gen number
	filename: Name of file where data will be stored in (default 'COM.npy')
		xtc_file: The traj_comp.xtc file for the frame
		gro_file: The .gro file for the run
		traj: The trajectory created by loading in the xtc_file and gro_file together
		data_path: The path to filename which holds all COM data for the project
	"""
   

	def __init__(self, proj_num, run_num, clone_num, gen_num, filename='COM.npy'):
	""" Constructor for COMCalculator class """
	Featurizer.__init__(self, proj_num, run_num, clone_num, gen_num, filename)
	if not os.path.isfile(self.data_path):
		create_npy_file(self)

	def compute_center_of_mass(self, traj=None, atom_indices=None):

		"""
	Compute the center of mass for each frame.

		Parameters:
			traj : Trajectory to compute center of mass for
			atom_indices : Atoms to compute center of mass for. 
		If None, will compute over all atoms

		Returns:     
			com : np.ndarray, shape=(n_frames, 3)
			 Coordinates of the center of mass for each frame
		"""
	if traj == None:
			traj = self.traj
		if atom_indices is None:
			atoms = traj.top.atoms
			coords = traj.xyz
		else:
			atoms = [traj.top.atom(i) for i in atom_indices]
			coords = np.take(traj.xyz, atom_indices, axis=1)

		com = np.zeros((traj.n_frames, 3))
		masses = np.array([a.element.mass for a in traj.top.atoms])
		masses = np.array([a.element.mass for a in atoms])
		masses /= masses.sum()

		#for i, x in enumerate(traj.xyz):
		for i, x in enumerate(coords):
			com[i, :] = x.astype('float64').T.dot(masses)
		return com

	def COM_atom_distance(self, traj=None, atom_indices=None, atom_index=None):
	"""
	Computes distance between the center of mass of the protein and ligand

	Parameters:
		traj: Trajectory for computing center of mass (default None: 
		   uses trajectory used by the class)
		atom_indices: the atom indices that ceneter of mass is calculated for
							(default None: evaluates a features file put
				 in the project directory by the user)
		atom_index: A single atom index that will be used to calculate
			the distance from the center of mass 
			(default None: reads in index from features file)
	Returns:
		distance_between_groups: The distances between ligand_indices
					 and atom_index
	"""	
	if traj == None:
		traj = self.traj

	if atom_indices == None or atom_index == None:
		feat_txt_dir = os.path.join(gro_dir, 'p' + str(self.proj_num), 'features', 'US')
		feat_txt = open(feat_txt_dir)
		feats = feat_txt.read()
		feats_spt = feats.split(',', 1)
		if atom_index == None:
		atom_index = int(feats_spt[0][1:])
		if atom_indices == None:
		atom_indices = eval(feats_spt[1][:-2])
		
		atoms_COM = self.compute_center_of_mass(traj, atom_indices)
		atom_xyz = np.take(traj.xyz, [atom_index], axis=1)

		for k in range(len(traj)):
			distance_between_groups = np.sqrt((atoms_COM[k][0] - atom_xyz[k][0][0])**2 +
				(atoms_COM[k][1] - atom_xyz[k][0][1])**2 +
				(atoms_COM[k][2] - atom_xyz[k][0][2])**2)
	return distance_between_groups

	def COM_atom_distances_for_gen(self, traj=None, save=True, atom_indices=None, atom_index=None):
	"""
	Uses COM_atom_distances to calculate COM distance over all frames in a gen
	
	Parameters:
		traj: Trajectory file (default None: uses the class's loaded in traj file)
		save: Saves to COM.npy if True (defualt True)
		atom_indices: Atom indices for center of mass (default None)
		atom_index: Atom index to compute distance to from center of mass (default None)

	Returns:
		arr: and array of all the COM to atom distances in that gen
	"""
	if traj == None:
			traj = self.traj
	arr = []
	for i in range(self.traj.n_frames):
		frame = traj[i]
		print('Calculating COM for frame %d of gen %d...' % (i, self.gen_num))
		data = self.COM_atom_distance(frame, atom_indices, atom_index)
		arr.append(data)
	if save:
		load_and_save(self, arr)
		print('COM saved to index [%d][%d][%d]' % (self.run_num, self.clone_num, self.gen_num))
	return arr	

	def load_COM(self):
	""" Returns the data located in the index of the specific gen """
	return np.load(self.data_path)[self.run_num][self.clone_num][self.gen_num]



class RMSDCalculator(Featurizer):
	"""
	This class is a child of the Featurizer class and is for calculating RMSD
	with a reference and target trajectory.

	Attributes:
		proj_num: The project number
		run_num: The run number
		clone_num: The clone number
		gen_num: The gen or results number
	filename: Name of file where data will be stored in (default 'RMSD.npy')
		proj_distances: Data of distances calculated for the project stored in distances.npy
		distances: Data of distances calculated for the frame
	data_path: path to filename which holds all RMSD calculated for the project
	"""
 
	def __init__(self, proj_num, run_num, clone_num, gen_num, filename='RMSD.npy'):
	""" Constructor for RMSDCalculator class """

	Featurizer.__init__(self, proj_num, run_num, clone_num, gen_num, filename)
	if not os.path.isfile(self.data_path):
		create_npy_file(self)
	
	def get_indices(self, atom_desc=None):
		""" 
		Returns atom indices for calculating distances
		
		parameters:
			atom_desc: Will retrieve atom indices with description atom_desc
					   (default None: will check for features.npy and if that
						does not exit it will return alpha and beta carbons)

		returns:
			indice_pairs: An itertools list of all combinations of atom indices

		 """
		print('Loading in atom indices...')
	feat_npy = '/home/server/server2/data/projects/Gromacs/p%d/features' % self.proj_num
	if atom_desc is not None:
		indices = self.traj.topology.select(atom_desc)
			indice_pairs = combinations(indices, 2)
			return indice_pairs
		elif os.path.isfile(feat_npy):
			indice_pairs = np.load(feat_npy)
			return indice_pairs
		else:
			indices = self.traj.topology.select('name CA or name CB')
			indice_pairs = combinations(indices, 2)
			return indice_pairs


	def calculate_RMSD(self, ref_struct=None, atom_desc=None, save=True):
	""" 
	Calculates RMSD from the reference of ref_struct

	Parameters:
		ref_struct: The reference structure file (default None: Uses the structure 
			file from the first trajectory)
		atom_desc: Description of atom indices (default None)
		save: Saves data to RMSD_data.npy if save is True (default True)	    

	Returns:
		RMSD_list: List of all RMSD values calculated for the gen
	"""

	gro_ref = md.load(self.gro_file)
	indice_pairs = self.get_indices(atom_desc)
	print('Calculating RMSD...')
	RMSD = md.rmsd(self.traj, gro_ref)
	if save:
		load_and_save(self, RMSD)
		print('Saved RMSD to index [%d][%d][%d]' % (self.run_num, self.clone_num, self.gen_num))
	return RMSD

	def load_RMSD(self):
		""" Loads RMSD data for the gen """

		return np.load(self.data_path)[self.run_num][self.clone_num][self.gen_num]


class Analysis(Featurizer):
	"""    
	The Analysis class is a child of the Featurizer class that is used for 
	analyzing calculated data

	Attributes:
	proj_num: The project number
		run_num: The run number
		clone_num: The clone number
		gen_num: The gen or results number
		proj_distances: Data of distances calculated for the project stored in distances.npy
		filename: The name of the .npy file which the data is stored in (default None)
		self.data: If a filename is specified, this will hold the full data for that project
	self.data_run: This holds all the data in the run number
	self.data_clone: This holds all the data in the clone number
	self.data_gen: This holds all the data in the gen number
	"""

	def __init__(self, proj_num, run_num, clone_num, gen_num, filename=None):
	""" Constructor for Analysis class """

		Featurizer.__init__(self, proj_num, run_num, clone_num, gen_num, filename)
		self.data = np.load(self.data_path)
		self.data_run = self.data[self.run_num]
		self.data_clone = self.data_run[self.clone_num]
		self.data_gen = self.data_clone[self.gen_num]
	plot_dir = '/home/server/server2/data/SVR166219/PROJ14137/plots'
	if not os.path.isdir(plot_dir):
		os.mkdir(plot_dir)

class COMAnalysis(Analysis):
	"""    
	The COMAnalysis class is a subclass of the analysis class that
	is used for analyzing Center of Mass data

	Attributes:
		proj_num: The project number
		run_num: The run number
		clone_num: The clone number
		gen_num: The gen or results number
		proj_distances: Data of distances calculated for the project stored in distances.npy
		filename: The name of the .npy file which the data is stored in (default COM.npy)
		self.data: This will hold the full data for that project
		self.data_run: This holds all the data in the run number
		self.data_clone: This holds all the data in the clone number
		self.data_gen: This holds all the data in the gen number
	"""


	def __init__(self, proj_num, run_num, clone_num, gen_num, filename='COM.npy'):
	""" Constructor for COMAnalysis class """

		Analysis.__init__(self, proj_num, run_num, clone_num, gen_num, filename)

	def plot_COM_atom_distances(self):
	""" plots a line graph of COM distances for a whole clone """

		try:
			arr = []
			plt.figure()
			plt.ylabel('COM Distance')
			plt.xlabel('Frame')
			for i, gen_data in enumerate(self.data_clone):
				if gen_data is not None:
						for item in gen_data:
								arr.append(item)
			plt.plot(arr)
		#plt.ylim(.5, 1.5)
			plt.savefig(os.path.join(base_dir, 'PROJ' + str(self.proj_num), 'plots',
										 'COM_distances_proj%d_run%d_clone%d.png'
										 % (self.proj_num, self.run_num, self.clone_num)))
			print('Plot made for run %d, clone %d of project %d' % (self.run_num, self.clone_num, self.proj_num))
			plt.close()
		except Exception as e:
			print(repr(e) + '\nThis function is meant to be run if all COM distances have been calculated for the clone this gen is in')



class RMSDAnalysis(Analysis):
	"""    
	The RMSDAnalysis class is a subclass of the Analysis class used
	for analyzing RMSD data

	Attributes:
		proj_num: The project number
		run_num: The run number
		clone_num: The clone number
		gen_num: The gen or results number
		proj_distances: Data of distances calculated for the project stored in distances.npy
		filename: The name of the .npy file which the data is stored in (default RMSD_data.npy)
		self.data: This will hold the full data for that project
		self.data_run: This holds all the data in the run number
		self.data_clone: This holds all the data in the clone number
		self.data_gen: This holds all the data in the gen number
	"""


	def __init__(self, proj_num, run_num, clone_num, gen_num, filename='RMSD.npy'):
	""" Constructor for RMSDAnalysis class """

	Analysis.__init__(self, proj_num, run_num, clone_num, gen_num, filename)

	def plot_RMSD(self, num_bins=500):
	""" 
	Plots all RMSD values of the project in a histrogram
	
	Parameters:
		num_bins: Number of desired bins in the histogram (default 500)
	"""

	plot_data = []
	for run, arr2d in enumerate(self.data):
		for clone, arr1d in enumerate(arr2d):
		for gen, nparray in enumerate(arr1d):
			if nparray is not None:
			for num in nparray:
				plot_data.append(num)
	RMSD_hist = plt.hist(plot_data, bins=num_bins)
	plt.ylabel('Bin Height')
	plt.xlabel('RMSD')
	plt.savefig(os.path.join(base_dir, 'PROJ' + str(self.proj_num), 'plots', 'RMSD_plot_' + str(self.proj_num) + '.png'))
	print('Plotted all RMSD values available for this project')

class TICA(Analysis):

	def __init__(self, proj_num, run_num, clone_num, gen_num, filename='distances.npy'):
	Analysis.__init__(self, proj_num, run_num, clone_num, gen_num, filename)
	self.tICA_dir= base_dir + 'PROJ%s/tICA' % self.proj_num
	if not os.path.isdir(self.tICA_dir):
		os.makedirs(self.tICA_dir)

	def calculate_stride_distances(self, stride, equil_steps):
	t = DistanceCalculator(self.proj_num, 0, 0, 0)
	indices = t.get_distance_indices()
	runs, clones, gens = self.data.shape[0], self.data.shape[1], self.data.shape[2]
	npy_path = '/home/server/git/fah-scripts/DataAnalysisScripts/stride_dist/stride_dist_%d.npy' % self.proj_num
	if not os.path.isfile(npy_path):
		arr = np.empty((runs, clones, gens), dtype=np.ndarray)
		np.save(npy_path, arr)
	arr = np.load(npy_path)
	for i in range(runs):
		for j in range(clones):
		for k in range(gens):
			if arr[i][j][k] is None and os.path.isdir('/home/server/server2/data/SVR166219/PROJ%d/RUN%d/CLONE%d/results%d' % (self.proj_num, i, j, k)):
				distances = []
				x = 0
				while True:
				try:
				if k == 0 and j == 0 and i == 0 and  x < equil_steps:
					x += stride
					continue
				t = DistanceCalculator(self.proj_num, i, j, k)
						distances.append(md.compute_distances(t.traj[x], indices))
				del t
				except:
					break
				x += stride
			if len(distances) > 0:
					arr[i][j][k] = np.asarray(np.concatenate(distances))
				print('Calculated distances of index[%d][%d][%d]' % (i, j, k))
			else:
			break
		print('Distances saved')
			np.save(npy_path, arr)


	

	def calculate_tica_components(self, cluster_method, calculate_strides=False, feats=None):

		'''Load in the features, calculate a given number of tICA components (tica_components) given a
		lagtime (lag_time), and save tICA coordinates and eigenvector data. It then creates and populates
		a list for each desired component, clusters the data, saving normalized populations as populations.dat
		and saving each cluster center as a .pdb. tICA plots are created and saved, and implied timescales are
		calculated, saved, and plotted.
		'''
	
	# tICA parameters
		tica_lagtime = 10 # determine from implied timescales
		tica_components = 8 # how many tICs to compute
		n_clusters = 100 # denotes number of microstates
		n_timescales = tica_components # plot all eigenvalues --> timescales
		md_time_step = 0.02 # ns
		subsampled_time_step = 1. # ns multiplier of timescales and lagtimes in implied timescale plot
		stride = int(subsampled_time_step / md_time_step)  #time step stride for sub-sampling
	equil_time = 1. # ns
		equil_steps = 1 #int(equil_time / md_time_step)  time steps to be removed from start
		lagtimes = np.array([1,2,4,8,16,32,64,128,256,512,1024])
		cluster_method = 'kcenters' # 'kcenters/kmeans'
		all_ticas = list(itertools.permutations(range(1,tica_components+1), 2)) # all combinations
		all_ticas = [[1,2]] # override: just show analysis for first two components
		cluster_percentage_cutoff = 5 # clusters with a relative population less than this
								  # number will not be labeled on plot i.e. 0 : all clusters labeled
		verbose = False

		print("\nCalculating tICA components...")

		# Load in feature files THIS WILL NEED TO BE CHANGED
	if feats == None:
		if calculate_strides:
				self.calculate_stride_distances(stride, equil_steps)
			data = np.load('/home/server/git/fah-scripts/DataAnalysisScripts/stride_dist/stride_dist_%d.npy' % self.proj_num)
		else:
			data = self.data
	else:
		data = np.load(feats)

	features = []
	for run in data:
		for clone in run:
			gen_seq = []
			for gen in clone:
				if gen is not None and gen[0] is not None:
						if calculate_strides or feats is not None:
				gen_seq.append(gen)
				else:
						gen_seq.append(gen[::stride])
			if len(gen_seq) > 0:
			gen_cat = np.concatenate(gen_seq)
			if calculate_strides:
				features.append(gen_cat)
			else:
				features.append(gen_cat[equil_steps:])
	features = np.asarray(features)
	print(features.shape)
	print(features[0].shape)
	tica_coordinates = tICA(lag_time=tica_lagtime,
			n_components=int(tica_components)).fit_transform(features)

		np.save('%s/lag_%d_coord_%d.npy' %(self.tICA_dir, tica_lagtime, tica_components), tica_coordinates)

		# Initiate and populate an array for each component
		for i in range(tica_components):
			exec('tica_' + str(i+1) + ' = []')

		for i in tqdm.tqdm(range(len(features))):
			for j in range(len(tica_coordinates[i])):
				for k in range(tica_components):
					exec('tica_' + str(k+1) + '.append(tica_coordinates[i][j][k])')

		# Perform clustering based on the cluster_method parameter.
		if cluster_method == 'kcenters':
			print("Clustering via KCenters...")
			clusters = KCenters(n_clusters)
		elif cluster_method == 'kmeans':
			print("Clustering via KMeans...")
			clusters = KMeans(n_clusters)
		else:
			sys.exit("Invalid cluster_method. Use kmeans or kcenters.")

		# Determine cluster assignment for each frame.
		sequences = clusters.fit_transform(tica_coordinates)
	
		np.save('%s/lag_%d_clusters_%d_sequences.npy' %(self.tICA_dir, tica_lagtime, n_clusters), sequences)
		np.save('%s/lag_%d_clusters_%d_center.npy' %(self.tICA_dir, tica_lagtime, n_clusters),
		clusters.cluster_centers_)

		# Determine cluster populations, normalize the counts, and save as percentages for
		# labeling if a cluster contains more than cluster_percentage_cutoff percent of the data.
		# Finally, save normalized counts.
		print("\nDetermining cluster populations...")

		if not os.path.exists('%s/cluster_centers' % self.tICA_dir):
			os.makedirs('%s/cluster_centers' % self.tICA_dir)
		counts = np.array([len(np.where(np.concatenate(sequences)==i)[0]) for i in range(n_clusters)])
		normalized_counts =  counts/float(counts.sum())
		percentages = [ i*100 for i in normalized_counts ]
		population_labels = [ [i,"%.2f"%percentages[i]] for i in range(len(percentages)) if percentages[i] > cluster_percentage_cutoff ]
		np.savetxt('%s/cluster_centers/populations.dat' % self.tICA_dir, normalized_counts)
	

		# Plot all unique combinations of tICA components
		print("\nPlotting tICA components with cluster centers...")
		all_ticas = list(itertools.permutations(range(1,tica_components+1), 2))
		for j in tqdm.tqdm(range(len(all_ticas))): # For each pair
			if all_ticas[j][0] < all_ticas[j][1]:
				plt.figure(j, figsize=(20,16))
				plt.hexbin(eval("tica_"+str(all_ticas[j][0])), eval("tica_"+str(all_ticas[j][1])), bins='log')
				x_centers = [clusters.cluster_centers_[i][all_ticas[j][0]-1] for i in range(len(clusters.cluster_centers_))]
				y_centers = [clusters.cluster_centers_[i][all_ticas[j][1]-1] for i in range(len(clusters.cluster_centers_))]
				high_pop_x_centers = [ x_centers[i] for i in range(len(x_centers)) if percentages[i] > cluster_percentage_cutoff ]
				high_pop_y_centers = [ y_centers[i] for i in range(len(y_centers)) if percentages[i] > cluster_percentage_cutoff ]
				plt.plot(x_centers, y_centers, color='y', linestyle="", marker="o")
				plt.plot(eval("tica_"+str(all_ticas[j][0])+'[0]'), eval("tica_"+str(all_ticas[j][1])+'[0]'), color='k', marker='*',markersize=24)
				plt.xlabel('tic'+str(all_ticas[j][0]))
				plt.ylabel('tic'+str(all_ticas[j][1]))
				plt.title(self.proj_num)
				# Add labels for high-population cluster centers
				for label, x, y in zip(population_labels, high_pop_x_centers, high_pop_y_centers):
					plt.annotate(
					  label,
					  xy = (x, y), xytext = (-15, 15),
					  textcoords = 'offset points', ha = 'right', va = 'bottom',
					  bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
					  arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
				plt.savefig('%s/tica_' % (self.tICA_dir) +str(all_ticas[j][0])+'_'+str(all_ticas[j][1])+'.png')
				plt.close()

###########################################################################
	for filename in os.listdir(self.tICA_dir + '/cluster_centers'):
		if filename.endswith('.pdb'):
				os.remove(self.tICA_dir + '/cluster_centers/' + filename)
	# Write out PDBs for each cluster center
		print("Performing cluster analytics and saving center PDBs...\n")
	runs, clones, gens = data.shape[0], data.shape[1], data.shape[2]
	x, y, z = 0, 0, 0
	for i in range(len(features)):
		if i % clones == 0 and i != 0:
		x += 1
		if i % gens == 0:
		y = 0
			n_snapshots = len(clusters.distances_[i])

			# Determine frames that are cluster centers
			cluster_indices = np.arange(n_snapshots)[ (clusters.distances_[i] < 1e-6) ]
			# Determine number of each cluster, correlates to populations.dat
			cluster_labels = sequences[i][cluster_indices]
			# Save each cluster center as a pdb
			if list(cluster_indices): # load center-containing xtcs to check length
		traj_cat = []
		print('x: %d, y: %d, z: %d' % (x, y, z))

		while True:
			try:
			traj = base_dir + 'PROJ%s/RUN%s/CLONE%s/results%s/traj_comp.xtc' % (self.proj_num, x, y, z)
					traj_cat.append(md.load(traj, top=self.gro_file))
			z += 1
			except:
			break
		if len(traj_cat) > 0:
			trajectory_file = md.join(traj_cat)
				xtc_len = len(trajectory_file)
		y += 1
			z = 0
			for j in range(len(cluster_indices)):
				frames = range(xtc_len) # map the strided frame number back to xtc frame number
				strided_frames = frames[equil_steps:][::stride]
				xtc_frame = frames.index(strided_frames[cluster_indices[j]])
				cluster_traj = trajectory_file[xtc_frame]
				cluster_traj.save_pdb('%s/cluster_centers/state_%d_%.3f.pdb'%(self.tICA_dir, cluster_labels[j],percentages[cluster_labels[j]]))
				if verbose:
					print('Successfully saved PDB for cluster: %d, (rel.pop: %.3f)'%(cluster_labels[j],percentages[cluster_labels[j]]))
					print('traj_file: %s (%d/%d)'%(trajectory_file,i,len(features)))
					print('frame: %d (%d/%d centers from this trajectory)'%(cluster_indices[j],j,len(cluster_indices)))
					print('strided: npy_frame/npy_len = %d/%d = %f'%(cluster_indices[j],n_snapshots,cluster_indices[j]/n_snapshots))
					print('re-mapped: orig_frame/xtc_len = %d/%d = %f\n'%(xtc_frame,xtc_len,xtc_frame/xtc_len))
	#if calculate_strides:    
		#os.remove('/home/server/git/fah-scripts/DataAnalysisScripts/temp_dist/stride_dist_%d.npy' % self.proj_num)



