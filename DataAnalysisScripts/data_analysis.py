from analysis_featurizer import *

base_dir = '/home/server/server2/data/SVR166219/'
gro_dir = '/home/server/server2/projects/Gromacs/'

def get_shape(proj_num):
    """ Given a project number, gets the numbers of runs, clones, gens of the project and returns them as a list """

    project_xml = open(os.path.join(gro_dir, 'p' + str(proj_num), 'project.xml'), 'r')
    xml = project_xml.read()
    soup = BeautifulSoup(xml, 'lxml')
    runs, clones, gens = int(soup.find('runs')['v']), int(soup.find('clones')['v']), int(soup.find('gens')['v'])
    shape_list = [runs, clones, gens]
    return shape_list

def distance_calc_proj(proj_num):
    """ Given a project number, calculates all distances for a whole project """

    shape = get_shape(proj_num)
    t = DistanceCalculator(proj_num, 0, 0, 0)
    dists = np.load(t.data_path)
    for i in range(shape[0]):
	for j in range(shape[1]):
	    for k in range(shape[2]):
		if os.path.isfile(os.path.join(base_dir, 'PROJ' + str(proj_num), 
					       'RUN' + str(i), 'CLONE' + str(j), 
					       'results' + str(k), 'traj_comp.xtc')): 
		    if dists[i][j][k] is None or dists[i][j][k] is 0:
			t = DistanceCalculator(proj_num, i, j, k)
		        t.calculate_distances()
			del t
		    else:
		        pass
		    
    print('Project finished')

    
def dihedral_calc_proj(proj_num, angle='omega'):
    """ Given a project number and an angle (default 'omega'), calculates dihedrals for a whole project """    

    shape = get_shape(proj_num)
    t = DihedralCalculator(proj_num, 0, 0, 0, angle)
    dih_npy = np.load(t.data_path)
    for i in range(shape[0]):
        for j in range(shape[1]):
            for k in range(shape[2]):
		if os.path.isfile(os.path.join(base_dir, 'PROJ' + str(proj_num), 
                                               'RUN' + str(i), 'CLONE' + str(j),
					       'results' + str(k), 'traj_comp.xtc')):
		    if dih_npy[i][j][k] is None or dih_npy[i][j][k] is 0:
			t = DihedralCalculator(proj_num, i, j, k, angle)
			t.calculate_dihedrals()
			del t
     
		   
	
def COM_calc_proj(proj_num):
    """ Given a project number, calculates COM for the whole project """
    shape = get_shape(proj_num)
    t = COMCalculator(proj_num, 0, 0, 0)
    COM = np.load(t.data_path)
    for i in range(10, 20):
        for j in range(shape[1]):
            for k in range(shape[2]):
		if os.path.isfile(os.path.join(base_dir, 'PROJ' + str(proj_num), 
                                               'RUN' + str(i), 'CLONE' + str(j), 
                                               'results' + str(k), 'traj_comp.xtc')):
		    if COM[i][j][k] is None or COM[i][j][k] is 0:
			t = COMCalculator(proj_num, i, j, k)
			t.COM_atom_distances_for_gen()
			del t

def RMSD_calc_proj(proj_num, ref_struct=None):
    """ Given a project number and a reference structure (default None: Reference is the first frame),
	calculates RMSD for an entire project """
    shape = get_shape(proj_num)
    t = RMSDCalculator(proj_num, 0, 0, 0)
    RMSD = np.load(t.data_path)
    for i in range(shape[0]):
        for j in range(shape[1]):
            for k in range(shape[2]):
                if os.path.isfile(os.path.join(base_dir, 'PROJ' + str(proj_num),
                                               'RUN' + str(i), 'CLONE' + str(j),
                                               'results' + str(k), 'traj_comp.xtc')):
		    if RMSD[i][j][k] is None or RMSD[i][j][k] is 0:
			t = RMSDCalculator(proj_num, i, j, k)
			t.calculate_RMSD(ref_struct)
			del t


def COM_plot_proj(proj_num):
    """ Given a project number, plots all available COM data for the entire project """
    shape = get_shape(proj_num)
    t = COMCalculator(proj_num, 0, 0, 0)
    COM = np.load(t.data_path)
    for i in range(shape[0]):
	for j in range(shape[1]):
	    if COM[i][j][0] is not None:
	        t = COMCalculator(proj_num, i, j, 0)
	        t.plot_COM_atom_distances()
		del t



