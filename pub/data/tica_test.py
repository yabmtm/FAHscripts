#!/home/server/anaconda2/bin/python

from __future__ import print_function
import glob, itertools, os, subprocess, re
import sys, time, tqdm, itertools, random
import mdtraj as md
import numpy as np
from itertools import groupby, count
import matplotlib, msmbuilder
%matplotlib inline
# Script: matplotlib.use('Agg')  | Notebook: %matplotlib inline
from matplotlib import pyplot as plt
from msmbuilder.cluster import KCenters, KMeans, KMedoids
from msmbuilder.decomposition import tICA
from msmbuilder.featurizer import AtomPairsFeaturizer
from msmbuilder.msm import MarkovStateModel
from msmbuilder.msm.core import *
from pyemma.thermo.util import get_averaged_bias_matrix as _get_averaged_bias_matrix
import pyemma

'''This script has lots of functionality and is based on analyzing Gromacs trajectories. A list of trajectory
   files is given as trajectory_files, as well as a general structure file. Other structure files should contain
   the same name as the corresponding trajectory file, e.g. traj_001.trr traj_001.gro.
'''

# this project represents a spiroligomer (1) from https://doi.org/10.1371/journal.pone.0045948
# bound to MDM2 (PDB: 1ycr)
# these runs represent the 20 ensembles: lam = false, d = {0.1, 2.0, 0.1}

# featurization parameters
project_title = 'PROJ14101' # creates sub-directory
structure_file = '%s/xtc.gro'%project_title
runs = len(glob.glob('%s/traj_data/RUN*'%project_title))
clones = max([len(glob.glob('%s/traj_data/RUN%d/*'%(project_title,x))) for x in range(runs)])
md_time_step = 0.02 # time in ns that trajectory files are saved (nstxtcout)
equil_time = 1. # time in ns to remove from beginning of each clone
subsampled_time_step = 0.5 # preferred subsampled time-step in ns
stride = int(subsampled_time_step / md_time_step) # time step stride for sub-sampling
equil_steps = int(equil_time / md_time_step) # time steps to be removed from start
custom_residues = ['B1A','B1B','B2A','B2B','B2C','B2D','B2E','B3A','B3B']
custom_residues += ['B4A','B4B','B4C','B4D','B5A','B5B','B5C','B6A'] # for spiroligomers
custom_residues += ['1MQ','20Q','20U','I18','I31','K23','NUT','YIN'] # for nutlins

# tICA/MSM parameters
unbiased_runs = [20]
tica_lagtime = 10 # determine from implied timescales / GMRQ
n_components = 8 # how many tICs to compute
n_clusters = 50 # denotes number of microstates
n_timescales = n_components # plot all eigenvalues --> timescales
lagtimes = np.array([1,2,4,8,16,32,64,128,256,512,1024]) # log scale
cluster_method = 'kcenters' # 'kcenters/kmeans/kmedoids'
all_ticas = list(itertools.permutations(range(1,n_components+1), 2)) # all combinations
all_ticas = [[1,2],[1,3],[2,3]] # override: just show analysis for first three components
cluster_percentage_cutoff = n_clusters/64 # clusters with a relative population less than this
                              # number will not be labeled on plot i.e. 0 : all clusters labeled


def compute_tica_components():
          
    '''Load in the features, calculate a given number of tICA components (tica_components) given a
       lagtime (lag_time), and save tICA coordinates and eigenvector data. It then creates and populates
       a list for each desired component, clusters the data, saving normalized populations as populations.dat
       and saving each cluster center as a .pdb. tICA plots are created and saved, and implied timescales are
       calculated, saved, and plotted.
    '''
        
    verbose = False
    save_pdb = True
    color_by = 'cluster'
    
    if verbose:
        print("\nCalculating tICA components...")
    if not os.path.exists(project_title + '/tica_%d'%n_clusters):
        os.mkdir(project_title + '/tica_%d'%n_clusters)
    
    # load in feature files and determine indices of unbiased ensembles
    feature_files = []
    for i in range(runs):
        run_files = sorted(glob.glob(/features/' + "P*R%d_*npy"%i))
        feature_files += run_files
        if i in unbiased_runs:
            unbiased_indices = [len(feature_files) - len(run_files),len(feature_files)]
    features = [np.load(x) for x in feature_files]
    
    # perform tICA calculation and extract score / eigenvectors
    tica_coordinates = tICA(lag_time=tica_lagtime,
        n_components=int(n_components)).fit_transform(features)
    tica_components = tICA(lag_time=tica_lagtime,
        n_components=int(n_components)).fit(features)
    eigenvectors = np.transpose(tica_components.eigenvectors_)
    tica_score = tica_components.score(features)
          
    np.save('%s/tica_%d/tica_coords-lag_%d-comp_%d.npy' %(
        project_title, n_clusters, tica_lagtime, n_components), tica_coordinates)
    np.save('%s/tica_%d/tica_comps-lag_%d-comp_%d.npy' %(
        project_title, n_clusters, tica_lagtime, n_components), tica_components)
    
    # Perform clustering based on the cluster_method parameter.
    if verbose:
        print('Clustering via %s'%cluster_method)
    if cluster_method == 'kcenters':
        clusters = KCenters(n_clusters)
    elif cluster_method == 'kmeans':
        clusters = KMeans(n_clusters)
    elif cluster_method == 'kmedoids':
        clusters = KMedoids(n_clusters)
    else:
        sys.exit('Invalid cluster_method. Use kcenters/kmeans/kmedoids.')
        
    # Cluster unbiased data and fit biased data to these centers
    new_assignments = []
    sequences = clusters.fit_transform(tica_coordinates[unbiased_indices[0]:unbiased_indices[1]])
    for i in tqdm.tqdm_notebook(range(unbiased_indices[0])):
        tica_traj = tica_coordinates[i]
        if isinstance(tica_traj, np.ndarray):
            if not (tica_traj.dtype == 'float32' or tica_traj.dtype == 'float64'):
                tica_traj = tica_traj.astype('float64')
        labels, inertia = msmbuilder.libdistance.assign_nearest(
            tica_traj, clusters.cluster_centers_, metric='euclidean')
        new_assignments.append(labels)

    new_assignments += sequences # tack the unbiased assignments back on to the end.


    np.save('%s/tica_%d/lag_%d_clusters_%d_assignments.npy' %(
        project_title, n_clusters, tica_lagtime, n_clusters), new_assignments)
    np.save('%s/tica_%d/lag_%d_clusters_%d_center.npy' %(
        project_title, n_clusters, tica_lagtime, n_clusters), clusters.cluster_centers_)

    # Determine cluster populations, normalize the counts, and save as percentages for
    # labeling if a cluster contains more than cluster_percentage_cutoff percent of the data.
    # Finally, save normalized counts.
    
    if verbose:
        print("\nDetermining cluster populations...")
    if not os.path.exists('%s/tica_%d/%s_clusters'%(project_title,n_clusters,cluster_method)):
        os.mkdir('%s/tica_%d/%s_clusters'%(project_title,n_clusters,cluster_method))
    if not os.path.exists('%s/tica_%d/plots'%(project_title,n_clusters)):
        os.mkdir('%s/tica_%d/plots'%(project_title,n_clusters))
        
    counts = np.array([len(np.where(np.concatenate(sequences)==i)[0]) for i in range(n_clusters)])
    normalized_counts =  counts/float(counts.sum())
    percentages = [ i*100 for i in normalized_counts ]
    population_labels = [ [i,"%.2f"%percentages[i]] for i in range(len(percentages)) if percentages[i] > cluster_percentage_cutoff ]
    np.savetxt('%s/tica_%d/%s_clusters/populations.dat'
               %(project_title,n_clusters,cluster_method), normalized_counts)

    # Plot all unique combinations of tICA components
    if verbose:
        print("\nPlotting tICA components...")
    tica_coordinates = np.concatenate(tica_coordinates)
    new_assignments = np.concatenate(new_assignments)
    cluster_colors = matplotlib.cm.rainbow(np.linspace(0,1,n_clusters))
    for j in tqdm.tqdm_notebook(range(len(all_ticas)),leave=False): # For each pair
        if all_ticas[j][0] < all_ticas[j][1]:
            plt.figure(j, figsize=(20,16))
            tICx, tICy = all_ticas[j][0]-1, all_ticas[j][1]-1
            plt.hexbin(tica_coordinates[:,tICx],tica_coordinates[:,tICy], bins='log')
            for l in tqdm.tqdm(range(len(tica_coordinates))[::stride*2]):
                if color_by == 'cluster':
                    plt.plot(tica_coordinates[l][tICx], tica_coordinates[l][tICy],
                        color=cluster_colors[new_assignments[l]], linestyle="", marker="o")
            x_centers = [clusters.cluster_centers_[i][tICx] for i in range(len(clusters.cluster_centers_))]
            y_centers = [clusters.cluster_centers_[i][tICy] for i in range(len(clusters.cluster_centers_))]
            high_pop_x_centers = [ x_centers[i] for i in range(len(x_centers)) if percentages[i] > cluster_percentage_cutoff ]
            high_pop_y_centers = [ y_centers[i] for i in range(len(y_centers)) if percentages[i] > cluster_percentage_cutoff ]
            plt.plot(x_centers, y_centers, color='y', linestyle="", marker="o")
            plt.plot(tica_coordinates[:,tICx][0],tica_coordinates[:,tICy][0], color='k', marker='*',markersize=24)
            plt.xlabel('tIC'+str(all_ticas[j][0]))
            plt.ylabel('tIC'+str(all_ticas[j][1]))
            plt.title(project_title)
            # Add labels for high-population cluster centers
            for label, x, y in zip(population_labels, high_pop_x_centers, high_pop_y_centers):
                plt.annotate(
                  label,
                  xy = (x, y), xytext = (-15, 15),
                  textcoords = 'offset points', ha = 'right', va = 'bottom',
                  bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
                  arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
            plt.savefig('%s/tica_%d/plots/tica_%d_%d.png'%(project_title,n_clusters,
                all_ticas[j][0], all_ticas[j][1]))
            plt.close()

    # Write out PDBs for each cluster center
    if verbose:
        print("Performing cluster analytics and saving center PDBs...\n")
    if save_pdb:
        trajectory_files, feature_files, cluster_features = [],[],[]
        for run in range(runs): # get only xtc files that correlate to cluster-center features
            trajectory_files += [re.sub('features',
                                    'traj_data/RUN%d'%run,re.sub('npy','xtc',x)
                                     ) for x in sorted(glob.glob('%s/features/*R%d_*npy'%(
                                        project_title,run)))]
            feature_files += sorted(glob.glob('%s/features/*R%d_*npy'%(project_title,run)))

        for i in tqdm.tqdm_notebook(range(len(trajectory_files)),leave=False):

                n_snapshots = len(clusters.distances_[i])

                # Determine frames that are cluster centers
                cluster_indices = np.arange(n_snapshots)[ (clusters.distances_[i] < 1e-6) ]

                # Determine number of each cluster, correlates to populations.dat
                cluster_labels = sequences[i][cluster_indices]

                # Save each cluster center as a pdb
                if list(cluster_indices): # load center-containing xtcs to check length
                    xtc_len = len(md.load(trajectory_files[i],top=structure_file))
                    
                # map strided frame number back to xtc frame number
                for j in range(len(cluster_indices)):
                        frames = range(xtc_len) 
                        strided_frames = frames[equil_steps:][::stride]
                        xtc_frame = frames.index(strided_frames[cluster_indices[j]])
                        cluster_traj = md.load_frame(trajectory_files[i], xtc_frame,
                                            top=structure_file)
                        cluster_features.append(np.load(feature_files[i])[cluster_indices[j]])
                        cluster_traj.save_pdb('%s/tica_%d/%s_clusters/state_%d.pdb'
                                            %(project_title,n_clusters,cluster_method,
                                            cluster_labels[j]))
                        
                        # save cluster information
                        with open('%s/tica_%d/cluster.dat'%(project_title,n_clusters),'w') as f:
                            f.write('\nSuccessfully saved PDB for cluster: %d, (rel.pop: %.3f)'%(
                                cluster_labels[j],percentages[cluster_labels[j]]))
                            f.write('traj_file: %s (%d/%d)'%(trajectory_files[i],i,len(features)))
                            f.write('frame: %d (%d/%d centers from this trajectory)'%(
                                cluster_indices[j],j,len(cluster_indices)))
                            f.write('strided: npy_frame/npy_len = %d/%d = %f'%(
                                cluster_indices[j],n_snapshots,cluster_indices[j]/n_snapshots))
                            f.write('re-mapped: orig_frame/xtc_len = %d/%d = %f\n'%(
                                xtc_frame,xtc_len,xtc_frame/xtc_len))
                            f.close()
                        
        # save features corresponding to each cluster center
        np.save('%s/tica_%d/cluster_features.npy'%(project_title,n_clusters),cluster_features)
                    
    return tica_score
