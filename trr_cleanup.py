import os, re, sys, glob, time, subprocess, matplotlib, socket
import numpy as np
from matplotlib import pyplot as plt
from tqdm import trange

# debug purposes on laptop
if socket.gethostname() == 'syzygy':
    debug_prefix = '/media/matt/ext'
else:
    debug_prefix = ''


script_name = sys.argv[0]

usage = f"""Description:

    This script will delete all unnecessasry *.trr files for all projects or a set number of projects.
    CAUTION: This is dangerous, and will irreparably destroy trajectory data -- ask yourself: are you
    sure EVERYTHING you need to analyze the data is written to the *.xtc files?

NOTE:

    This should only be used for a GROMACS project, with trajectory data in a nested structure like this:
    /array1/server2/data/SVR166220/PROJ14627/RUN0/CLONE0/results9/frame9.trr

Usage:

    $ python {script_name} [test|delete] PROJDIR  

    To perform a "dry-run" that will assess how much data is to be deleted, use "test"
    Example: python {script_name} test /array1/server2/data/SVR166220/PROJ14627

    To delete the data, use "delete"
    Example: python {script_name} delete /array1/server2/data/SVR166220/PROJ14627

"""

if len(sys.argv) < 3 or sys.argv[1] not in ['delete', 'test']:
    print(usage)
    sys.exit(1)

script_name = sys.argv[0]
action      = sys.argv[1]
projdir     = os.path.abspath( sys.argv[2] )

if sys.argv[1] == 'delete':
    delete = True
if sys.argv[1] == 'test':
    delete = False


#### change these ####
verbose = False
count, sizes = 0, []

projects = [ projdir ] 

for i in projects:

    project_data_path = i
          
    # How many RUNs are there?
    runs = int( subprocess.check_output('ls -1 %s | grep RUN | wc -l'%project_data_path, shell=True).split()[0].decode('utf-8') )
    print('runs =', runs)

    # How many CLONEs are there?
    run0_path = os.path.join(project_data_path, 'RUN0')
    clones = int( subprocess.check_output('ls -1 %s | grep CLONE | wc -l'%run0_path, shell=True).split()[0].decode('utf-8') )
    print('clones =', clones)

    print("Processing project: ", i)
    for j in trange(runs, desc='Runs',leave=True,smoothing=1):
        for k in range(clones):

                trr_wildcard = os.path.join(project_data_path, 'RUN%d/CLONE%d/results*/frame*.trr'%(j, k) )
                trr_files = glob.glob(trr_wildcard)
                print('Removing %d *.trr files)'%(len(trr_files)), '...')


                # find the maximum frame index
                basenames = [os.path.basename(s) for s in trr_files]
                frame_indices = [int(s.replace('frame','').replace('.trr','')) for s in basenames]
                Ind = np.argsort(frame_indices)
                #print(Ind)

                sorted_trr_files  = [trr_files[m] for m in Ind]
                #print(sorted_trr_files)
                
                ### delete all but the last TWO *.trr files
                for trr_file in sorted_trr_files[0:-2]:

                    size_cmd = 'du -sh %s' %trr_file
                    size = subprocess.check_output('du -sh ' + trr_file, shell=True).split()[0].decode('utf-8')
                    if 'G' in size:
                        size = float(size[:-1])
                    elif 'M' in size:
                        size = float(size[:-1])/1024
                    elif 'K' in size:
                        size = float(size[:-1])/1024/1024
                    sizes.append(size)
                        
                    if verbose:
                        print("Removing: " + trr_file)
                        
                    if delete:
                        os.remove(trr_file)
                        
                    count += 1

                for trr_file in sorted_trr_files[-2:]:
                    print("\tKeeping: " + trr_file)

                 

if delete:
    print("\nA total of %s files were deleted, freeing %s GB" %(str(count),str(sum(sizes))))
else:
    print("\nA total of %s files could be deleted, freeing %s GB" %(str(count),str(sum(sizes))))
