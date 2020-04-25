import os, re, sys, glob, time, subprocess, matplotlib, socket
from matplotlib import pyplot as plt
from tqdm import trange

# debug purposes on laptop
if socket.gethostname() == 'syzygy':
    debug_prefix = '/media/matt/ext'
else:
    debug_prefix = ''

if len(sys.argv) != 2 or sys.argv[1] not in ['delete', 'test']:
    print('This script will delete all unnecessasry tpr files for all projects or a set number of projects.')
    print('Usage:')
    print('Delete: python %s delete'%sys.argv[0])
    print('Dry-run: python %s test'%sys.argv[0]) 
    sys.exit(0)

if sys.argv[1] == 'delete':
    delete = True
if sys.argv[1] == 'test':
    delete = False

#### CHANGE THESE ####
delete = delete # (True/False)
custom_project_numbers = [] #range(14051,14055) # [14100, 14101, 14102]
######################

verbose = False
count, sizes, project_data_dirs = 0,[],[]

# Hostname definitions
hostname = (subprocess.check_output('hostname', shell=True).decode()).split('.')[0]   # need a decode here
print('hostname', hostname)

# Path definitions
if hostname == 'vav3':
    array0 = debug_prefix + '/array0/data/'
    array1 = debug_prefix + '/array1/server2/data/SVR166219/'
    arrays = [array1, array0]
    results_dir = debug_prefix + '/array1/server2/.results'
elif hostname == 'vav4':
    array0 = debug_prefix + '/array0/projectdata/'
    array1 = debug_prefix + '/array1/server2/data/SVR166220/'
    arrays = [array1, array0]
    results_dir = debug_prefix + '/array1/server2/.results'
elif hostname == 'vav15':
    array1 = debug_prefix + '/data/SVR2616698069'
    arrays = [array1]
    results_dir = debug_prefix + '/home/server/server2/.results'
elif hostname == 'vav16':
    array1 = debug_prefix + '/data/SVR2616698070'
    arrays = [array1]
    results_dir = debug_prefix + '/home/server/server2/.results'
else:
    print('Hostname not recognized. Define hostname and paths in slackbot.py.')


project_config_prefix = debug_prefix + '/array1/server2/projects/Gromacs/'
timestr = time.strftime("%Y%m%d")
logfile=open( os.path.join(results_dir,"log.txt"), "a+")

for array in arrays:
    if custom_project_numbers != []:
        projects = [array + "PROJ%s"%str(i) for i in custom_project_numbers]
    else:
        projects = sorted(glob.glob(array + "PROJ*"))

    for i in projects:
        project_data_path = i
        project_number = re.sub(r'^.*PROJ', "", i)
        project_config_path = project_config_prefix + 'p' + project_number
        project_xml = project_config_path + '/project.xml'
        try: # parse xml for stuff 
            
            with open(project_xml, 'r') as file:
                lines = [line.strip().split() for line in file.readlines()]
            file.close()

            for j in range(len(lines)):
                if '<runs' in lines[j]:
                    runs = int(lines[j][1][3:-3])
                if '<clones' in lines[j]:
                    clones = int(lines[j][1][3:-3])
                if '<gens' in lines[j]:
                    gens = int(lines[j][1][3:-3])

            print("Processing project: ", i)
            for j in trange(runs, desc='Runs',leave=True,smoothing=1):
                for k in range(clones):
                    for l in range(gens+1):
                        tpr_file = project_data_path + '/RUN' + str(j) + '/CLONE' + str(k) + '/frame%s.tpr' %str(l)
                        next_tpr_file = project_data_path + '/RUN' + str(j) + '/CLONE' + str(k) + '/frame%s.tpr' %str(l+1)
                        size_cmd = 'du -sh %s' %tpr_file
                        if l == 0:
                            if verbose:
                                print("Not removing: " + tpr_file)
                        elif os.path.isfile(tpr_file) and os.path.isfile(next_tpr_file):
                            
                            size = subprocess.check_output('du -sh ' + tpr_file, shell=True).split()[0].decode('utf-8')
                            if 'G' in size:
                                size = float(size[:-1])
                            elif 'M' in size:
                                size = float(size[:-1])/1024
                            elif 'K' in size:
                                size = float(size[:-1])/1024/1024
                            sizes.append(size)
                                
                            if verbose:
                                print("Removing: " + tpr_file)
                                
                            if delete:
                                os.remove(tpr_file)
                                
                            count += 1

                        elif os.path.isfile(tpr_file) and not os.path.isfile(next_tpr_file):
                            if verbose:
                                print("Not removing: " + tpr_file)
                        else:
                            pass
            
        except IOError:
            print("\nproject.xml not found for project: " + project_number)
            pass
if delete:
    print("\nA total of %s files were deleted, freeing %s GB" %(str(count),str(sum(sizes))))
else:
    print("\nA total of %s files could be deleted, freeing %s GB" %(str(count),str(sum(sizes))))
