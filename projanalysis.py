## CODE FOR CALCULATING DISTRIBUTION OF TRAJECTORY LENGTHS FOR ALL PROJECTS PER DAY

import os, re, sys, glob, time, subprocess, matplotlib
matplotlib.use('Agg')
# %matplotlib inline
from matplotlib import pyplot as plt

# server specifics, set once and don't change
projects, project_numbers, project_data_dirs, count = [],[],[], 0

# User and Hostname definitions
hostname = subprocess.check_output('hostname', shell=True).split('.')[0]

# Path definitions
debug_prefix = ''
if hostname == 'vav3':
    array0 = debug_prefix + '/array0/data/'
    array1 = debug_prefix + '/array1/server2/data/SVR166219/'
elif hostname == 'vav4':
    array0 = debug_prefix + '/array0/projectdata/'
    array1 = debug_prefix + '/array1/server2/data/SVR166220/'
else:
    print('Hostname not recognized. Define hostname and paths in slackbot.py.')


project_config_prefix = '/array1/server2/projects/Gromacs/'
results_prefix = '/array1/server2/.results/'
timestr = time.strftime("%Y%m%d")
logfile = open(results_prefix + "log.txt", "a+")
sizefile = results_prefix + "sizes.txt"
find_missing_projects = True

# choose to update trajectory length distribution or project size, one must be True
TrajLenDistribution = True
ProjSize = True

# use this to update a particular (set of) project(s)
custom_project_numbers = [] # e.g. [9001] ['13337','6969'] range(1000,10000)

# make sure a calculation is selected
if not TrajLenDistribution and not ProjSize:
    print('Must choose either TrajLenDistribution, ProjSize, or both.')
    sys.exit(0)

# checks for projects in data directories that do not yet have results processed
if find_missing_projects and custom_project_numbers == []:
    for i in sorted(glob.glob(results_prefix + "*/*_log.txt")):
        project_file = re.sub(results_prefix + ".*/", '', i)
        project_number = re.sub("_log.txt", '', project_file)
        project_numbers.append(project_number)
        project_numbers = sorted(project_numbers)
    for array in [array1,array0]:
        for i in sorted(glob.glob(array + "PROJ*")):
            projects.append(re.sub(r'^.*PROJ', "",i))

for i in projects:
    if i not in project_numbers:
        custom_project_numbers.append(i)


# runs parsing script for projects, will default to those without results procesed
# update this to determine array first
for array in [array1, array0]:

    if custom_project_numbers != []:
        print("Attempting to parse %s for missing projects:"%array,custom_project_numbers)
        projects = [array + "PROJ%s"%str(i) for i in custom_project_numbers]
    else:
        print("No new projects. Updating all project summaries.")
        projects = sorted(glob.glob(array + "PROJ*"))

    for i in projects:
        if not os.path.exists(i): # not on that array
            continue

        print("Processing project: ", i)
        project_data_path = i
        project_number = re.sub(r'^.*PROJ', "", i)
        project_config_path = project_config_prefix + 'p' + project_number
        project_xml = project_config_path + '/project.xml'
        try:
            project_size = subprocess.check_output('du -sh ' + project_data_path, shell=True).split()[0].decode('utf-8')
        except Exception as e:
            print('g',e)
            break
    
        try: # get the mdp path
            mdp_path = ''
            for file in os.listdir(project_config_path):
                if 'mdout' not in file and 'mdp' in file and 'mdp' != file: # we need to standardize mdp loc/names
                    mdp_path = project_config_path + '/' + file
            if mdp_path == '':
                for file in os.listdir(project_config_path + '/RUN0'):
                    if 'mdout' not in file and 'mdp' in file and 'mdp' != file:
                        mdp_path = project_config_path + '/RUN0/' + file
        except Exception as e:
            print('error a', ': empty directory')
            logfile.write('\n' + timestr + "-- No config directory for project: " + project_number)
            continue
    
        try: # parse mdp for frame length
            with open(mdp_path, 'r') as file:
                lines = [line.strip().split() for line in file.readlines()]
            file.close()
        
            for j in range(len(lines)):
                if 'nsteps' in lines[j]:
                    ns_per_frame = float(lines[j][2])/500000.
        except Exception as e:
            print('error b',mdp_path,e)
            logfile.write('\n' + timestr + "-- No mdp file found for project: " + project_number)
            continue
        
        try: # parse xml for stuff 
            with open(project_xml, 'r') as file:
                lines = [line.strip().split() for line in file.readlines()]
            file.close()

            for j in range(len(lines)):
                if '<contact' in lines[j]:
                    user = re.sub(r'@.*', "", lines[j][1][3:])
                if '<runs' in lines[j]:
                    runs = int(lines[j][1][3:-3])
                if '<clones' in lines[j]:
                    clones = int(lines[j][1][3:-3])
                if '<gens' in lines[j]:
                    gens = int(lines[j][1][3:-3])

        except Exception as e:
            print('error c',e)
            logfile.write('\n' + timestr + "-- No xml file found for project: " + project_number)
            continue
        
        user_results_dir = results_prefix + '/' + user
        if not os.path.exists(user_results_dir):
            os.makedirs(user_results_dir)
        project_log=open(user_results_dir + '/' + project_number + "_log.txt", "a+")

# HONGBIN CODE BELOW for getting frame length distribution

        if TrajLenDistribution:
            Continue = True
            lengths, traj_len_distribution_x, traj_len_distribution_y = [],[],[]

            frame = 0
            maxframe = gens-1
            total_ns = 0.
        
            while Continue:
            
                try:
                    cmd = 'ls -1 %s | wc -l'%(project_data_path + '/RUN*/CLONE*/results%d/traj_comp.xtc'%frame)
                    n = int(subprocess.check_output(cmd, shell=True, stderr=open(os.devnull, 'w')))
                except ValueError as e:
                    print(e)
                    n = 0
                
                total_ns += ns_per_frame*n
                if (n < 1) or frame >= maxframe:
                    Continue = False
                else:
                    frame +=1

                traj_len_distribution_x.append((frame+1)*ns_per_frame)
                traj_len_distribution_y.append(n)

            bins = int(ns_per_frame*(frame+1)/2)

            plt.figure(count)
            plt.bar(traj_len_distribution_x, traj_len_distribution_y, width=ns_per_frame)
            plt.xlabel('Clone Length (ns)')
            plt.ylabel('Number of Trajectories')
            plt.title('Trajectory length distribution for project: ' + project_number)
            plt.savefig(user_results_dir + '/' + project_number +'.png')
            plt.close()
           
            count += 1

            project_log.write('\n' + timestr + ' ' + str(total_ns) + ' ' + project_size)

            project_log.close()

# log total frames/size

        else:
            project_log.write('\n' + timestr + ' ' + 'N/A' + ' ' + project_size)

sizes, lengths = [],[]
data_log = results_prefix + 'all_data.txt'

if not os.path.exists(data_log):
    f = open(data_log, 'w')
    f.write('YYYYMMDD, total_size, total_length, change_size, change_len\n')
    f.close()

for i in sorted(glob.glob(results_prefix + '*/*_log.txt')):

    with open(i, 'r') as file:
        lines = [line.strip().split() for line in file.readlines()]
    file.close()
    
    if lines == []:
        break
    
    size = lines[-1][2]
    if 'T' in size:
        size = float(size[:-1])*1024
    elif 'G' in size:
        size = float(size[:-1])
    elif 'M' in size:
        size = float(size[:-1])/1024
    else:
        size = 0

    sizes.append(size)

    if TrajLenDistribution:
        lengths.append(float(lines[-1][1]))

with open(data_log) as file:
    lines = [line.strip().split() for line in file.readlines()]
file.close()

try:
    del_size = sum(sizes) - float(lines[-1][1])
    if TrajLenDistribution:
        del_length = sum(lengths) - float(lines[-1][2])
    else:
        del_length = 0

except Exception as e:
    del_size, del_length = 0, 0


with open(data_log, 'a+') as file:
    file.write(timestr + " " + str(sum(sizes))+ " " + str(sum(lengths)) + " " + str(del_size) + " " + str(del_length))
file.close()

# plotting stuff

y_delsize, y_dellength = [],[]

with open(data_log) as file:
    lines = [line.strip().split() for line in file.readlines()][1:]
file.close()

for i in range(len(lines)):
    y_delsize.append(lines[i][3])
    y_dellength.append(lines[i][4])

plt.figure(1337)
plt.plot(range(len(y_delsize)), y_delsize)
plt.savefig(results_prefix + "overall_delsize.png")
plt.close()
plt.figure(1338)
plt.plot(range(len(y_dellength)), y_dellength)
plt.savefig(results_prefix + "overall_dellength")
plt.close()

logfile.close()

