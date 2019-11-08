#!/usr/bin/env python

""" This script should be run in the project directory, e.g. p1337. It will parse the contents of the directory
    for all the necessary simulation files and parameters. Once these are extracted, it can build a project.xml
    to be used for the FaH job. It will also create a tpr_test directory and build a test tpr, which can be run
    with "gmx mdrun -nt X -deffnm tpr_test.tpr" where X is the number of cores to use; Default vav3/gmx5: 24.
"""

from __future__ import print_function
import os, re, subprocess, sys

# Check that current directory is of type ~/server2/projects/Gromacs/pXXXXX
try:
    cmd = 'echo $(pwd) | sed "s/^.*Gromacs\/p//" | sed "s/\/[A-Z].*[0-9].*//"'
    ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    proj_id = int(ps.communicate()[0])
except ValueError, e:   
    print("Cannot parse project ID from current directory.")

# Currently only Gromacs project directories are supported. Raise an issue for this to change.
project_type = "gromacs"

# Prints contents of project directory
print("\nContents of current directory: (" + os.getcwd() + ")")
for i in os.listdir(os.getcwd()):
    print(i, end="  ")
print("\n",proj_id)

# Parse directory for coordinate file. Prompt if more than one exists.
structure=[]
for file in os.listdir(os.getcwd()):
    if file.endswith(".gro"):
        structure.append(file)
while len(structure) > 1:
    input = raw_input("There is more than one structure file in this directory. Please input structure file: ")
    if input in os.listdir(os.getcwd()):
        structure=[input]
    else:
        print("That file is not in this directory.")
if len(structure) == 0:
    print("There is no structure file in this directory."); sys.exit(0)
structure = structure[0]

# Parse directory for topology file. Prompt if more than one exists.
topology=[]
for file in os.listdir("."):
    if file.endswith(".top"):
        topology.append(file)
while len(topology) > 1:
    input = raw_input("There is more than one topology file in this directory. Please input topology file: ")
    if input in os.listdir(os.getcwd()):
        topology=[input]
    else:
        print("That file is not in this directory.")
if len(topology) == 0:
    print("There is no topology file in this directory."); sys.exit(0)
topology = topology[0]

# Parse directory for mdp file. Prompt if more than one exists.
mdp=[]
for file in os.listdir("."):
    if file.endswith(".mdp") and file != "mdout.mdp" and file != "mdpOut.mdp":
        mdp.append(file)
while len(mdp) > 1:
    input = raw_input("There is more than one mdp file in this directory. Please input mdp file: ")
    if input in os.listdir(os.getcwd()):
        mdp=[input]
    else:
        print("That file is not in this directory.")
if len(mdp) == 0:
    print("There is no mdp file in this directory."); sys.exit(0)
mdp = mdp[0]

# Parse directory for index file. Prompt if more than one exists.
index_=[]
for file in os.listdir("."):
    if file.endswith(".ndx"):
        index_.append(file)
while len(index_) > 1:
    input = raw_input("There is more than one index file in this directory. Please input index file: ")
    if input in os.listdir(os.getcwd()):
        index_=[input]
    else:
        print("That file is not in this directory.")
if len(index_) == 0:
    print("There is no index file in this directory. Creating a new one...")
    os.system("echo q | gmx make_ndx -f " + structure)
    index_ = "index.ndx"
    print("Index file was missing. A new one was created based on" + structure + ".")
index_ = index_[0]

try: # extracting number of steps from mdp file
    cmd = 'cat ' + mdp + ' | grep nsteps | grep -o "[0-9].*" | sed "s/ *;//" | cut -d" " -f1'
    ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    steps = int(ps.communicate()[0])
except ValueError, e:
    print("Cannot parse mdp file for number of steps.")

try: # extracting timestep from mdp file
    cmd = 'cat ' + mdp + ' | grep dt | grep -o "[0-9].*" | sed "s/ *//" | cut -d" " -f1'
    ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    timestep = re.findall("\d+\.\d+", ps.communicate()[0])
    timestep = float(timestep[0])
except ValueError, e:
    print("Cannot parse mdp file for timestep.")

try: # extracting number of atoms from structure file
    cmd = 'head -n 2 ' + structure + ' | sed "1d"'
    ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    atoms = int(ps.communicate()[0])
except ValueError, e:
    print("Structure file does not have number of atoms on second line.")
    
try: # calculating time per frame, points, timeout, and deadline
    simulation_time = steps * timestep / 1000
    points_per_ns = 0.006285*atoms+6.588
    points = int(points_per_ns * simulation_time)
    timeout_per_ns = 0.00010024*atoms-0.034292
    timeout = float(timeout_per_ns * simulation_time)
    deadline_per_ns = 0.00021931*atoms-0.065796 
    deadline = float(deadline_per_ns * simulation_time)

    
except ValueError, e:
    print("Something... went wrong.")
    sys.exit(0)

# Print project summary
print("\nFile Parameters")
print("project ID:", proj_id)
print("structure:", structure)
print("topology:", topology)
print("mdp:", mdp)
print("index:", index_)
print("\nSimulation Parameters")
print("atoms:",atoms)
print("steps:",steps)
print("timestep:", timestep, "ns")
print("time:", simulation_time, "ns")
print("\nFAH Parameters")
print("points:", points)
print("deadline:", deadline)
print("timeout:", timeout)

# Prompt to generate project.xml based on project summary, and write if yes

fout = open('project.xml', 'w')
if project_type == 'gromacs': # re-add frame*.trr under the <return> flag
    fout.write("""<project type="GRO_A7" id="%s">
  <min-core-version v="0.0.17"/>

  <!-- project settings -->
  <runs v="10"/>
  <clones v="100"/>
  <gens v="100"/>
  <atoms v="%d"/>
  <delete-old-tpr/>
  <delete-old-trr/>

  <!-- simulation time == %f ns -->

  <!-- stats -->
  <stats_credit v="%d"/>
  <timeout v="%f"/>
  <deadline v="%f"/>

  <contact v="voelz@temple.edu"/>
  <accept-mode v="assign"/>

  <send>
    frame$gen.tpr
  </send>

  <return>
    md.log
    *.xtc
    *.xvg
  </return>

  <core-args>
    -s frame$gen.tpr
    -o frame$gen.trr
  </core-args>

  <create-command>
    /usr/local/bin/gromacs-5.0.4/bin/grompp -c $home/%s -f $home/%s -p $home/%s -n $home/%s -o $jobdir/frame0.tpr -po $jobdir/mdout.mdp -maxwarn 1
  </create-command>

  <next-gen-command>
    /usr/local/bin/gromacs-5.0.4/bin/convert-tpr -s $jobdir/frame0.tpr -f $results/frame$prev-gen.trr -o $jobdir/frame$gen.tpr -extend %d
  </next-gen-command>
</project>""" % (proj_id, proj_id, atoms, simulation_time, points, timeout, deadline, proj_id, structure, mdp, topology, index_, simulation_time*1000))

if project_type == 'openmm':
            fout.write("""<project type="OPENMM_21" id="1170x">
  <gpu/>
  <min-core-version v="0.0.18"/>

  <runs v="1"/>
  <clones v="500"/>
  <gens v="100"/>
  <atoms v="%d"/>

  <!-- simulation time == %f ns --> 

  <!-- stats -->
  <timeout v="%f"/>
  <deadline v="%f"/>
  <k-factor v="0.75"/>
  <stats-credit v="%d"/>

  <description v="%s %s"/>
  <contact v="voelz@temple.edu"/>
  <accept-mode v="assign"/>

  <send>
    $home/core.xml
    $home/SET8_Cat_Cycle/RUN$run/system.xml
    $home/SET8_Cat_Cycle/RUN$run/integrator.xml
    state.xml
  </send>

  <create-command>
    ln -s -f $home/SET8_Cat_Cycle/RUN$run/state.xml $jobdir/state.xml
  </create-command>

  <next-gen-command>
    mv -f $results/checkpointState.xml $jobdir/state.xml
  </next-gen-command>"""  % (proj_id, atoms, simulation_time, timeout, deadline, points, proj_id, description))
print("\nproject.xml generated. Please check that the number of runs, clones, gens is correct!\n")

# Prompt to test tpr by running a brief simulation on half the available cores
tpr_input = raw_input("\nWould you like to run a tpr test? (y/n): ")
if tpr_input == 'y':
    os.system("mkdir test_a_tpr; cp {" + structure + "," + mdp + "," + topology + "," + index_ + "} test_a_tpr; cp *itp test_a_tpr")
    for i in os.listdir(os.getcwd()):
        if "amber" in i:
            os.system("cp -r " + i + " test_a_tpr")
    os.chdir("test_a_tpr")
    os.system("/usr/local/bin/gromacs-5.0.4/bin/grompp -f " + mdp + " -c " + structure + " -p " + topology + " -n " + index_ + " -o " + " tpr_test.tpr -maxwarn 1")

    run_input = raw_input("\n Would you like to start a test run of your tpr? (y/n): ")
    if run_input == 'y':
        os.system("procs=$(nproc); gmx mdrun -v -nt $((procs/2)) -deffnm tpr_test")
        os.chdir("..")
        
