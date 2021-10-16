import os, sys, glob
import numpy as np

# import xvg_tools
import argparse, textwrap

parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent('''\

    Shows progress bar of how much sampling the simulations have reached

    EXAMPLE
    $ python progress_bar.py ~/server2/data/SVR2166411725 18433 -r 0,1,3-4

      ''' ))

parser.add_argument('datadir', type=str, help='The direcory where projdata is saved')
parser.add_argument('projnum', type=int, help='The project number')
#parser.add_argument('--agg', dest='agg', action='store_true',
#                    help='Specify the structure comes from an AGG simulation.  The resnum in conf.gro is off by +1!')
parser.add_argument('-r', action='append',
                    dest='specific_runs',
                    default=[],
                    help='A string (with no spaces! e.g. 0,1,3-4) denoting a specific set of runs')

args = parser.parse_args()
print('args.datadir', args.datadir)
print('args.projnum', args.projnum)
print('args.specific_runs', args.specific_runs)

## If the user has chosen specific runs, parse this string
myruns = None
if len(args.specific_runs) > 0:
    use_specific_runs = True
    myruns = []
    fields = args.specific_runs[0].split(',')
    for field in fields:
        if field.count('-') == 0:
            myruns.append(int(field))
        else:
            ends = [int(s) for s in field.split('-')]
            myruns += list(range(ends[0],ends[1]+1)) 
print('myruns', myruns)



def project_length_in_ns(projnum):
    """Returns a float with the project length in nanoseconds (ns)"""

    w = {}

    # sampl9 challenges
    w[18433] = 10.0   # RL -  8500 atoms
    w[18434] = 10.0   # L  -  4091 atoms

    # sampl9 challenges
    w[16958] = 10.0   # chignolin in different FFs 

    return w[projnum]



# find the number of runs from the project.xml   file
project_xml_file = '/home/server/server2/projects/p%d/project.xml'%args.projnum
if not os.path.exists(project_xml_file):
    # ns338286 projects do not have a leading "p" .... try that
    project_xml_file = '/home/server/server2/projects/%d/project.xml'%args.projnum
fin = open(project_xml_file, 'r')
lines  = fin.readlines()
for line in lines:
    if line.count('runs'):
        fields = line.split('"')
        nruns = int(fields[1])


# Find the viable clones 
if myruns == None:
    myruns = range(0,nruns)
rundirs = [ os.path.join(args.datadir, 'PROJ%d/RUN%d'%(args.projnum,run)) for run in myruns ]

print()
print('PROJECT', args.projnum)

for i in range(len(myruns)):

    run = myruns[i]
    rundir = rundirs[i]
    clonedirs = glob.glob( os.path.join(rundir, 'CLONE*') )
    clones = [ int(os.path.basename(clonedir).replace('CLONE','')) for clonedir in clonedirs]

    # sort those clones!
    Isort = np.argsort(np.array(clones))
    clonedirs = [ clonedirs[j] for j in Isort]
    clones = [ clones[j] for j in Isort ]

    # Create a header ruler bar to fit the length of the WU
    run_label = 'RUN%d'%run
    run_header  = '%-8s | '%run_label
    nchar_limit = 120
    tick = 1
    while len(run_header) < nchar_limit:  
        # put ticks every multiple of 10
        tick_value = int(10*project_length_in_ns(args.projnum))
        if tick_value >= 100:
            if tick == 1:
                run_header += '......%3d|'%(tick*tick_value)
            else:
                run_header += 'ns....%3d|'%(tick*tick_value)
        else:
            if tick == 1:
                run_header += '.......%2d|'%(tick*tick_value)
            else:
                run_header += 'ns.....%2d|'%(tick*tick_value)
        tick += 1
    print(run_header)

    # for each clone, count the number of results dirs
    #     TODO grab the data in the dhdl.xvg and pullx.xvg in each gen

    for j in range(len(clones)):
        clone = clones[j]
        clonedir = clonedirs[j]

        ### Find all the result0, result1, etc. dirs !
        resultdirs = []
        resultdirs1 = glob.glob( os.path.join(clonedir, 'results?') )
        resultdirs1.sort()
        resultdirs += resultdirs1

        resultdirs2 = glob.glob( os.path.join(clonedir, 'results??') )  
        resultdirs2.sort()
        resultdirs += resultdirs2

        resultdirs3 = glob.glob( os.path.join(clonedir, 'results???') )  
        resultdirs3.sort()
        resultdirs += resultdirs3

        empty_last_results_dir = False
        if empty_last_results_dir:
            print('CLONE %04d'%clone, '*'*(int((len(resultdirs)-1))  ) )  # added -1 to account for empty results dir at end
        else:
            print('CLONE %04d'%clone, '*'*(int((len(resultdirs)))  ) ) 

    #    dhdl_xvgfiles = os.path.join(clonedir
    #try:
        
    





