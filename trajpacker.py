import os, sys, glob, subprocess, shutil
import numpy


class TrajectoryPacker(object):
    """A class to make easier the packaging of FAH trajectory data."""

    def __init__(self, projnum, run, outdir, xtc_grofile, xtc_topfile,
                       gmx_path='/usr/local/gromacs-2020/bin/gmx',
                       data_dir='/home/server/server2/data/SVR2166411725',
                       setup_dir=None,
                       verbose=False, debug=False):
        """Initialize the TrajectoryPacker class."""


        self.projnum = projnum
        self.outdir  = outdir
        self.run     = run
        self.xtc_grofile = xtc_grofile
        self.xtc_topfile = xtc_topfile
        self.gmx_path = gmx_path
        self.data_dir = data_dir
        if setup_dir == None:
            self.setup_dir = f'/home/server/server2/projects/p{projnum}'
        else:
            self.setup_dir = setup_dir

        self.writedir = os.path.join(self.outdir, f'P{projnum}_R{run}')
        self.verbose = verbose
        self.debug = debug

    def run_cmd(self, cmd):
        """Run a shell command."""
        if self.verbose:
            print('>>', cmd )
        subprocess.check_output(cmd, stderr=subprocess.STDOUT,shell=True).decode().split('\n')


    def pack(self):
        """Pack up a concatenated trajectory."""

        # create output directories
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)   
        if not os.path.exists(self.writedir):
            os.mkdir(self.writedir)

        # copy the xtc_grofile to this directory so we can use it visualize the trajectories later
        copied_xtc_grofile  = os.path.join(self.writedir, os.path.basename(self.xtc_grofile) )
        shutil.copyfile(self.xtc_grofile, copied_xtc_grofile)

        ### Step 1: Build a custom *.tpr for the subset of atoms  in the *.xtc trajectories
            
        ## write a dummy *.mdp for minimization (we will need a tpr for trjconv)
        write_mdp_cmd = f'echo "integrator          = steep" > {self.writedir}/xtc.mdp'
        self.run_cmd(write_mdp_cmd)
            
        # gmx grompp to make a fake *.tpr
        tpr_file = f'{self.writedir}/xtc.tpr'
        make_xtctpr_cmd = f'{self.gmx_path} grompp -f {self.writedir}/xtc.mdp -c {self.xtc_grofile} -p {self.xtc_topfile} -o {tpr_file} -maxwarn 2'
        if self.verbose:
            print(f'### Building a fake *.tpr {tpr_file}...')
        self.run_cmd(make_xtctpr_cmd)

        ### Step 2:  Compile all the CLONE indices with results in them
        data_rundir = os.path.join(self.data_dir, f'PROJ{self.projnum}/RUN{self.run}')
        clonedirs = glob.glob( os.path.join(data_rundir, 'CLONE*') )
        if self.debug:
            print('clonedirs', clonedirs)
        clones_with_results = []
        for clonedir in clonedirs:
            if os.path.exists( os.path.join(clonedir, 'results0')):
                clone = int(os.path.basename(clonedir).replace('CLONE',''))
                clones_with_results.append(clone) 
        if self.verbose:
            print('clones_with_results', clones_with_results)

        # Concatenate all the *.xtc into a single *.xtc
        for clone in clones_with_results:
            clonedir = os.path.join(data_rundir, f'CLONE{clone}')
            concat_xtcfile   = os.path.join(self.writedir, f'CLONE{clone}.xtc')
            xtc_filelist = []
            xtc_filelist += glob.glob(os.path.join(clonedir, 'results?/*.xtc'))
            xtc_filelist += glob.glob(os.path.join(clonedir, 'results??/*.xtc'))
            xtc_filelist += glob.glob(os.path.join(clonedir, 'results???/*.xtc'))
            xtc_filelist += glob.glob(os.path.join(clonedir, 'results????/*.xtc'))
            ngens = len(xtc_filelist)
            print(f'Concatening {ngens} gens in {clonedir}...', end='')
            # print('xtc_filelist', xtc_filelist)
            xtc_string = ' '.join(xtc_filelist) 
            trjcat_cmd = f'{self.gmx_path} trjcat -o {concat_xtcfile} -f {xtc_string}'
            self.run_cmd(trjcat_cmd)

            # gmx trjconv for PBC correction
            pbc_concat_xtcfile = concat_xtcfile.replace('.xtc', '_pbc.xtc')
            pbc_concat_cmd = f'echo "0\n0\n" | {self.gmx_path} trjconv -f {concat_xtcfile} -s {tpr_file} -pbc mol -center -o {pbc_concat_xtcfile}'
            self.run_cmd(pbc_concat_cmd)
            print(f'...wrote {pbc_concat_xtcfile}')
        

def parse_projectxml(xmlfile):
    """Parses key/value information from the project.xml tags.

    RETURNS
    result      - a dictionary of keys and values."""

    # Example text from the project.xml
    """<project type="GRO_A8" id="18433">

  <!-- project settings -->
  <runs v="68"/>
  <clones v="50"/>
  <gens v="1000"/>
  <atoms v="8500"/>
  <job-priority v="-gen clone run"/>
  <stats_credit v="5400"/>
  <timeout v="2.0"/>
  <deadline v="3.5"/>
  <description v="SAMPL9 host-guest"/>
  <contact v="jason.pattis@temple.edu"/>
    """

    result = {}

    fin = open(xmlfile, 'r')
    lines = fin.readlines()
    fin.close()    

    # NOTE: we *only* want the tags with 'v="' in them ! **Not** the multi-line tags

    for line in lines:
        if line.count('v=') > 0:
            key, value = xmltag_to_keyval(line)
            result[key] = value

    return result



def xmltag_to_keyval(tag):
    stripped_tag = tag.strip().replace('<', '').strip()   # should give just 'runs v="68">', e.g.
    key = stripped_tag.split()[0]
    value  = (tag.split('"'))[1]
    return key, value

    

if __name__ == '__main__':

    usage = """ Usage: $  python trajpacker.py [projnum] [run] [OUTDIR] [xtc-grofile] [xtc-topfile]

    INPUTS
        projnum              the project number
        run	             the run number
        OUTDIR               the directory output path
        xtc-grofile          a *.gro file containing ONLY the atoms in the *.xtc
        xtc-topfile          a *.top file that works for this grofile

    OPTIONS
        -g GMXPATH           e.g. /usr/local/gromacs-2020/bin/gmx 
        -d DATADIR           e.g. ~/server2/data/SVR2166411725 
        -s SETUPDIR          e.g. ~/server2/projects/p18433
        -v                   verbose flag -- extra print statements for each step

    EXAMPLE
        $ python trajpacker.py 18433 1 /data/packed_trajectories ~/server2/projects/p18433/RUN1/xtc.gro ~/server2/projects/p18433/RUN1/xtc.top 

    DESCRIPTION

        This script will find all CLONES in the  specified runs and create a series of concatenated trajectories
        [OUTDIR]/P[projnum]_RUN[run]/CLONE*.xtc....

    """

    if len(sys.argv) < 5:
        print(usage)
        sys.exit(1)

    import argparse, textwrap
    
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=usage)
    parser.add_argument('projnum', type=int, help='The project number')
    parser.add_argument('run', type=int, help='The RUN number')
    parser.add_argument('outdir', type=str, help='The directory to save the trajecory data')
    parser.add_argument('xtc_grofile', type=str, help='a *.gro file containing ONLY the atoms in the *.xtc')
    parser.add_argument('xtc_topfile', type=str, help='a *.top file corresponding ONLY the atoms in the *.xtc')
    parser.add_argument('-g', action='append',
                    dest='gmx_path',
                    default=None,
                    help='The path of your gmx installation, e.g. /usr/local/gromacs-2020/bin/gmx')
    parser.add_argument('-d', action='append',
                    dest='data_dir',
                    default=None,
                    help='The path of your data directory')
    parser.add_argument('-s', action='append',
                    dest='setup_dir',
                    default=None, 
                    help='The path of your setup directory')
    parser.add_argument('-v', action='store_true',
                    dest='verbose',
                    help='verbose flag' )


    args = parser.parse_args()
    print('args.projnum', args.projnum)
    print('args.run', args.run)
    print('args.outdir', args.outdir)
    print('args.xtc_grofile', args.xtc_grofile)
    print('args.xtc_topfile', args.xtc_topfile)
    print('args.gmx_path', args.gmx_path)
    print('args.data_dir', args.data_dir)
    print('args.setup_dir', args.setup_dir)


    ### parse and/or look up all pathnames based on hostname 

    projnum = args.projnum
    run     = args.run
    outdir  = args.outdir
    gmx_path = args.gmx_path
    data_dir = args.data_dir
    setup_dir = args.setup_dir

    # Hostname definitions
    hostname = (subprocess.check_output('hostname', shell=True).decode()).split('.')[0]   # need a decode here
    hostname = hostname.strip()
    print('hostname', hostname)

    # Path definitions
    if hostname == 'vav22':
        if gmx_path == None:
            gmx_path = '/usr/local/gromacs-2020/bin/gmx'
        if setup_dir == None:
            setup_dir = f'/home/server/server2/projects/p{projnum}'
        if data_dir == None:
            data_dir = '/home/server/server2/data/SVR2166411725'

    elif hostname == 'vav21':
        if gmx_path == None:
            gmx_path = '/usr/local/gromacs-2020/bin/gmx'
        if setup_dir == None:
            setup_dir = f'/home/server/server2/projects/p{projnum}'
        if data_dir == None:
            data_dir = '/home/server/server2/data/SVR2166411724'

    else:
        print('This script does not have default paths for hostname {hostname}.  Try specifying paths with flags.  Exiting.')
        sys.exit(1)

    print(f'Using default paths for hostname = {hostname}:')
    print('gmx_path', gmx_path)
    print('data_dir', data_dir)
    print('setup_dir', setup_dir)


    # Instantiate a TrajectoryPacker() object
    t = TrajectoryPacker(projnum, run, outdir, args.xtc_grofile, args.xtc_topfile, verbose=args.verbose)
    t.pack()
