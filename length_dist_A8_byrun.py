import os, sys, glob, commands

# Written by Hongbin Wan #
# Updated by Steven Goold and Vincent Voelz 8/20/2021

if len(sys.argv) < 5:
    print """Usage: python2 length_dist_A8_byrun.py PROJDIR ns_per_frame nruns outprefix

    INPUTS

        PROJDIR       - the name of the data directory for the project (e.g. PROJ16969)
        ns_per_frame  - the number of nanoseconds for each WU (consult the *.mdp file in ~/server2/projects/pXXXXX)
        nruns         - the number of runs, i.e.  RUN0 ... RUN(nruns-1), to analyze
        outprefix     - Output files will be written  <outprefix>.run*.txt

    NOTE

        This script should only be run python v2

    EXAMPLE

        Try this: $ python2 length_dist_A8_byrun.py PROJ16969 5.0 5 p16969

    """

    sys.exit(1)

projdir = sys.argv[1]
ns_per_frame = float(sys.argv[2])
nruns = int(sys.argv[3])
outprefix = sys.argv[4]
lengths = []

DEBUG = False

for run in range(nruns):

    Continue = True

    # Create and open the output file for this RUN
    outfile = '%s.run%d.txt'%(outprefix,run)
    file = open(outfile, 'w')
    print '# Writing to', outfile, '...' 
    print '*** RUN %d ***'%run

    header = '#frame\tlength(ns)\tN(frame)'
    print header
    file.write(header+'\n') 

    # For each gen (i.e. frame0.xtc, frame1.xtc, etc), count how many WU have been returned  
    frame = 0
    maxframe = 1000
    total_ns = 0.
    while Continue:
        try:
            cmd = 'ls -1 %s | wc -l'%os.path.join(projdir, 'RUN%d/CLONE*/results%d/frame%d.xtc'%(run,frame,frame))  
            if DEBUG:
		print '>>', cmd
            n = int(commands.getoutput(cmd))
            if DEBUG:
		print '\t%d clones have reached gen %d'%(n,frame)
        except Exception as e:
            n = 0 
        outstring = '%d\t%4.2f\t%d'%(frame, (frame+1)*ns_per_frame, n) 
        print outstring
        file.write(outstring+'\n')

        # Keep track of the total about of simulation time
        total_ns += ns_per_frame*n

        # if n==0, then there are no more gens to find
        if (n < 1) or frame >= maxframe:
            Continue = False
        else:
            frame +=1
    file.close()
     
    print '...Wrote file %s.'%outfile

print 'Total ns in the dataset:', total_ns

