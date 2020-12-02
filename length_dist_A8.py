import os, sys, glob, commands

if len(sys.argv) < 4:
    print "Usage: python2 length_dist_A8.py PROJDIR ns_per_frame nruns"
    sys.exit(1)

projdir = sys.argv[1]
ns_per_frame = float(sys.argv[2])
nruns = int(sys.argv[3])
Continue = True 
lengths = []

print '#frame\tlength(ns)\tN(frame)'
frame = 0
maxframe = 1000
total_ns = 0.
while Continue:

    n = 0
    # print 'run',
    for run in range(nruns): 
        try:
            n_thisrun = int(commands.getoutput('ls -1 %s | wc -l'%os.path.join(projdir, 'RUN%d/CLONE*/results%d/frame%d.xtc'%(run,frame,frame)) ))
        except:
            n_thisrun = 0 
        n += n_thisrun
        # print run,  
    # print
    print '%d\t%d\t%d'%(frame, (frame+1)*ns_per_frame, n)
    total_ns += ns_per_frame*n
    if (n < 1) or frame >= maxframe:
        Continue = False
    else:
        frame +=1

print 'Total ns in the dataset:', total_ns

