import os, sys, glob, commands

# Written by Hongbin Wan #

if len(sys.argv) < 3:
    print "Usage: python length_dist.py PROJDIR ns_per_frame"
    sys.exit(1)

projdir = sys.argv[1]
ns_per_frame = float(sys.argv[2])
Continue = True 
lengths = []

print '#frame\tlength(ns)\tN(frame)'
frame = 0
maxframe = 100
total_ns = 0.
while Continue:
    try:
        n = int(commands.getoutput('ls -1 %s | wc -l'%os.path.join(projdir, 'RUN*/CLONE*/results%d/traj_comp.xtc'%frame) ))
    except:
        n = 0 
    print '%d\t%d\t%d'%(frame, (frame+1)*ns_per_frame, n)
    total_ns += ns_per_frame*n
    if (n < 1) or frame >= maxframe:
        Continue = False
    else:
        frame +=1

print 'Total ns in the dataset:', total_ns

