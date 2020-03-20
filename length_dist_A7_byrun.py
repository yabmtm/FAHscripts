import os, sys, glob, commands

# Written by Hongbin Wan #

if len(sys.argv) < 4:
    print "Usage: python length_dist.py PROJDIR ns_per_frame nruns"
    sys.exit(1)

projdir = sys.argv[1]
ns_per_frame = float(sys.argv[2])
nruns = int(sys.argv[3])
lengths = []

for run in range(nruns):

    Continue = True
    outfile = 'run%d.txt'%run
    file = open(outfile, 'w')
    print '*** RUN %d ***'%run
    print '#frame\tlength(ns)\tN(frame)'
    frame = 0
    maxframe = 100
    total_ns = 0.
    while Continue:
      try:
        n = int(commands.getoutput('ls -1 %s | wc -l'%os.path.join(projdir, 'RUN%d/CLONE*/results%d/traj_comp.xtc'%(run,frame)) ))
      except Exception as e:
        n = 0 
      file.write('%d\t%d\t%d\n'%(frame, (frame+1)*ns_per_frame, n))
      total_ns += ns_per_frame*n
      if (n < 1) or frame >= maxframe:
        Continue = False
      else:
        frame +=1
    file.close()
    print('wrote run%d.txt'%run)
print 'Total ns in the dataset:', total_ns

