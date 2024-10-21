import os, sys, glob, subprocess

usage = """Usage: python length_dist_A8.py PROJDIR ns_per_frame nruns

    Calculates the distribution of trajectory lengths for a project.
    Example:
        python length_dist_A8.py ~/server2/data/PROJ12462 10.0 2
"""

if len(sys.argv) < 4:
    print(usage)
    sys.exit(1)

projdir = sys.argv[1]
ns_per_frame = float(sys.argv[2])
nruns = int(sys.argv[3])
Continue = True 
lengths = []

print('#frame\tlength(ns)\tN(frame)')
frame = 0
maxframe = 1000
total_ns = 0.
while Continue:

    n = 0
    # print('run'),
    for run in range(nruns): 
        try:
            xtcs_with_this_frame = glob.glob(os.path.join(projdir, f'RUN{run}/CLONE*/results{frame}/frame{frame}.xtc'))
            #print('xtcs_with_this_frame', xtcs_with_this_frame)
            n_thisrun = len(xtcs_with_this_frame)
        except:
            n_thisrun = 0 
        n += n_thisrun
        # print(run),  
    print('%d\t%d\t%d'%(frame, (frame+1)*ns_per_frame, n))
    total_ns += ns_per_frame*n
    if (n < 1) or frame >= maxframe:
        Continue = False
    else:
        frame +=1

print('Total ns in the dataset:', total_ns)

