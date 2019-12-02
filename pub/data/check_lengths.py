#!/home/server/anaconda2/bin/python

import glob, re, subprocess, tqdm, os, datetime
import matplotlib
matplotlib.use('agg')
import numpy as np
from matplotlib import pyplot as plt

for project in sorted(glob.glob('/home/server/server2/pub/data/Gromacs/p*')):
    if not os.path.exists('%s/lengths'%project):
        os.makedirs('%s/lengths'%project)

    proj_id = re.sub('.*p','',project)
    cmd = 'grep Runs %s/index.html | sed "s/.*Runs: //" | sed "s/<.*//"'%(project)
    runs = int(subprocess.check_output(cmd, shell=True))
    cmd = 'grep Clones %s/index.html | sed "s/.*Clones: //" | sed "s/<.*//"'%(project)
    clones = int(subprocess.check_output(cmd, shell=True))
    cmd = 'grep Gens %s/index.html | sed "s/.*Gens: //" | sed "s/<.*//"'%(project)
    gens = int(subprocess.check_output(cmd, shell=True))
    cmd = 'grep "Time/WU" %s/index.html | sed "s/.*WU: //" | sed "s/ ps.*//"'%(project)
    ns_per_frame = float(subprocess.check_output(cmd, shell=True))/1000.
    output, total_ns = '',0
    X,Y = [],[]

    for gen in tqdm.tqdm(range(gens)):
        cmd = 'ls -1 /home/server/server2/data/SVR166219/PROJ%s/RUN*/CLONE*/results%d/traj_comp.xtc 2>/dev/null | wc -l'%(proj_id,gen)
        count = int(subprocess.check_output(cmd, shell=True))
        if count:
            output += '    <li>%s ns: %s</li>\n'%(((gen+1)*ns_per_frame), count)
            X.append((gen+1)*ns_per_frame)
            Y.append(count)
            total_ns += ns_per_frame * count

    plt.bar(X,Y)
#    plt.xticks([str(x) for x in X])
    plt.xlabel('Trajectory length (ns)')
    plt.ylabel('Counts')
    plt.title('%s length distribution: %s%% complete'%(proj_id,float((100*np.sum(Y)/(runs*clones*gens)))))
    plt.savefig('%s/lengths/%s.png'%(project,'lengths')) #datetime.date.today().strftime("%m_%d_%Y")))
    plt.close()

    with open('%s/lengths.html'%project,'w') as f:
        f.write('''<!DOCTYPE html>
<html lang="en">
<head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, user-scalable=no, minimum-scale=1.0, maximum-scale=1.0">
        <title>p14086</title>
        <link rel="stylesheet" href="../../../css/style.css">
  <link href="https://fonts.googleapis.com/css?family=Montserrat" rel="stylesheet">
</head>
<body>''')
        f.write('<li><b><a href="index.html" style="color:#474747";>Back</a></b></li>\n'%proj_id)
        f.write(output)
        f.write('    <li>Total ns in dataset: %s</li>\n</body>\n</html>'%total_ns)

    print('Total ns in %s: %s\n'%(proj_id,total_ns))
