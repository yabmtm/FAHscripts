#!/home/server/anaconda2/bin/python

import glob, subprocess, re, shutil, os, sys

# get list of active projects
active_projects = sorted(glob.glob('/home/server/server2/projects/p*'))
active_project_numbers = [re.sub('.*projects/','',x) for x in active_projects]
gro_projects, omm_projects = [],[]

# separate into Gromacs and OpenMM
for i in active_projects:
    with open('%s/project.xml'%i) as f:
        lines = f.readlines()[0]
    if 'GRO' in lines:
        gro_projects.append(i)
    elif 'OPENMM' in lines:
        omm_projects.append(i)

# archive old project pages
web_projects = glob.glob('/home/server/server2/pub/data/Gromacs/p*') + glob.glob(
    'home/server/server2/pub/data/OpenMM/p*')
for i in web_projects:
    project_number = re.sub('.*data/.*/','',i)
    if project_number not in active_project_numbers:
        print('Archiving %s'%project_number)
        shutil.move(i,'/home/server/server2/pub/data/Archive/%s'%project_number)

# make sure ngl.html exists for all active projects
for i in web_projects:
    if not os.path.exists('%s/ngl.html'%i):
        shutil.copy2('/home/server/server2/pub/data/ngl.html',i)

# set number of projects per column in dropdown menu
entries_per_column = 5
gro_cols = len(gro_projects)/entries_per_column + 1
omm_cols = len(omm_projects)/entries_per_column + 1
gro_str, omm_str  = '',''

# write out html for gromacs dropdown
for i in range(gro_cols):
    gro_str += '\t\t\t\t\t\t<ul>\n'
    for j in range(entries_per_column):
        try:
            proj_id = re.sub('.*projects/','',gro_projects[(entries_per_column*i)+j])
            if not os.path.exists('/home/server/server2/pub/data/Gromacs/%s'%proj_id):
                os.makedirs('/home/server/server2/pub/data/Gromacs/%s'%proj_id)
            gro_str += '\t\t\t\t\t\t\t<li><a href="/pub/data/Gromacs/%s/index.html">%s</a>\n'%(proj_id,proj_id)
        except Exception as e:
            break
    gro_str += '\t\t\t\t\t\t</ul>\n'

# write out html for openmm dropdown
for i in range(omm_cols):
    omm_str += '\t\t\t\t\t\t<ul>\n'
    for j in range(entries_per_column):
        try:
            proj_id = re.sub('.*projects/','',omm_projects[(entries_per_column*i)+j])
            if not os.path.exists('/home/server/server2/pub/data/OpenMM/%s'%proj_id):
                os.makedirs('/home/server/server2/pub/data/OpenMM/%s'%proj_id)
            omm_str += '\t\t\t\t\t\t\t<li><a href="/pub/data/OpenMM/%s/index.html">%s</a>\n'%(proj_id,proj_id)
        except Exception as e:
            break
    omm_str += '\t\t\t\t\t\t</ul>\n'

active_projects = ''.join(['      <a href="#">%s</a>\n'%re.sub('.*projects/','',i) for i in active_projects])

# extract project statistics from project.xml
for project in sorted(glob.glob('/home/server/server2/pub/data/Gromacs/p*')):
    for page in ['index.html','plots.html','structure.html']:
        proj_id, footer = re.sub('.*p','',project), ''
        html_str = "<h1 style='padding-top:9em; margin-left:20px'>PROJ%s</h1>\n  <div style='display:table-cell; float:left; margin:50px'>\n  <ul>\n"%proj_id

        cmd = 'grep runs /home/server/server2/projects/Gromacs/p%s/project.xml | sed "s/.*=.//" | sed "s/\\".*//"'%proj_id
        html_str += '    <li>Runs: %s</li>\n'%(subprocess.check_output(cmd, shell=True).strip('\n'))
        cmd = 'grep clones /home/server/server2/projects/Gromacs/p%s/project.xml | sed "s/.*=.//" | sed "s/\\".*//"'%proj_id
        html_str += '    <li>Clones: %s</li>\n'%(subprocess.check_output(cmd, shell=True).strip('\n'))
        cmd = 'grep gens /home/server/server2/projects/Gromacs/p%s/project.xml | sed "s/.*=.//" | sed "s/\\".*//"'%proj_id
        html_str += '    <li>Gens: %s</li>\n'%(subprocess.check_output(cmd, shell=True).strip('\n'))
        cmd = 'grep atoms /home/server/server2/projects/Gromacs/p%s/project.xml | sed "s/.*=.//" | sed "s/\\".*//"'%proj_id
        html_str += '    <li>Atoms: %s</li>\n'%(subprocess.check_output(cmd, shell=True).strip('\n'))
        cmd = 'grep stats_credit /home/server/server2/projects/Gromacs/p%s/project.xml | sed "s/.*=.//" | sed "s/\\".*//"'%proj_id
        html_str += '    <li>Credit: %s</li>\n'%(subprocess.check_output(cmd, shell=True).strip('\n'))
        cmd = 'grep contact /home/server/server2/projects/Gromacs/p%s/project.xml | sed "s/.*=.//" | sed "s/\\".*//"'%proj_id
        html_str += '    <li>Contact: %s</li>\n'%(subprocess.check_output(cmd, shell=True).strip('\n'))

# find mdp file (this could be done better)
        try:
            for file in os.listdir("/home/server/server2/projects/p%s/RUN0"%proj_id):
                if file.endswith(".mdp") and not file.startswith('mdout'):
                    file = '/home/server/server2/projects/p%s/RUN0/%s'%(proj_id,file)
                    break
        except Exception as e:
            for file in os.listdir("/home/server/server2/projects/p%s"%proj_id):
                if file.endswith(".mdp") and not file.startswith('mdout'):
                    file = '/home/server/server2/projects/p%s/%s'%(proj_id,file)
                    break

# extract project statistics from mdp file
        cmd = 'grep dt %s | sed "s/.*= //" | sed "s/;.*//"'%file
        timestep = float(subprocess.check_output(cmd, shell=True).strip('\n'))*1000
        html_str += '    <li>Timestep: %s fs</li>\n'%timestep
        cmd = 'grep nsteps %s | sed "s/.*= //" | sed "s/;.*//"'%file
        html_str += '    <li>Time/WU: %s ps</li>\n'%(int(subprocess.check_output(cmd, shell=True).strip('\n'))*timestep/1000)
        cmd = 'grep ^nstxout %s | sed "s/.*= //" | sed "s/;.*//"'%file
        html_str += '    <li>TRR Frequency: %s ps</li>\n'%(int(subprocess.check_output(cmd, shell=True).strip('\n'))*timestep/1000)
        cmd = 'grep nstxtcout %s | sed "s/.*= //" | sed "s/;.*//"'%file
        html_str += '    <li>XTC Frequency: %s ps</li>\n'%(int(subprocess.check_output(cmd, shell=True).strip('\n'))*timestep/1000) 

# extract current distribution of project lengths (and add links)
        cmd = 'grep Total %s/lengths.html 2>/dev/null'%project
    
        try:
            html_str += '    <ul>\n        <li><a href="index.html" style="color:#474747";>Initial Structure</a></li>\n'
            html_str += '        <li><a href="plots.html" style="color:#474747";>Plots</a></li>\n'
            html_str += '        <li><a href="ngl.html" style="color:#474747";>NGLView</a></li>\n'
            html_str += re.sub('<li>','    <li><a href="lengths.html" style="color:#474747";>',re.sub('</li>','</a></li>',subprocess.check_output(cmd, shell=True)))
            big_images = '    <img id="default" src="lengths/lengths.png" />\n'

        except Exception as e:
            big_images = ''

        thumbnails = ''
        html_str += '  </ul>\n  </div>\n'

# write html for ngl viewer
        if page == 'index.html':
            try:
                cmd = 'find /home/server/server2/projects/Gromacs/p%s | grep xtc.gro | head -n 1'%proj_id
                gro_file = subprocess.check_output(cmd, shell=True).strip('\n')
                shutil.copy2(gro_file, '/home/server/server2/pub/data/Gromacs/p%s/p%s.gro'%(proj_id,proj_id))
                html_str += '\n  <script src="../../node_modules/ngl/dist/ngl.js"></script>\n  <script>\n    var stage;\n'
                html_str += '    document.addEventListener("DOMContentLoaded", function () {\n'
                html_str += '      stage = new NGL.Stage("viewport");\n'
                html_str += '      stage.loadFile("p%s.gro", {defaultRepresentation: true});\n    });\n  </script>\n'%proj_id
                html_str += '  <div id="viewport" style="position:relative; width:1000px; height:500px; border: 3px green; float:left; margin-left:64;"></div>\n'
        
            except Exception as e:
                print(e)
                footer += '<li><b>Unable to find xtc.gro. Please put this structure file in your project directory.</b></li>\n'
                pass

# process tica plots into gallery if they exist
        if page == 'plots.html':
            if os.path.exists('/home/server/server2/pub/data/Gromacs/p%s/tica'%proj_id):
                tica_images = sorted(glob.glob('/home/server/server2/pub/data/Gromacs/p%s/tica/plots/*png'%proj_id))
                tica_images = [re.sub('.*plots/','',x) for x in tica_images]

                for i in range(len(tica_images)):
                    big_images += '    <img id="image%d" src="tica/plots/%s" />\n'%(i+1,tica_images[i])
                    thumbnails += '    <li><a href="#image%d"><img src="tica/plots/%s" /></a></li>\n'%(i+1,tica_images[i])

            if thumbnails:
                thumbnails += '    <li><a href="#image%d"><img src="lengths/lengths.png" /></a></li>\n'%(i+2)

            html_str += '<link rel="stylesheet" href="../../../css/gallery.css">\n<div class="title"></div>\n<div class="image-gallery">\n'
            html_str += '  <div class="big-image">\n'
            html_str += '%s\n  </div>\n  '%big_images #+ '<a class="button prev"><</a>\n  <a class="button next">></a>\n'
            html_str += '  <div class="thumbs">\n  <ul>\n%s\n  </ul>\n  </div>\n</div>\n'%thumbnails

# write main body of each html for each project
        with open('/home/server/server2/pub/data/Gromacs/p%s/%s'%(proj_id,page), 'w') as f:
            f.write("""<!DOCTYPE html>
<html lang="en">
<head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, user-scalable=no, minimum-scale=1.0, maximum-scale=1.0">
        <title>p%s</title>
        <link rel="stylesheet" href="../../../css/style.css">
  <link href="https://fonts.googleapis.com/css?family=Montserrat" rel="stylesheet">
</head>
<body style="background-color:#DFDFDF;">
        <body link="#000000" vlink="#808080" alink="#FF0000">
        <header>
                <nav>
                        <ul id="main-menu"><!-- ul for the top menu items -->
                                <li><a href="../../../dashboard.html">Dashboard</a></li>
                                <li><a href="">Gromacs Projects</a>
                                        <div class="fw-dropdown"><!-- start of mega menu dropdown that appears on hover, fw = full width -->
%s
                                        </div><!-- end of the mega menu dropdown that appears on hover -->
                                </li>
                                <li><a href="">OpenMM Projects</a>
                                        <div class="fw-dropdown">
                                                <ul class="dropdown-column-list"><!-- start of mega menu dropdown that appears on hover -->
%s
                                        </div><!-- end of mega menu dropdown that appears on hover -->
                                </li>
                                <li><a href="index.html">Logout</a></li>
                        </ul>
                </nav>
        </header>
%s

<footer>
%s
</footer>
</body>
</html>"""%(proj_id,gro_str,omm_str,html_str,footer))


