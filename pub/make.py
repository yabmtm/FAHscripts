#!/home/server/anaconda2/bin/python

import glob, re, subprocess, os

active_projects = sorted(glob.glob('/home/server/server2/projects/p*'))
gro_projects, omm_projects = [],[]

for i in active_projects:
    with open('%s/project.xml'%i) as f:
        lines = f.readlines()[0]
    if 'GRO' in lines:
        gro_projects.append(i)
    elif 'OPENMM' in lines:
        omm_projects.append(i)
        
gro_cols = len(gro_projects)/5 + 1
omm_cols = len(omm_projects)/5 + 1
gro_str, omm_str  = '',''

for i in range(gro_cols):
    gro_str += '\t\t\t\t\t\t<ul>\n'
    for j in range(5):
        try:
            proj_id = re.sub('.*projects/','',gro_projects[(5*i)+j])
            if not os.path.exists('/home/server/server2/pub/data/Gromacs/%s'%proj_id):
                os.makedirs('/home/server/server2/pub/data/Gromacs/%s'%proj_id)
            gro_str += '\t\t\t\t\t\t\t<li><a href="data/Gromacs/%s/index.html">%s</a>\n'%(proj_id,proj_id)
        except Exception as e:
            break
    gro_str += '\t\t\t\t\t\t</ul>\n'

for i in range(omm_cols):
    omm_str += '\t\t\t\t\t\t<ul>\n'
    for j in range(5):
        try:
            proj_id = re.sub('.*projects/','',omm_projects[(5*i)+j])
            if not os.path.exists('/home/server/server2/pub/data/OpenMM/%s'%proj_id):
                os.makedirs('/home/server/server2/pub/data/OpenMM/%s'%proj_id)
            omm_str += '\t\t\t\t\t\t\t<li><a href="data/OpenMM/%s/%s.html">%s</a>\n'%(proj_id,proj_id,proj_id)
        except Exception as e:
            break
    omm_str += '\t\t\t\t\t\t</ul>\n'

active_projects = ''.join(['      <a href="#">%s</a>\n'%re.sub('.*projects/','',i) for i in active_projects])
dashboard = subprocess.check_output('/home/server/server2/pub/stats.sh', shell=True)
hostname = subprocess.check_output('hostname', shell=True).strip('\n')
log = subprocess.check_output('tail -n 500 /home/server/server2/fah-work.log| sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[mGK]//g"', shell=True)

with open('/home/server/server2/pub/dashboard.html', 'w') as f:
    f.write("""<!DOCTYPE html>
<html lang="en">
<head>
	<meta charset="UTF-8">
	<title>%s</title>
	<link rel="stylesheet" href="css/style.css">
  <link href="https://fonts.googleapis.com/css?family=Montserrat" rel="stylesheet">
  <meta name="viewport" content="width=device-width">
</head>
<body style="background-color:#DFDFDF;">
	<header>
		<nav>
			<ul id="main-menu"><!-- ul for the top menu items -->
                                <li><a href="dashboard.html">Dashboard</a></li>
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
<div class="container" style="position:fixed; top:34%%; left:27%%;">

  <section class="container__snippet snippet--regex">
    <p style="color:#FFFFFF">LOG</p>
    <pre>
      <code style="color:#FFFFFF">
        <!-- snippet  -->
      </code>
    </pre>
  </section>

</div>

<script src="js/log.js"></script>

</body>
</html>"""%(hostname,gro_str,omm_str,dashboard))








with open('/home/server/server2/pub/js/log.js', 'w') as f:
    f.write("""// target all code elements and add the example snippet to both
const snippets = document.querySelectorAll("pre code");

const code = `
%s
`;

snippets.forEach(snippet => snippet.textContent = code);


// target the snippet in the regex section
const snippetRegex = document.querySelector(".snippet--regex pre code");

// include the regular expressions to find specific sections of text
// include them in an array of object, each detailing the expression and the class which needs to be applied upon the strings matching the connecting expression
const regularExpressions = [
        { 
          expression: /([0-9]+:[0-9]+:[0-9]+:..:|(REQ|OUT)[0-9]+:|^............#[0-9]+:)/gi,
          class: "timestamps"
        },
        { 
          expression: /([0-9]+\.[0-9]+\.[0-9]+\.[0-9]+:[0-9]+.+|Connecting to|Client: '.+[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+'|(to |from ).+@.+|assign.+foldingathome.+|POST.+)/g,
          class: "ip_addrs"
        },
	{
          expression: /Executing.+/gi,
          class: "command"
	},
        { 
          expression: /(PRCG.[0-9]+.[0-9]+,[0-9]+,[0-9]+.|P[0-9]+.R[0-9]+.C[0-9]+.G[0-9]+|Job)/g,
          class: "projects"
        },
	{
	  expression: /Request.WORK_REQUEST|Request.WORK_RESULTS|WORK_REQUEST|WORK_RESULTS|WORK_ASSIGNMENT|WORK_ACK|Response|accepted|assigned|Updating.AS|KEY_ACCEPTED|Retrying f..... job/gi,
	  class: "good_response"
	},
	{
	  expression: /Request WORK_FAILED|WORK_FAILED|WORK_FAULTY|Client reported Failed Assignment|#[0-9]+|Request WORK_FAULTY|Core.+Assignment|gcq#.+/gi,
		class: "bad_response"
	},
	{
	  expression: /(Empty WU.+|Retrying failed job|.[0-9]+ previous retries\)|Secure.+|Registering.+CS.at|PLEASE_WAIT)/gi,
		class: "neutral_response"
	},
        {
          expression: /.+requesting shutdown|.+Clean exit/gi,
                class: "server_down"
        }
];

// loop through the array of regex and update the HTML structure of the snippet wrapping the strings matching the expressions in span elements, bearing a class paired to the expression itself
// class defined in the stylesheet to alter the appearance of the matching strings
regularExpressions.forEach((regularExpression) => snippetRegex.innerHTML = snippetRegex.innerHTML.replace(regularExpression.expression, `<span class=${regularExpression.class}>$&</span>`));
"""%log)
