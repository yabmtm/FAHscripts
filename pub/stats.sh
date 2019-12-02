#!/bin/bash
# Shell script to get disk usage, cpu usage, RAM usage,system load,etc.
# from multiple Linux servers and output the information on a single
# server in html format. Read below for usage/installation info
# *--------------------------------------------------------------------*
# * ORIGINAL WORK BY:
# * dig_remote_linux_server_information.bash,v0.1
# * Last updated on 25-Jul-2005*
# * Copyright (c) 2005 nixCraft project *
# * Comment/bugs: http://cyberciti.biz/fb/ *
# * Ref url: http://cyberciti.biz/nixcraft/forum/viewtopic.php?t=97 *
# * This script is licensed under GNU GPL version 2.0 or above *
# *--------------------------------------------------------------------*
# * Installation Info *
# ---------------------------------------------------------------------*
# You need to setup ssh-keys to avoid password prompt, see url how-to setup:
# ssh-keys:
# cyberciti.biz/nixcraft/vivek/blogger/2004/05/ssh-public-key-based-authentication.html
#
# [1] You need to setup correct VARIABLES script:
#
# (a) Change Q_HOST to query your host to get information
# Q_HOST="192.168.1.2 127.0.0.1 192.168.1.2"
#
# (b) Setup USR, who is used to connect via ssh and already setup to connect
# via ssh-keys
# USR="nixcraft"
#
# (c)Show warning if server load average is below the limit for last 5 minute.
# setup LOAD_WARN as per your need, default is 5.0
#
# LOAD_WARN=5.0
#
# (d) Setup your network title using MYNETINFO
# MYNETINFO="My Network Info"
#
# (e) Save the file
#
# Please refer to forum topic on this script:
#
# ----------------------------------------------------------------------
# Execute script as follows:
# this.script.name > /var/www/html/info.html
# If you generate a cron task (using "crontab -e") you can create
# a webpage every x minutes to show always updated information
# On the crontab file you should have something similar to:
# */5 * * * * /usr/local/bin/gen-statistics.sh /var/www/info.html
# ======================================================================
# This script is part of nixCraft shell script collection (NSSC)
# Visit http://bash.cyberciti.biz/ for more information.
# ----------------------------------------------------------------------
# Modified by Sergio Oller
#   - Use HTML5 and CSS instead of tables
#   - Use more colors for DISK, RAM and CPU thresholds.
#   - Get number of cores to know maximum CPU load value.
#   - Show only local mounted disks.
#   - Use LANG=C before "df" command to parse properly the results
# TODO:
#   - Use a single SSH connection for each HOST. This will make 
#     the script MUCH faster.
#   - Use cgi on a web server to regenerate the page on request, not
#     on crontab (idea: if no one is looking, do not update the info)
# SSH SERVER HOST IPS, setup me
# Change this to query your host
Q_HOST=$(hostname)

# SSH USER, change me
USR="server"

# Show warning if server load average is above the limit (5 min avg)
LOAD_GREEN=0.2
LOAD_RED=0.8

# Your network info
MYNETINFO=$(hostname)
#
PBY='This information is updated every minute'

# Local path to ssh and other bins
PING="/bin/ping"
NOW="$(date -R)"

## functions ##
writeHead(){
echo '
<!DOCTYPE html>
<html>
   <head><title>WS Status</title>
   <meta http-equiv="refresh" content="60">
   <style type="text/css">
table.gridtable {
	font-family: verdana,arial,sans-serif;
	font-size:11px;
	color:#333333;
	border-width: 1px;
	border-color: #666666;
	border-collapse: collapse;
}
table.gridtable th {
	border-width: 1px;
	padding: 8px;
	border-style: solid;
	border-color: #666666;
	background-color: #dedede;
}
table.gridtable td {
	border-width: 1px;
	padding: 8px;
	border-style: solid;
	border-color: #666666;
	background-color: #ffffff;
}
tr.redrow {
   color: red;
}
tr.orangerow {
   color: orange;
}
</style>
   </head>
   <body>'
echo "<h1 style='padding-top:8.5em'>Dashboard</h1>"
echo "<p>Generated on $NOW.</p>"
echo "<p>Hostname $MYNETINFO</p>"
}

writeFoot(){
echo "<p>This information is updated every minute</p>"
echo "</body></html>"
}

process_host(){
myhost="$1"
echo '<div style="display: table-cell; left-margin:24px;">'

######
_CMD=""
rhostname="$(hostname)"

ruptime="$($_CMD uptime)"
ruptime2="$ruptime"
if $(echo $ruptime | grep -E "min|days" >/dev/null); then
    x=$(echo $ruptime | awk '{ print $3 $4}')
else
    x=$(echo $ruptime | sed s/,//g| awk '{ print $3 " (hh:mm)"}')
fi
ruptime="$x"

rloadout="$(echo ${ruptime2} |awk -F'average:' '{ print $2}')"
rload=1
x="$(echo $rloadout | sed s/,//g | awk '{ print $2}')"
numcores=$($_CMD cat /proc/cpuinfo | grep 'model name' | wc -l)
ygreen="$(echo "$x <= ($LOAD_GREEN*$numcores)" | bc)"
yred="$(echo "$x >= ($LOAD_RED*$numcores)" | bc)"

[ "$ygreen" == "1" ] && rload="<span style=\"color:green\">$rloadout (max: $numcores)</span>"
[ "$yred" == "1" ] &&   rload="<span style=\"color:red\">$rloadout (max: $numcores)</span>"
[ "$rload" == "1" ] &&  rload="<span style=\"color:orange\">$rloadout (max: $numcores)</span>"

rclock="$(date +"%r")"
rtotalprocess="$($_CMD ps axue | grep -vE "^USER|grep|ps" | wc -l)"
rfs="$(df -hTl 2>/dev/null | grep -vE "^Filesystem|shm" |
 awk ' \
  {w=sprintf("%d",$6); \
   if (w > 90) {
      print "<tr class=\"redrow\"><td>" $7 "</td><td>" $6 "<td>" $2 "</td><td>" $4"/"$3 "</td></tr>" \
   } else if (w > 80) {
      print "<tr class=\"orangerow\"><td>" $7 "</td><td>" $6 "<td>" $2 "</td><td>" $4"/"$3 "</td></tr>" \
   } else {
      print "<tr><td>" $7 "</td><td>" $6 "<td>" $2 "</td><td>" $4"/"$3 "</td></tr>" \
   }
  }')"

rusedram="$(LANG=C free -mt | grep 'buffers/cache' | awk '{print $3}')"
rfreeramout="$(LANG=C free -mt | grep 'buffers/cache' | awk '{print $4}')"
rtotalram="$(LANG=C free -mto | grep Mem: | awk '{ print $2 }')"

rfreeram="1"
rfreeramred="$(echo "$rfreeramout < 0.01*${rtotalram}" | bc)"
rfreeramyellow="$(echo "$rfreeramout < 0.15*${rtotalram}" | bc)"
[ "$rfreeramyellow" == "1" ] &&   rfreeram="<span style=\"color:orange\">$rfreeramout MB</span>"
[ "$rfreeramred" == "1" ] && rfreeram="<span style=\"color:red\">$rfreeramout MB</span>"
[ "$rfreeram" == "1" ] &&  rfreeram="<span style=\"color:green\">$rfreeramout MB</span>"

rusedram="$rusedram MB"
rtotalram="$rtotalram MB"

rusedswap="$(LANG=C free -mto | grep Swap: | awk '{ print $3 }')"
rfreeswapout="$(LANG=C free -mto | grep Swap: | awk '{ print $4  }')"
rtotalswap="$(LANG=C free -mto | grep Swap: | awk '{ print $2 }')"

rfreeswap="1"
rfreeswapred="$(echo "$rfreeswapout < 0.01*${rtotalswap}" | bc)"
rfreeswapyellow="$(echo "$rfreeswapout < 0.15*${rtotalswap}" | bc)"
[ "$rfreeswapyellow" == "1" ] &&   rfreeswap="<span style=\"color:orange\">$rfreeswapout MB</span>"
[ "$rfreeswapred" == "1" ] && rfreeswap="<span style=\"color:red\">$rfreeswapout MB</span>"
[ "$rfreeswap" == "1" ] &&  rfreeswap="<span style=\"color:green\">$rfreeswapout MB</span>"


rusedswap="$rusedswap MB"
rtotalswap="$rtotalswap MB"

user_processes_body=$(for cc in $(LANG=C ps hax -o uid,user | sort | uniq | awk -F' ' '{if ($1 > 1000 && $1 < 10000) print $2}'); do percent=$(top -b -n 1 -u "$cc" | awk 'NR>7 { sum += $9; } END { print sum/'"${numcores}"'; }') ; percentram=$(top -b -n 1 -u "$cc"  | awk 'NR>7 {sum += $10;} END {print sum;}') ; echo "<tr><td>$cc</td><td>$percent %</td><td>${percentram} %</td></tr>"; done)

$PING -c1 $myhost>/dev/null
if [ "$?" != "0" ] ; then
rping="<span style=\"color:red\"> Failed </span>"
else
rping="<span style=\"color:green\"> Ok </span>"

echo -e "<ul>\n<li>Ping status: $rping</li>"
echo "<li>Time: $rclock</li>"
echo "<li>Uptime: $ruptime </li>"
echo "<li>CPU Load average (1m, 5m, 15m): $rload </li>"
echo "<li>Number of cores: $numcores</li>"
echo "<li>Total running processes: $rtotalprocess </li>"
echo "<li style='padding-top:1em'>Memory status:"
echo "<table class=\"gridtable\"><thead><tr><th>Type</th><th>Used</th><th>Free</th><th>Total</th></tr></thead>"
echo "<tbody><tr><td>RAM</td><td>$rusedram</td><td>$rfreeram</td><td>$rtotalram</td></tr>"
echo "       <tr><td>SWAP</td><td>$rusedswap</td><td>$rfreeswap</td><td>$rtotalswap</td></tr></tbody></table>"
echo "</li>"
echo "<li style='padding-top:1em'>Disk status:"
echo "<table class=\"gridtable\"><thead><tr><th>Mount point</th><th>Percent used</th><th>File System</th><th>Used/Total</th></tr></thead>"
echo "<tbody>"
echo "$rfs"
echo "</tbody></table>"
echo "</li>"
echo "</ul>"
fi
######


echo '</div>'
}
## main ##

main() {
writeHead
echo '<div style="width: 100%; display: table;">'
echo '    <div style="display: table-row">'
for host in $Q_HOST; do
   process_host $host
done
echo '    </div>'
echo '</div>'
writeFoot
}

OUTFILE="$1"
mytext=$(main)
echo "$mytext" #> "$OUTFILE"

exit 0
