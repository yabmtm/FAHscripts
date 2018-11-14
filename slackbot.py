import glob, os, re, socket, time, io, subprocess, sys
from slackclient import SlackClient
import numpy as np

# Check if Slackbot is already running
cmd = 'ps aux | grep slack | wc -l'
output = subprocess.check_output(cmd, shell=True)
if int(output) > 3: # Seems silly, but it's true
    print("Slackbot is already running! Kill process and try again.")
    sys.exit(1)

# instantiate Slack client
slack_client = SlackClient(os.environ.get('SLACK_BOT_TOKEN'))

# Slackbot's user ID in Slack: value is assigned after the bot starts up
starterbot_id = None

# Constants
RTM_READ_DELAY = 1 # 1 second delay between reading from RTM
MENTION_REGEX = "^<@(|[WU].+?)>(.*)"
VERBOSE = False

# Used for debugging on matt's laptop
if socket.gethostname() == 'syzygy':
    debug_prefix = '/media/matt/ext'
else:
    debug_prefix = ''

# User and Hostname definitions
hostname = 'vav3'
user_list = ['tuf74040', 'rraddi', 'voelz', 'matthew.hurley']

# Path definitions
array0 = debug_prefix + '/array0/data/'
array1 = debug_prefix + '/array1/server2/data/SVR166219/'
project_config_prefix = debug_prefix + '/array1/server2/projects/Gromacs/'
results_prefix = debug_prefix + '/array1/server2/.results/'
timestr = time.strftime("%Y%m%d-%H%M%S")
logfile=open(results_prefix + "slack_log.txt", "a+")
project_numbers = []

# Populate list of project numbers for all projects parsed by projanalysis.py
for i in sorted(glob.glob(results_prefix + "*/*_log.txt")):
    project_file = re.sub(results_prefix + ".*/", '', i)
    project_number = re.sub("_log.txt", '', project_file)
    project_numbers.append(project_number)
    project_numbers = sorted(project_numbers)    

def parse_bot_commands(slack_events):
    """
        Parses a list of events coming from the Slack RTM API to find bot commands.
        If a bot command is found, this function returns a tuple of command and channel.
        If its not found, then this function returns None, None.
    """
    for event in slack_events:
        if VERBOSE:
            print(event)
        if event["type"] == "message" and not "subtype" in event:
            user_id, message = parse_direct_mention(event["text"])
            if user_id == starterbot_id:
                return message, event["channel"]
    return None, None


def parse_direct_mention(message_text):
    """
        Finds a direct mention (a mention that is at the beginning) in message text
        and returns the user ID which was mentioned. If there is no direct mention, returns None
    """
    matches = re.search(MENTION_REGEX, message_text)
    # the first group contains the username, the second group contains the remaining message
    return (matches.group(1), matches.group(2).strip()) if matches else (None, None)


def handle_command(command, channel):
    """
        Executes bot command if the command is known
    """

    # Default response is help text for the user
    default_response = "Not sure what you mean. Try *{}*.".format("help")

    # Defaults to above response if invalid command is entered
    response = None

    # FAH-specific commands

    if command.startswith('help'):
        response = "\n*FAHbot is running and able to give information about %s.*\nTry one of the following:" %hostname
        response += "\n>`@" + hostname + " project 12345` : Prints project information on project 12345"
        response += "\n>`@" + hostname + " user voelz` : Prints user information on user voelz"
        response += "\n>`@" + hostname + " usage` : Prints data usage for " + hostname
        response += "\n>`@" + hostname + " restart` : Restarts " + hostname +" server daemon"
        response += "\n\n To see all available users use: `@" + hostname + " user list`"
        response += "\n To see all available projects use: `@" + hostname + " project list`"
        response += "\n\n To change points for a project to a value or by percent, use one of:"
        response += "\n>`@" + hostname + " change 12345 +15%`"
        response += "\n>`@" + hostname + " change 12345 1337`"
        response += "\nNote that the server daemon must be restarted for stats changes to take effect."
        response += "\n\n *Future additions:*"
        response += "\n Daily server stat plots for performance and space"
        response += "\n Other commands if you have ideas."
        response += "\n\n Report bugs to Matt."

    if command.startswith('user'):
        if command.split()[1] in user_list:
            response = get_user_info(user_list[user_list.index(command.split()[1])])
        else:
            response = "Error: No such user, " + command.split()[1] 
        if command.split()[1] == 'list':
            response = "\n*" + hostname + " User List:*"
            for i in user_list:
                response += '\n' + i

    if command.startswith('project'):
        if command.split()[1] in project_numbers:
            response = get_project_info(command.split()[1])
        if command.split()[1] == 'list':
            response = "\n*" + hostname + " Project List:*"
            for project_number in project_numbers:
                project_xml = project_config_prefix + 'p' + project_number + "/project.xml"

                try: # parse xml for stuff
                    with open(project_xml, 'r') as file:
                        lines = [line.strip().split() for line in file.readlines()]
                    file.close()

                    for j in range(len(lines)):
                        if '<description' in lines[j]:
                            description = " ".join(lines[j][1:])[3:-3]

                    response += '\n' + project_number + ": " + description
                except:
                    response += '\n' + project_number + ": No project.xml for this project."
    if command.startswith('usage'):
        response = get_data_usage()

    if command.startswith('change'):
        response = change_stats(command.split()[1], command.split()[2])    

    if command.startswith('say'):
        cmd = 'cowsay "' + ' '.join(command.split()[1:]) + '"'
        response = subprocess.check_output(cmd, shell=True)

    if command == 'restart':
        response = daemon_restart()        

    if command == 'lunch':
        response = pick_my_lunch()
    # Sends the response back to the channel
    slack_client.api_call(
        "chat.postMessage",
        channel=channel,
        text=response or default_response
    )

#    else:
#        with open('/home/matt/Pictures/wallhaven-673983.png', 'rb') as f:
#            slack_client.api_call(
#                "files.upload",
#                channels=channel,
#                filename='wallpaper.png',
#                title='sampletitle',
#                initial_comment='sampletext',
#                file=f)
            
#### fah functions

def get_project_info(project_number):
    if project_number in project_numbers:
        project_xml = project_config_prefix + 'p' + project_number + "/project.xml"
        
        try: # parse xml for stuff 
            with open(project_xml, 'r') as file:
                lines = [line.strip().split() for line in file.readlines()]
            file.close()

            for j in range(len(lines)):
                if '<contact' in lines[j]:
                    user = re.sub(r'@.*', "", lines[j][1][3:])
                if '<runs' in lines[j]:
                    runs = lines[j][1][3:-3]
                if '<clones' in lines[j]:
                    clones = lines[j][1][3:-3]
                if '<gens' in lines[j]:
                    gens = lines[j][1][3:-3]
                if '<stats-credit' in lines[j] or '<stats_credit' in lines[j]:
                    credit = lines[j][1][3:-3]
                if '<timeout' in lines[j]:
                    timeout = lines[j][1][3:-3]
                if '<deadline' in lines[j]:
                    deadline = lines[j][1][3:-3]
                if '<description' in lines[j]:
                    description = " ".join(lines[j][1:])[3:-3]
                    
        except:
            response = "Error: No project.xml file found for project: " + project_number
            pass
            
        try: # Parse results dir for most recent statistics
            results = results_prefix + user + '/' + project_number + '_log.txt'

            with open(results, 'r') as file:
                lines = [line.strip().split() for line in file.readlines()]
            file.close()
            
            size = lines[-1][2]
            length = lines[-1][1]

            try:
                growth = float(lines[-1][1]) - float(lines[-2][1])

            except:
                growth = 0            
        
            response = "*Project: " + project_number + "*\n\tUser: " + user + "\n\tDescription: " + description
            response += "\n\tSize: " + size + "\n\tTotal Length: " + length + "ns\n\t24-Hour Growth: " + str(growth)
            response += "\n\tRuns-Clones-Gens: " + runs + '-' + clones + '-' + gens + '\n\tCredit: ' + credit + '\n\tTimeout/Deadline: '
            response += timeout + '/' + deadline
                    
        except:
            response = "Error: Results file does not exist for project: " + project_number
            pass
        
    else:
        response = "Invalid Project Number: " + project_number + ".  No data directory found."
    
    return response

def get_user_info(user):
    
    sizes, lengths, growths = [],[],[]
    
    project_data = sorted(glob.glob(results_prefix + user + '/*_log.txt'))
    filenames = [ re.sub(r'.*%s/'%user, '', i) for i in project_data ]
    project_numbers = [ re.sub(r'_log.txt', '', i) for i in filenames ]
    response = "*Summary for user " + user + ":* \n" + str(len(project_numbers)) + " Projects:"
    
    for i in range(len(project_numbers)):
        try:
            with open(project_data[i], 'r') as file:
                lines = [line.strip().split() for line in file.readlines()]
            file.close()
            try:
                sizes.append(lines[-1][2])
                lengths.append(lines[-1][1])
            except:
                pass
            try:
                growths.append(float(lines[-1][1]) - float(lines[-2][1]))

            except IndexError:
                growths.append(0)
        
            response += "\n\n\t*" + project_numbers[i] + "*\n\tSize: " + sizes[-1] + "\n\tLength: "
            response += str(lengths[-1]) + "ns\n\t24-Hour Growth: " + str(growths[-1]) + "ns"

        except IOError:
            pass
    
    try:
        for i in reversed(range(len(sizes))):
            if 'M' in sizes[i] or 'K' in sizes[i]:
                sizes.remove(sizes[i])

        standardized_sizes = [ float(re.sub('G','',i)) for i in sizes ]
        response += "\n\n*Total size:* " + str(sum(standardized_sizes))
        response += "GB\n*Total Length:* " + str(sum([float(i) for i in lengths]))
        response += "ns\n*Total 24-hour Growth:* " + str(sum([float(i) for i in growths])) + "ns"

    except NameError:
        pass
    
    return response

def get_data_usage():
    
    response = "\t*%s Disk Usage:*\n" %hostname
    header = subprocess.check_output('df -h | head -n 1', shell=True).decode("utf-8")
    array0_cmd = subprocess.check_output('df -h | grep array0', shell=True)
    array1_cmd = subprocess.check_output('df -h | grep array1', shell=True)
    
    for line in [header, array0_cmd, array1_cmd]:
        response += line + '\n'
    
    return response    
    

def change_stats(project_number, change):

    if project_number in project_numbers:
        project_xml = project_config_prefix + 'p' + project_number + "/project.xml"

        try: # parse xml for stuff 
            with open(project_xml, 'r') as file:
                lines = [line.strip().split() for line in file.readlines()]
            file.close()

            for j in range(len(lines)):
                if '<stats-credit' in lines[j] or '<stats_credit' in lines[j]:
                   credit = lines[j][1][3:-3]
        except:
            response = "Error: No project.xml file found for project: " + project_number
            pass

        if '%' in change:
            try:
                adjustment = 1 + float(re.sub('%','',change))/100
                new_credit = int(float(credit) * adjustment)

            except:
                response = "Invalid stats adjustment argument: " + change
                pass

            try:
                slack_client.api_call(
                  "chat.postMessage",
                  channel=channel,
                  text = "Are you sure you want to change stats for project %s by %s from %s to %s?\nTo reply, type either:\n\t`@%s yes`\n\t`@%s no`"
                    %(project_number, change, credit, str(new_credit), hostname, hostname))

                choice = wait_for_response()

                if choice == 'y':
                    backup = project_config_prefix + 'p' + project_number + '/.' + timestr + '_project.xml'
                    os.rename(project_xml, backup)
                    cmd = 'sed "s/  <stats.*/  <stats_credit v=\\"%s\\"\/>/" %s > %s' %(new_credit, backup, project_xml)
                    print(cmd)
                    subprocess.check_output(cmd, shell=True)
                    response = "Stats changed successfully.\nBackup xml created at:\n\t`%s`\nTo restart server daemon, type:\n\t`@%s restart`" %(backup,hostname)

                elif choice == 'n':
                    response = "Change cancelled. No changes made."

                else:
                    response = "Timeout. Try again."

            except IOError:
                response = "Error: Try Adjusting Manually."
                pass

        else:
            try:
                change = int(change)
                try:
                    slack_client.api_call(
                      "chat.postMessage",
                      channel=channel,
                      text="Are you sure you want to change stats for project %s from %s to %s?\nTo reply, type either:\n\t`@%s yes`\n\t`@%s no`"
                        %(project_number, credit, str(change), hostname, hostname))

                    choice = wait_for_response()

                    if choice == 'y':
                        backup = project_config_prefix + 'p' + project_number + '/.' + timestr + '_project.xml'
                        os.rename(project_xml, backup)  
                        cmd = 'sed "s/  <stats.*/  <stats_credit v=\\"%s\\"\/>/" %s > %s' %(str(change), backup, project_xml)
                        print(cmd)
                        subprocess.check_output(cmd, shell=True)
                        response = "Stats changed successfully.\nBackup xml created at:\n\t`%s`\nTo restart server daemon, type:\n\t`@%s restart`" %(backup,hostname)

                    elif choice == 'n':
                        response = "Change cancelled. No changes made."

                    else:
                        response = "Please respond with either:\n\t`@%s yes`\n\t`@%s no`"
                        
                except IOError:
                    response = "Error: Try Adjusting Manually."
            
            except:
                response = "Invalid arguments. See @" + hostname + " help"
    else:
        response = "Invalid project number: " + project_number

    return response


def wait_for_response():
    timeout, count = 10, 0
    while count <= 20:

        command, channel = parse_bot_commands(slack_client.rtm_read())
        print(command, channel)  
        if command:
            if command.split()[0] in ['yes', 'Yes', 'y', 'Y']:
                return('y')
            elif command.split()[0] in ['no', 'No', 'n', 'N']:
                return('n')
        time.sleep(RTM_READ_DELAY)
        count += 1
    return(0)

def daemon_restart():
    output = subprocess.check_output('sudo /etc/init.d/fah-work -n ' + hostname + ' restart', shell=True)
    return "Server daemon restarting..."

def pick_my_lunch():
    food = ['Subway\'s', 'Burger Tank', 'The Indian Food Truck', 'Mexican Grill Stand', 'Top Bap', 'The Teppanyaki Truck', 'Halal', 'Blaze Pizza', 'Chipotle']
    line = ['You should get something from ', 'Try a tasty treat from ', 'I hear there\'s some yummy stuff at ', 'Maybe you should try ', 'I think you should eat some food from ', 'It\'s been a while since you went to ']
    return  line[np.random.randint(len(line))] + food[np.random.randint(len(food))] + '!'


if __name__ == "__main__":

    if slack_client.rtm_connect(with_team_state=False):
        print("Slack Bot for %s is connected and running!"%hostname)
        # Read bot's user ID by calling Web API method `auth.test`
        starterbot_id = slack_client.api_call("auth.test")["user_id"]

        while True:
            command, channel = parse_bot_commands(slack_client.rtm_read())

            if command:
                handle_command(command, channel)
            time.sleep(RTM_READ_DELAY)

    else:
        print("Connection failed. See Exception Traceback Above.")
