# figures with rearrangements, scenarios from dcj2hp


import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

import sys,random

parameters = sys.argv[1:]

if len(parameters) < 3 or len(parameters) > 4:
    print("marker_file source_file scenario_file [centromere]")
    exit()
if len(parameters) == 4:
    centromere = int(parameters[3])
else:
    centromere = 0

file_markers = open(parameters[0],"r").readlines()
file_ancestor = open(parameters[1],"r").readline().strip()
file_scenarios = open(parameters[2],"r")

# markers
strata = {0:0}
for line in file_markers[1:]:
    words = line.strip().split()
    marker = words[20]
    strata[int(marker)] = words[21]
strata[len(file_markers)] = 0 # PAR
strata[centromere] = 0 # centromere

diverse_colors = {0:"yellow","S1a":"#ae017e","S3":"#000000","S1b":"#fbb4b9", "S1c":"#f768a1", "S2":"#0570b0", "DF":"blue"}

# read ancestor
permutation = []
#position = {}
words = file_ancestor.split()
for w in range(len(words)):
    permutation.append(int(words[w]))
    #position[abs(int(words[w]))] = w*abs(int(words[w]))/(int(words[w]))


# read scenario : output a list of couples of breakpoints (couple of consecutive markers)    
line = file_scenarios.readline()
stop = False
scenario = []
circular = True
while len(line) > 0 and not stop:
    #print(line)
    if len(line.split()) > 12 and line.split()[11] == "circular:":
        if line.split()[12] == "0.0":
            circular = False
        else:
            circular = True
        #print("circularity",line,circular)
    if line.strip() == "Path 0" and not circular: # one HP scenario
        #print("debut_scenario")
        scenario = []
        nb_both = 0
        nb_one = 0
        premier = False
        line = file_scenarios.readline()
        pos = 1
        while len(line.strip()) > 0:
            words = line.strip()[1:-1].split("|")
            bp1 = words[0].split(",")
            bp1[0] = int((int(bp1[0])+1)/2)
            bp1[1] = int((int(bp1[1])+1)/2)
            bp2 = words[1].split(",")
            bp2[0] = int((int(bp2[0])+1)/2)
            bp2[1] = int((int(bp2[1])+1)/2)
            if strata[bp1[0]] != strata[bp1[1]]:
                if strata[bp2[0]] != strata[bp2[1]]:
                    nb_both = nb_both + 1
                    if pos == 1:
                        premier = True
                else:
                    nb_one = nb_one + 1
            elif strata[bp2[0]] != strata[bp2[1]]:
                nb_one = nb_one + 1
            scenario.append([bp1,bp2])
            line = file_scenarios.readline()
            pos = pos + 1
 
        if random.random() < 0.01:
            stop = True
 
    line = file_scenarios.readline()
    


# Create a figure and axis
fig, ax = plt.subplots()


# Set limits for the plot
ax.set_xlim(0, len(list(strata.keys()))+2)
ax.set_ylim(-(len(scenario)*2+2),3)

def plot_genome(permutation,position):
    ax.axhline(y=-2*position-0.17, color='black', linewidth=1)
    for e in range(len(permutation)):
        element = permutation[e]
        direction = abs(element) / element
        number = abs(element)
        if number == centromere:
            arrow = plt.Circle(( e+1.25, -(2*position)-0.15 ), 0.5 ,color="yellow")
        elif direction == -1:
            arrow = plt.Arrow(e+1.5, -(2*position+0.3), -0.5, 0, width=1.5, color=diverse_colors[strata[number]])
        else:
            arrow = plt.Arrow(e+1, -(2*position), 0.5, 0, width=1.5, color=diverse_colors[strata[number]])
        ax.add_patch(arrow)

def apply_inversion(permutation,inversion,position):
    bp1 = inversion[0]
    bp1.sort()
    bp2 = inversion[1]
    bp2.sort()
    if bp1[0] == 0:
        if bp1[1] == abs(permutation[0]):
            first_breakpoint = -1
        elif bp1[1] == abs(permutation[-1]):
            first_breakpoint = len(permutation) - 1
        else:
            print("error",bp1,permutation)
    else:
        i = 0
        while abs(permutation[i]) != bp1[0] and abs(permutation[i]) != bp1[1]:
            i = i + 1
        first_breakpoint = i
    if bp2[0] == 0:
        if bp2[1] == abs(permutation[0]):
            second_breakpoint = -1
        elif bp2[1] == abs(permutation[-1]):
            second_breakpoint = len(permutation) - 1
        else:
            print("error",bp2,permutation)
    else:
        i = 0
        while abs(permutation[i]) != bp2[0] and abs(permutation[i]) != bp2[1]:
            i = i + 1
        second_breakpoint = i
    #print(first_breakpoint,second_breakpoint,permutation[first_breakpoint],permutation[second_breakpoint])
    ax.plot([first_breakpoint+2, second_breakpoint+1.5], [-(2*position-1), -(2*position-1)], color='green', linewidth=2)
    result = []
    temp = first_breakpoint
    first_breakpoint = min(first_breakpoint,second_breakpoint)
    second_breakpoint = max(temp,second_breakpoint)
    #print(first_breakpoint,second_breakpoint)
    i = 0
    while i <= first_breakpoint:
        result.append(permutation[i])
        i = i + 1
    i = second_breakpoint
    while i > first_breakpoint:
        result.append(-permutation[i])
        i = i - 1
    i = second_breakpoint + 1
    while i < len(permutation):
        result.append(permutation[i])
        i = i + 1
    #print("result",result)
    return result

# Draw an arrow from (0, 0) to (2, 2)
#for marker in position:
plot_genome(permutation,0)
for s in range(len(scenario)):
    permutation = apply_inversion(permutation,scenario[s],s+1)
    #print(permutation)
    plot_genome(permutation,s+1)
    

# Add labels
#plt.text(0.5, 1, "Arrow Example", fontsize=12, ha='center')
#plt.text(2, 2.2, "(2, 2)", fontsize=10, ha='center')
#plt.text(-0.2, -0.2, "(0, 0)", fontsize=10, ha='center')

# Set axis labels
plt.xlabel("X")
plt.ylabel("Y")

# Show the plot
plt.show()

