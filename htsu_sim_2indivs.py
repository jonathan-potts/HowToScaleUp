###############################################################################
# Name: htsu_sim_2indivs.py
#
# Purpose: Simulate 2 interacting individuals moving through a resource layer with mutual avoidance
#
# Usage: python htsu_sim_2indivs.py random_field_100.inp 0.5 1000 20000 15000 10000 imgtemp.png > htsu_sim_2indivs.out
###############################################################################

import sys, math, random
from matplotlib import pyplot as plt

# File containing layer
curr_arg = 1
infile = open(sys.argv[curr_arg],'r')
# beta_r-value (strength of resource effect)
curr_arg += 1
beta_r = float(sys.argv[curr_arg])
# beta-value relating to strength of avoidance
curr_arg += 1
beta_avoid = float(sys.argv[curr_arg])
# Number of steps
curr_arg += 1
step_no = int(sys.argv[curr_arg])
# Burn in no of steps
curr_arg += 1
burn_in = int(sys.argv[curr_arg])
# Number of steps that OD is calculated over
curr_arg += 1
no_steps_for_od = int(sys.argv[curr_arg])
# Get file for saving plot
curr_arg += 1
savefile = sys.argv[curr_arg]

# Amount to decrease OD by each step
dec_inc = 1/float(no_steps_for_od)

# Start the random number generator
random.seed()

# Get the R-values and also exp(beta_R*R(x))
r_array = [[]]
exp_layer = [[]]
infile.seek(0)
y_val = -1
for line in infile:
    y_val += 1
    if y_val != 0:
        r_array += [[]]
        exp_layer += [[]]
    split_line = line.rsplit()
    for x_val in range(len(split_line)):
        r_array[y_val] += [float(split_line[x_val])]
        exp_layer[y_val] += [math.exp(beta_r*r_array[y_val][x_val])]
infile.close()

# Box size
box_width = len(r_array[0])
box_height = len(r_array)

# Set up array to store occurrence distributions
no_indivs = 2
occurrence_dist = [[[]]]
for individual in range(no_indivs):
  if individual != 0:
    occurrence_dist += [[[]]]
  # Loop through space
  for spacey in range(box_height):
    if spacey != 0:
      occurrence_dist[individual] += [[]]
    for spacex in range(box_width):
      occurrence_dist[individual][spacey] += [0]  

# Start locations: in the centre of the left- and right-hand halves of the domain
loc_x = [[]]
loc_y = [[]]
xc = [box_width/4,3*box_width/4]
yc = [box_height/2,box_height/2]
for indiv in range(no_indivs):
  if indiv != 0:
    loc_x += [[]]
    loc_y += [[]]
  loc_x[indiv] += [xc[indiv]]
  loc_y[indiv] += [yc[indiv]]
  occurrence_dist[indiv][yc[indiv]][xc[indiv]] = 1

# Simulate paths
for step in range(1,step_no):
  for indiv in range(no_indivs):
    # Draw a random number
    random_no = random.random()
    # Calculate the probabilities of moving up/down/left/right
    if loc_y[indiv][step-1] == 0:
      # Can't move up if at the top of the domain
      up_prob = 0
    else:
      # Include the effect of the resource layer
      up_prob = exp_layer[loc_y[indiv][step-1]-1][loc_x[indiv][step-1]]
      # Include the effect of the occurrence distribution of the other individual
      for indiv2 in range(no_indivs):
        if indiv2 != indiv:
          up_prob *= math.exp(-beta_avoid*occurrence_dist[indiv2][loc_y[indiv][step-1]-1][loc_x[indiv][step-1]])
    if loc_y[indiv][step-1] == box_height - 1:
      # Can't move down if at the bottom of the domain
      down_prob = 0
    else:
      # Include the effect of the resource layer
      down_prob = exp_layer[loc_y[indiv][step-1]+1][loc_x[indiv][step-1]]
      # Include the effect of the occurrence distribution of the other individual
      for indiv2 in range(no_indivs):
        if indiv2 != indiv:
          down_prob *= math.exp(-beta_avoid*occurrence_dist[indiv2][loc_y[indiv][step-1]+1][loc_x[indiv][step-1]])
    if loc_x[indiv][step-1] == 0:
      # Can't move left if at the far left of the domain
      left_prob = 0
    else:
      # Include the effect of the resource layer
      left_prob = exp_layer[loc_y[indiv][step-1]][loc_x[indiv][step-1]-1]
      # Include the effect of the occurrence distribution of the other individual
      for indiv2 in range(no_indivs):
        if indiv2 != indiv:
          left_prob *= math.exp(-beta_avoid*occurrence_dist[indiv2][loc_y[indiv][step-1]][loc_x[indiv][step-1]-1])
    if loc_x[indiv][step-1] == box_width - 1:
      # Can't move right if at the far right of the domain
      right_prob = 0
    else:
      # Include the effect of the resource layer
      right_prob = exp_layer[loc_y[indiv][step-1]][loc_x[indiv][step-1]+1]
      # Include the effect of the occurrence distribution of the other individual
      for indiv2 in range(no_indivs):
        if indiv2 != indiv:
          right_prob *= math.exp(-beta_avoid*occurrence_dist[indiv2][loc_y[indiv][step-1]][loc_x[indiv][step-1]+1])
    tot_prob = up_prob + down_prob + left_prob + right_prob
    # Normalise probabilities 
    up_prob /= tot_prob
    down_prob /= tot_prob
    left_prob /= tot_prob
    right_prob /= tot_prob

    # Move individual
    if random_no < up_prob:
      # Move on up
      loc_x[indiv] += [loc_x[indiv][step-1]]
      loc_y[indiv] += [loc_y[indiv][step-1]-1]
    elif random_no < up_prob + down_prob:
      # Get down
      loc_x[indiv] += [loc_x[indiv][step-1]]
      loc_y[indiv] += [loc_y[indiv][step-1]+1]
    elif random_no < up_prob + down_prob + left_prob:
      # Go left
      loc_x[indiv] += [loc_x[indiv][step-1]-1]
      loc_y[indiv] += [loc_y[indiv][step-1]]
    else:
      # To the right
      loc_x[indiv] += [loc_x[indiv][step-1]+1]
      loc_y[indiv] += [loc_y[indiv][step-1]]

    # Write locations to file
    if step > burn_in:
      if indiv == no_indivs - 1:
        sys.stdout.write("%i\t%i\n" % (loc_x[indiv][step], loc_y[indiv][step]))
      else:
        sys.stdout.write("%i\t%i\t" % (loc_x[indiv][step], loc_y[indiv][step]))
   
  # Update occurrence distributions 
  for indiv in range(no_indivs):
    if step < no_steps_for_od:
      # We are within the first no_steps_for_od steps, so decrease the OD at the first step
      occurrence_dist[indiv][loc_y[indiv][0]][loc_x[indiv][0]] -= dec_inc
    else: 
      # Deprecate the occurrence distribution at the first position used to calculate this distribution
      occurrence_dist[indiv][loc_y[indiv][step - no_steps_for_od]][loc_x[indiv][step - no_steps_for_od]] -= dec_inc
    # Increment the occurrence distribution at the current position
    occurrence_dist[indiv][loc_y[indiv][step]][loc_x[indiv][step]] += dec_inc

# Add random jitter to account for the fact that locations may be at any point within a pixel
for indiv in range(no_indivs):
  for point in range(step_no):
    loc_x[indiv][point] += random.random()-0.5
    loc_y[indiv][point] += random.random()-0.5
    
# Plot resource layer with locations on top
fig = plt.figure()
fig.set_size_inches(6,6)
fig.add_subplot(1,1,1)
plt.hold('on')
bottomcontour = -4
topcontour = 4
contourres = 1
filled_contours = plt.contourf(r_array, origin='lower', extent=[0,box_width-1,0,box_height-1], levels=[float(x)/contourres for x in range(bottomcontour,topcontour)],cmap=plt.cm.Greens)
colors = ['k','r','b','m']
for indiv in range(no_indivs):
  plt.scatter(loc_x[indiv],loc_y[indiv],s=3, c=colors[indiv], marker='o', edgecolors=colors[indiv])
  
# Save and show figure
plt.savefig(savefile)
plt.show()