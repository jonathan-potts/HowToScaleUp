###############################################################################
# Name: htsu_sim_path_ex2.py
#
# Purpose: Simulate path through a resource layer with correlation
#
# Usage: To run this, use the following command
#   python htsu_sim_path_ex2.py random_field_100.inp 1.5 2 0.2 100 imgtemp.png > temp.out
#
# The first parameter is random_field_100.inp, which gives the resource layer
# The second parameter (1.5) corresponds to beta_1 from Equation (4) in Supplementary Appendix A
# The third parameter (2) corresponds to kappa from Equation (4) in Supplementary Appendix A
# The fourth parameter (0.2) corresponds to lambda from Equation (4) in Supplementary Appendix A
# The fifth parameter (100) is the number of steps to be simulated
# The sixth parameter (imgtemp.png) is a file storing a plot of the simulated locations, over the resource layer
# The seventh parameter (temp.out) is a file storing the locations of the simulated animal
# 
# Author: Jonathan R. Potts
#
# Feel free to share and adapt, but giving appropriate credit
###############################################################################

import sys, math, random
from matplotlib import pyplot as plt

# File containing layer
curr_arg = 1
infile = open(sys.argv[curr_arg],'r')
# beta_r-value (strength of resource effect)
curr_arg += 1
beta_r = float(sys.argv[curr_arg])
# kappa-value (strength of correlation)
curr_arg += 1
kappa_val = float(sys.argv[curr_arg])
# lambda-value (parameter of step length distribution)
curr_arg += 1
lambda_val = float(sys.argv[curr_arg])
# Number of steps
curr_arg += 1
step_no = int(sys.argv[curr_arg])
# Get file for saving plot
curr_arg += 1
savefile = sys.argv[curr_arg]

# Start the random number generator
random.seed()

# Get the Z-values corresponding to the resource layer, and also exp(beta*Z(x)) 
z_array = [[]]
exp_layer = [[]]
infile.seek(0)
y_val = -1
for line in infile:
    y_val += 1
    if y_val != 0:
        z_array += [[]]
        exp_layer += [[]]
    split_line = line.rsplit()
    for x_val in range(len(split_line)):
        z_array[y_val] += [float(split_line[x_val])]
        exp_layer[y_val] += [math.exp(beta_r*z_array[y_val][x_val])]
infile.close()

# Box size
box_width = len(z_array[0])
box_height = len(z_array)

# Central point
xc = box_width/2
yc = box_height/2

# Start location and bearing
loc_x = [xc]
loc_y = [yc]
alpha_x = 0

# Find locations
for step in range(1,step_no):
    # Draw a random number
    random_no = random.random()
    # Calculate the cumulative probability distribution for the movement kernel
    cum_dist = [[]]
    cum_val = 0
    for y_val in range(len(exp_layer)):
        if y_val != 0:
            cum_dist += [[]]
        for x_val in range(len(exp_layer[y_val])):
            alpha_z = math.atan2(y_val-loc_y[step-1],x_val-loc_x[step-1])
            cum_val += (math.exp(-lambda_val*math.sqrt((float(x_val-loc_x[step-1]))**2+(float(y_val-loc_y[step-1]))**2))*
                        exp_layer[y_val][x_val]*
                        math.exp(kappa_val*math.cos(alpha_x-alpha_z)))
            cum_dist[y_val] += [cum_val]
            alpha_x = alpha_z
    # Normalise probabilities 
    for y_val in range(len(cum_dist)):
        for x_val in range(len(cum_dist[y_val])):
            cum_dist[y_val][x_val] /= cum_val
    
    # Flag that we have not yet found the location given by the random draw
    found = 0
    # Find the x- and y-values corresponding to this random draw
    for y_val in range(len(cum_dist)):
        for x_val in range(len(cum_dist[y_val])):
            if cum_dist[y_val][x_val] > random_no:
                loc_x += [x_val]
                loc_y += [y_val]
                sys.stdout.write("%i\t%i\n" % (x_val, y_val))
                # Flag that we have found the location
                found = 1
                break
        if found == 1:
            break

# Plot resource layer and locations
fig = plt.figure()
fig.set_size_inches(7,6)
fig.add_subplot(1,1,1)
plt.hold('on')
# Plot layer
bottomcontour = -4
topcontour = 4
contourres = 1
filled_contours = plt.contourf(z_array, origin='lower', extent=[0,box_width-1,0,box_height-1], levels=[float(x)/contourres for x in range(bottomcontour,topcontour)],cmap=plt.cm.Greens)
plt.ylabel('Northing',fontsize=20)
plt.xlabel('Easting',fontsize=20)
# Make a colorbar 
cbar = plt.colorbar(filled_contours)
# Plot locations as black dots
plt.scatter(loc_x, loc_y, s=3, c='k', marker='o', edgecolors=None)

# Save figure
plt.savefig(savefile)
plt.show()