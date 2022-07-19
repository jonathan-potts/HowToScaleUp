###############################################################################
# Name: htsu_me_ex1.py
#
# Purpose: Master equation solution for an individual moving through a resource layer
#
# Usage: To run this, use the following command
#   python htsu_me_ex1.py random_field_50.inp 1.5 0.2 100 imgtemp.png 
#
# The first parameter is random_field_50.inp, which gives the resource layer
# The second parameter (1.5) corresponds to beta_R from Equation (S.2) in Supplementary Appendix A
# The third parameter (0.2) corresponds to lambda from Equation (S.2) in Supplementary Appendix A
# The fourth parameter (100) is the number of steps to be simulated
# The fifth parameter (imgtemp.png) is a file storing a plot of the simulated locations, over the resource layer
# 
# Author: Jonathan R. Potts
#
# Feel free to share and adapt, but giving appropriate credit
###############################################################################

import sys, math
from matplotlib import pyplot as plt

# File containing layer
curr_arg = 1
infile = open(sys.argv[curr_arg],'r')
# beta_r-value (strength of resource effect)
curr_arg += 1
beta_r = float(sys.argv[curr_arg])
# lambda-value (parameter of step length distribution)
curr_arg += 1
lambda_val = float(sys.argv[curr_arg])
# Number of steps
curr_arg += 1
step_no = int(sys.argv[curr_arg])
# Get file for saving plot
curr_arg += 1
savefile = sys.argv[curr_arg]

# Box size
box_width = 50
box_height = 50

# Get the Z-values corresponding to the resource layer, exp(beta*Z(x)), the initial conditions for 
# the utilisaton distribution (ud_array; U(s,t) from Eqn (5)), and set up an array for the 
# sum in the master equation (sumval)
z_array = [[]]
exp_layer = [[]]
ud_array = [[]]
sumval = [[]]
infile.seek(0)
y_val = -1
for line in infile:
    y_val += 1
    if y_val != 0:
        z_array += [[]]
        exp_layer += [[]]
        ud_array += [[]]
        sumval += [[]]
    split_line = line.rsplit()
    for x_val in range(len(split_line)):
        z_array[y_val] += [float(split_line[x_val])]
        exp_layer[y_val] += [math.exp(beta_r*z_array[y_val][x_val])]
        # Initial condition at middle of box
        if x_val == box_width/2 and y_val == box_height/2:
            ud_array[y_val] += [1]
        else:
            ud_array[y_val] += [0]
        sumval[y_val] += [0]
infile.close()

# Calculate the movement kernel for each possible step from s'=(from_x,from_y) to s=(to_x,to_y).  This
# is P(s|s') from Equation (5)
kernel = [[[[]]]]
for from_y in range(box_height):
    if from_y != 0:
        kernel += [[[[]]]]
    for from_x in range(box_width):
        if from_x != 0:
            kernel[from_y] += [[[]]]
        kernel_sum = 0
        for to_y in range(box_height):
            if to_y != 0:
                kernel[from_y][from_x] += [[]]
            for to_x in range(box_width):
                kernel[from_y][from_x][to_y] += [math.exp(-lambda_val*math.sqrt((float(from_x-to_x))**2+
                                                 (float(from_y-to_y))**2))*exp_layer[to_y][to_x]]
                # For calculating integral over z of k(z|x)                
                kernel_sum += kernel[from_y][from_x][to_y][to_x]
        # Normalise probabilities
        for to_y in range(box_height):
            for to_x in range(box_width):
                kernel[from_y][from_x][to_y][to_x] /= kernel_sum

for step in range(1,step_no):
    for to_y in range(box_height):
        for to_x in range(box_width):
            sumval[to_y][to_x] = 0
            for from_y in range(box_height):
                for from_x in range(box_width):
                    # P(s|s',t)U(s',t) is the object being summed in Equation (5) 
                    sumval[to_y][to_x] += kernel[from_y][from_x][to_y][to_x]*ud_array[from_y][from_x]
    # Update the probability distribution array
    for y_val in range(box_height):
        for x_val in range(box_width):
            ud_array[y_val][x_val] = sumval[y_val][x_val]

# Do the plot
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
# Plot the utilisation distribution after step_no steps
pd_levels = [0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.1]
plt.contour(ud_array, origin='lower',colors='k', levels=pd_levels,linewidths=2)

# Save figure
plt.savefig(savefile)