###############################################################################
# Name: htsu_fig2.py
#
# Purpose: Plots Figure 2 of "How to scale up from animal movement decisions to spatio-temporal
#          patterns: an approach via step selection" by JR Potts and L Borger
#
# Usage: python htsu_fig2.py random_field_100.inp 1.5 0.2 0.2 1000 -4 4 1 fig2.png 
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
# beta_c-value (strength of central point attraction)
curr_arg += 1
beta_c = float(sys.argv[curr_arg])
# lambda-value (parameter of step length distribution)
curr_arg += 1
lambda_val = float(sys.argv[curr_arg])
# Number of steps
curr_arg += 1
step_no = int(sys.argv[curr_arg])
# Get value of topcontour * contourres and bottomcontour*contourres
curr_arg += 1
bottomcontour = int(sys.argv[curr_arg])
curr_arg += 1
topcontour = int(sys.argv[curr_arg])
curr_arg += 1
contourres = int(sys.argv[curr_arg])
# Get file for saving plot
curr_arg += 1
savefile = sys.argv[curr_arg]

############################
# Simulate the path
############################

# Start the random number generator
random.seed()

# Get the Z-values and also exp(beta*Z(x)) and exp(2*beta*Z(x))
z_array = [[]]
exp_layer = [[]]
exp_2layer = [[]]
infile.seek(0)
y_val = -1
for line in infile:
    y_val += 1
    if y_val != 0:
        z_array += [[]]
        exp_layer += [[]]
        exp_2layer += [[]]
    split_line = line.rsplit()
    for x_val in range(len(split_line)):
        z_array[y_val] += [float(split_line[x_val])]
        exp_layer[y_val] += [math.exp(beta_r*z_array[y_val][x_val])]
        exp_2layer[y_val] += [math.exp(beta_r*2*z_array[y_val][x_val])]
infile.close()

# Box size
box_width = len(z_array[0])
box_height = len(z_array)

# Central point
xc = box_width/2
yc = box_height/2

# Start location
loc_x = [xc]
loc_y = [yc]

# Find locations
for step in range(1,step_no):
    # Draw a random number
    random_no = random.random()
    # Calculate the cumulative distribution
    cum_dist = [[]]
    cum_val = 0
    for y_val in range(len(exp_layer)):
        if y_val != 0:
            cum_dist += [[]]
        for x_val in range(len(exp_layer[y_val])):
            cum_val += (math.exp(-lambda_val*math.sqrt((float(x_val-loc_x[step-1]))**2+(float(y_val-loc_y[step-1]))**2))*
                        exp_layer[y_val][x_val]*
                        math.exp(-beta_c*math.sqrt((float(x_val-xc))**2+(float(y_val-yc))**2)))
            cum_dist[y_val] += [cum_val]
    # Normalise probabilities 
    for y_val in range(len(cum_dist)):
        for x_val in range(len(cum_dist[y_val])):
            cum_dist[y_val][x_val] /= cum_val
    
    # Flag that we have not yet found the position
    found = 0
    # Find the x- and y-values corresponding to this random draw
    for y_val in range(len(cum_dist)):
        for x_val in range(len(cum_dist[y_val])):
            if cum_dist[y_val][x_val] > random_no:
                loc_x += [x_val]
                loc_y += [y_val]
                # Flag that we have found the position
                found = 1
                break
        if found == 1:
            break

# Add random jitter to account for the fact that locations may be at any point within a pixel
for point in range(len(loc_x)):
    loc_x[point] += random.random()-0.5
    loc_y[point] += random.random()-0.5

###############################
# Calculate estimated UDs 
###############################

# Estimations from Equation (18) and (19) 
ud1 = [[]]
ud2 = [[]]
ud1_tot = 0
ud2_tot = 0
for y_val in range(len(exp_layer)):
    if y_val != 0:
        ud1 += [[]]
        ud2 += [[]]
    for x_val in range(len(exp_layer[0])):
        ud1[y_val] += [exp_layer[y_val][x_val]*math.exp(-beta_c*math.sqrt((float(x_val-xc))**2+(float(y_val-yc))**2))]
        ud2[y_val] += [exp_2layer[y_val][x_val]*math.exp(-2*beta_c*math.sqrt((float(x_val-xc))**2+(float(y_val-yc))**2))]
        ud1_tot += ud1[y_val][x_val] 
        ud2_tot += ud2[y_val][x_val]
for y_val in range(len(exp_layer)):
    for x_val in range(len(exp_layer[0])):
        ud1[y_val][x_val] /= ud1_tot
        ud2[y_val][x_val] /= ud2_tot

# Estimations from Equation (17) requires a numerical integration (i.e. a sum)
ud = [[]]
ud_tot = 0 # Sum of the UD values
for y_val in range(len(exp_layer)):
    if y_val != 0:
      ud += [[]]
    for x_val in range(len(exp_layer[0])):
        ud[y_val] += [0]
        # Calculate the integral of phi(|x-z|)*exp(beta.Z(z)) over all z.  Here x=(y_val,x_val), z=(int_x,int_y)
        for int_y in range(len(exp_layer)):
            for int_x in range(len(exp_layer[0])):
                ud[y_val][x_val] += (math.exp(-lambda_val*math.sqrt((float(x_val-int_x))**2+(float(y_val-int_y))**2))*
                                     exp_layer[int_y][int_x]*
                                     math.exp(-beta_c*math.sqrt((float(int_x-xc))**2+(float(int_y-yc))**2)))
        # Multiply by exp(beta.Z(x))
        ud[y_val][x_val]*=exp_layer[y_val][x_val]*math.exp(-beta_c*math.sqrt((float(x_val-xc))**2+(float(y_val-yc))**2))
        # Increment sum of UD values
        ud_tot += ud[y_val][x_val]
# Normalise
for y_val in range(len(exp_layer)):
    for x_val in range(len(exp_layer[0])):
        ud[y_val][x_val] /= ud_tot

###############################
# Plot locations and contours
###############################
        
# Plot resource layer and locations
fig = plt.figure()
fig.set_size_inches(12,12)
fig.add_subplot(2,2,1)
plt.hold('on')
filled_contours = plt.contourf(z_array, origin='lower', extent=[0,box_width-1,0,box_height-1], levels=[float(x)/contourres for x in range(bottomcontour,topcontour)],cmap=plt.cm.Greens)
plt.ylabel('Northing',fontsize=20)
# Make a colorbar 
#cbar = plt.colorbar(filled_contours)

# Plot locations as black dots
plt.scatter(loc_x, loc_y, s=3, c='k', marker='o', edgecolors=None)
plt.text(2,90,'a)',fontsize=26)

# Plot estimated UD and locations: Equation (17)
fig.add_subplot(2,2,2)
plt.hold('on')
# Plot resource layer
filled_contours = plt.contourf(z_array, origin='lower', extent=[0,box_width-1,0,box_height-1], levels=[float(x)/contourres for x in range(bottomcontour,topcontour)],cmap=plt.cm.Greens)
# Plot UD contours
pc = 95.0 # kernel percentage 
max_found = 2
curr_pc = 0
# Set up array for storing contour locations.  This is an array of zeros initially.  Later,
# the places within the pc% kernel of the UD will be populated with ones
hr_contour = [[]]
for y_val in range(len(ud)):
    if y_val != 0:
        hr_contour += [[]]
    for x_val in range(len(ud[0])):
        hr_contour[y_val] += [0]
# Populate the places within the pc% kernel of the UD with ones
while curr_pc < pc/100.0:
    curr_max = 0
    max_x = 0
    max_y = 0
    # Find the highest value in the array that is less than max_found.  This then becomes a 
    for y_val in range(len(ud)):
        for x_val in range(len(ud[0])):
            if (ud[y_val][x_val] > curr_max) and (ud[y_val][x_val] < max_found):
                curr_max = ud[y_val][x_val]
                max_x = x_val
                max_y = y_val
    # Set the output array at position where this max is found to 1
    hr_contour[max_y][max_x] = 1
    max_found = curr_max
    curr_pc += max_found
plt.contour(hr_contour, origin='lower',colors='b', extent=[0,len(ud[0]),0,len(ud)],levels=[0,1],linestyles='solid')
plt.scatter(loc_x, loc_y, s=3, c='k', marker='o', edgecolors=None)
plt.text(2,90,'b)',fontsize=26)

# Plot estimated UD and locations: Equation (18)
fig.add_subplot(2,2,3)
plt.hold('on')
# Plot resource layer
filled_contours = plt.contourf(z_array, origin='lower', extent=[0,box_width-1,0,box_height-1], levels=[float(x)/contourres for x in range(bottomcontour,topcontour)],cmap=plt.cm.Greens)
# Plot UD contours
pc = 95.0 # kernel percentage 
max_found = 2
curr_pc = 0
# Set up array for storing contour locations.  This is an array of zeros initially.  Later,
# the places within the pc% kernel of the UD will be populated with ones
for y_val in range(len(hr_contour)):
    for x_val in range(len(hr_contour[0])):
        hr_contour[y_val][x_val] = 0
# Populate the places within the pc% kernel of the UD with ones
while curr_pc < pc/100.0:
    curr_max = 0
    max_x = 0
    max_y = 0
    # Find the highest value in the array that is less than max_found.  This then becomes a 
    for y_val in range(len(ud1)):
        for x_val in range(len(ud1[0])):
            if (ud1[y_val][x_val] > curr_max) and (ud1[y_val][x_val] < max_found):
                curr_max = ud1[y_val][x_val]
                max_x = x_val
                max_y = y_val
    # Set the output array at position where this max is found to 1
    hr_contour[max_y][max_x] = 1
    max_found = curr_max
    curr_pc += max_found
plt.contour(hr_contour, origin='lower',colors='b', extent=[0,len(ud1[0]),0,len(ud1)],levels=[0,1],linestyles='solid')
plt.scatter(loc_x, loc_y, s=3, c='k', marker='o', edgecolors=None)
plt.xlabel('Easting',fontsize=20)
plt.ylabel('Northing',fontsize=20)
plt.text(2,90,'c)',fontsize=26)

# Plot estimated UD and locations: Equation (19)
fig.add_subplot(2,2,4)
plt.hold('on')
# Plot resource layer
filled_contours = plt.contourf(z_array, origin='lower', extent=[0,box_width-1,0,box_height-1], levels=[float(x)/contourres for x in range(bottomcontour,topcontour)],cmap=plt.cm.Greens)
# Plot UD contours
pc = 95.0 # kernel percentage 
max_found = 2
curr_pc = 0
# Set up array for storing contour locations.  This is an array of zeros initially.  Later,
# the places within the pc% kernel of the UD will be populated with ones
for y_val in range(len(hr_contour)):
    for x_val in range(len(hr_contour[0])):
        hr_contour[y_val][x_val] = 0
# Populate the places within the pc% kernel of the UD with ones
while curr_pc < pc/100.0:
    curr_max = 0
    max_x = 0
    max_y = 0
    # Find the highest value in the array that is less than max_found.  This then becomes a 
    for y_val in range(len(ud2)):
        for x_val in range(len(ud2[0])):
            if (ud2[y_val][x_val] > curr_max) and (ud2[y_val][x_val] < max_found):
                curr_max = ud2[y_val][x_val]
                max_x = x_val
                max_y = y_val
    # Set the output array at position where this max is found to 1
    hr_contour[max_y][max_x] = 1
    max_found = curr_max
    curr_pc += max_found
plt.contour(hr_contour, origin='lower',colors='b', extent=[0,len(ud2[0]),0,len(ud2)],levels=[0,1],linestyles='solid')
plt.scatter(loc_x, loc_y, s=3, c='k', marker='o', edgecolors=None)
plt.xlabel('Easting',fontsize=20)
plt.text(2,90,'d)',fontsize=26)

# Save and show figure
plt.savefig(savefile)
plt.show()