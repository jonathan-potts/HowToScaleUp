###############################################################################
# Name: htsu_fig1.py
#
# Purpose: Plots Figure 1 of "How to scale up from animal movement decisions to spatio-temporal
#          patterns: an approach via step selection" by JR Potts and L Borger
#
# Usage: python htsu_fig1.py random_field_100.inp 2 0.25 0.2 1 1000 -4 4 1 fig1.png
#
# Author: Jonathan R. Potts
#
# Feel free to share and adapt, but giving appropriate credit
###############################################################################

import sys, math, random, numpy
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
# Von-mises parameter for turning angles
curr_arg += 1
vm = float(sys.argv[curr_arg])
# Von-mises parameter for turning angles
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
# Start the random number generator
random.seed()

# Get the Z-values and their exponent
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
start_x = xc-5
start_y = yc-5
alpha_x = -math.pi/4

# Find distribution of next location
prob_dist = [[]]
cum_val = 0
for y_val in range(len(exp_layer)):
    if y_val != 0:
        prob_dist += [[]]
    for x_val in range(len(exp_layer[y_val])):
        alpha_z = math.atan2(y_val-start_y,x_val-start_x)
        prob = (math.exp(-lambda_val*math.sqrt((float(x_val-start_x))**2+(float(y_val-start_y))**2))*
                exp_layer[y_val][x_val]*
                math.exp(-beta_c*math.sqrt((float(x_val-xc))**2+(float(y_val-yc))**2))*
                math.exp(vm*math.cos(alpha_x-alpha_z)))
        prob_dist[y_val] += [prob]
        cum_val += prob
# Normalise probabilities
mk_dist = [[]]
for y_val in range(len(prob_dist)):
    if y_val != 0:
        mk_dist += [[]]
    for x_val in range(len(prob_dist[y_val])):
        mk_dist[y_val] += [prob_dist[y_val][x_val]/cum_val]

# Find locations
loc_x = [start_x]
loc_y = [start_y]
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
            alpha_x = alpha_z
            alpha_z = math.atan2(y_val-loc_y[step-1],x_val-loc_x[step-1])
            cum_val += (math.exp(-lambda_val*math.sqrt((float(x_val-loc_x[step-1]))**2+(float(y_val-loc_y[step-1]))**2))*
                    exp_layer[y_val][x_val]*
                    math.exp(-beta_c*math.sqrt((float(x_val-xc))**2+(float(y_val-yc))**2))*
                    math.exp(vm*math.cos(alpha_x-alpha_z)))
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

# Calculate the smoothing parameters (h)
h_val = math.sqrt((numpy.var(loc_x)+numpy.var(loc_y))/2)*(float(len(loc_x))**(-1.0/6.0))

# Calculate the KDE
unit = 1
lat_space = 1
kde =[[]]
for y_val in range(0,box_height):
    if y_val != 0:
        kde += [[]]
    for x_val in range(0,box_width):
        kde[y_val] += [0]
        for data in range(len(loc_x)):
            kde[y_val][x_val] += (1/(unit**2))*math.exp(-((x_val*lat_space-loc_x[data])**2+(y_val*lat_space-loc_y[data])**2)/(2*(h_val**2)))/(2*math.pi*len(loc_x)*(h_val**2))

# Plot MK and UD
fig = plt.figure()
fig.set_size_inches(12,6)
fig.add_subplot(1,2,1)
plt.hold('on')
# Plot resource layer
filled_contours = plt.contourf(z_array, origin='lower', extent=[0,box_width*lat_space*unit,0,box_height*lat_space*unit], levels=[float(x)/contourres for x in range(bottomcontour,topcontour)],cmap=plt.cm.Greens)
plt.ylabel('Northing',fontsize=20)
plt.xlabel('Easting',fontsize=20)
# Plot MK contour
pc1 = 50.0 # kernel percentage 
pc2 = 95.0 # kernel percentage 
max_found = 2
curr_pc = 0
# Set up array for storing contour locations.  This is an array of zeros initially.  Later,
# the places within the pc% kernel of the UD will be populated with ones
hr_contour = [[]]
for y_val in range(len(mk_dist)):
    if y_val != 0:
        hr_contour += [[]]
    for x_val in range(len(mk_dist[0])):
        hr_contour[y_val] += [0]
# Populate the places within the pc% kernel of the UD with ones
while curr_pc < pc2/100.0:
    curr_max = 0
    max_x = 0
    max_y = 0
    # Find the highest value in the array that is less than max_found.  This then becomes a 
    for y_val in range(len(mk_dist)):
        for x_val in range(len(mk_dist[0])):
            if (mk_dist[y_val][x_val] > curr_max) and (mk_dist[y_val][x_val] < max_found):
                curr_max = mk_dist[y_val][x_val]
                max_x = x_val
                max_y = y_val
    # Set the output array at position where this max is found to 1, 2, or 3
    if curr_pc < pc1/100:
        hr_contour[max_y][max_x] = 2
    else:
        hr_contour[max_y][max_x] = 1
    max_found = curr_max
    curr_pc += max_found
plt.axis([20,80,30,80])
plt.contour(hr_contour, origin='lower',colors=['b','m'], extent=[0,len(mk_dist[0]),0,len(mk_dist)],levels=[0,1,2],linestyles='solid',linewidths=2)
plt.scatter([xc],[yc], s=70, c='k', marker='o', edgecolors=None)
plt.scatter([start_y],[start_y], s=70, c='k', marker='o', edgecolors=None)
plt.plot([start_x-10,start_x-0.5],[start_y+10,start_y+0.5], linewidth=2,color='k')
plt.plot([start_x-10,start_x-10, start_x-10, start_x-10, start_x-10,start_x-10],[start_y+10,start_y+12,start_y+14,start_y+16,start_y+18,start_y+20], 'k--',linewidth=1)
plt.plot([start_x-0.5,start_x-0.5],[start_y+0.5,start_y+2.5], linewidth=2,color='k')
plt.plot([start_x-2.5,start_x-0.5],[start_y+0.5,start_y+0.5], linewidth=2,color='k')
plt.text(21,77,'a)',fontsize=26)
plt.text(48.5,50.8,r'${\bf x}_C$',fontsize=26)
plt.text(46,44.2,r'${\bf x}$',fontsize=26)
plt.text(37.1,54.0,r'$\alpha_{\bf x}$',fontsize=26)

fig.add_subplot(1,2,2)
plt.hold('on')
# Plot resource layer
filled_contours = plt.contourf(z_array, origin='lower', extent=[0,box_width*lat_space*unit,0,box_height*lat_space*unit], levels=[float(x)/contourres for x in range(bottomcontour,topcontour)],cmap=plt.cm.Greens)
plt.ylabel('Northing',fontsize=20)
plt.xlabel('Easting',fontsize=20)
# Plot MK contour
pc1 = 50.0 # kernel percentage 
pc2 = 95.0 # kernel percentage 
max_found = 2
curr_pc = 0
# Set up array for storing contour locations.  This is an array of zeros initially.  Later,
# the places within the pc% kernel of the UD will be populated with ones
hr_contour = [[]]
for y_val in range(len(mk_dist)):
    if y_val != 0:
        hr_contour += [[]]
    for x_val in range(len(mk_dist[0])):
        hr_contour[y_val] += [0]
# Populate the places within the pc% kernel of the UD with ones
while curr_pc < pc2/100.0:
    curr_max = 0
    max_x = 0
    max_y = 0
    # Find the highest value in the array that is less than max_found.  This then becomes a 
    for y_val in range(len(kde)):
        for x_val in range(len(kde[0])):
            if (kde[y_val][x_val] > curr_max) and (kde[y_val][x_val] < max_found):
                curr_max = kde[y_val][x_val]
                max_x = x_val
                max_y = y_val
    # Set the output array at position where this max is found to 1, 2, or 3
    if curr_pc < pc1/100:
        hr_contour[max_y][max_x] = 2
    else:
        hr_contour[max_y][max_x] = 1
    max_found = curr_max
    curr_pc += max_found
plt.axis([20,80,30,80])
plt.contour(hr_contour, origin='lower',colors=['b','m'], extent=[0,len(kde[0]),0,len(kde)],levels=[0,1,2],linestyles='solid',linewidths=2)
plt.scatter(loc_x, loc_y, s=3, c='k', marker='o', edgecolors=None)
plt.text(21,77,'b)',fontsize=26)

# Save and show figure
plt.savefig(savefile)
plt.show()