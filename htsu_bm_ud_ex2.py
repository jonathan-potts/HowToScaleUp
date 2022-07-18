###############################################################################
# Name: htsu_bm_ud_ex2.py
#
# Purpose: An example of calculating the UD from a movement kernel, using the
#          Barnett-Moorcroft formula.  The animal's movement is biased towards 
#          higher quality resources, given by a single resource layer, as well
#          as towards a central point.  Also plots the approximations given in 
#          Equations (12) and (13).
#
# Usage: To run this, use the following command
#   python htsu_bm_ud_ex2.py random_field_100.inp 1.5 0.2 0.2 imgtemp.png 
#
# The first parameter is random_field_100.inp, which gives the resource layer
# The second parameter (1.5) corresponds to beta_R from the paper and Supplementary Appendices
# The third parameter (0.2) corresponds to beta_C 
# The fourth parameter (0.2) corresponds to lambda 
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
# beta_c-value (strength of central point attraction)
curr_arg += 1
beta_c = float(sys.argv[curr_arg])
# lambda-value (parameter of step length distribution)
curr_arg += 1
lambda_val = float(sys.argv[curr_arg])
# Get file for saving plot
curr_arg += 1
savefile = sys.argv[curr_arg]

# Get the R-values (resource layer) and also exp(beta_R*R(x)) and exp(2*beta_R*R(x))
r_array = [[]]
exp_layer = [[]]
exp_2layer = [[]]
infile.seek(0)
y_val = -1
for line in infile:
    y_val += 1
    if y_val != 0:
        r_array += [[]]
        exp_layer += [[]]
        exp_2layer += [[]]
    split_line = line.rsplit()
    for x_val in range(len(split_line)):
        r_array[y_val] += [float(split_line[x_val])]
        exp_layer[y_val] += [math.exp(beta_r*r_array[y_val][x_val])]
        exp_2layer[y_val] += [math.exp(beta_r*2*r_array[y_val][x_val])]
infile.close()

# Box size
box_width = len(r_array[0])
box_height = len(r_array)

# Central point
xc = box_width/2
yc = box_height/2

###############################
# Calculate estimated UDs 
###############################

# Barnett-Moorcroft UD from Equation (11).  This requires an integration.  The integral becomes 
# a sum when translated into computer code.
ud = [[]]
ud_tot = 0 # Sum of the UD values
for y_val in range(box_height):
    if y_val != 0:
      ud += [[]]
    for x_val in range(box_width):
        ud[y_val] += [0]
        # Calculate the integral of phi(|x-z|)*exp(beta_R*R(z)-beta_C*|z-x|) over all z.  Here x=(y_val,x_val), z=(int_x,int_y)
        for int_y in range(box_height):
            for int_x in range(box_width):
                ud[y_val][x_val] += (math.exp(-lambda_val*math.sqrt((float(x_val-int_x))**2+(float(y_val-int_y))**2))*
                                     exp_layer[int_y][int_x]*
                                     math.exp(-beta_c*math.sqrt((float(int_x-xc))**2+(float(int_y-yc))**2)))
        # Multiply by exp(beta_R*R(z)-beta_C*|z-x|)
        ud[y_val][x_val]*=exp_layer[y_val][x_val]*math.exp(-beta_c*math.sqrt((float(x_val-xc))**2+(float(y_val-yc))**2))
        # Increment sum of UD values
        ud_tot += ud[y_val][x_val]
# Normalise
for y_val in range(box_height):
    for x_val in range(box_width):
        ud[y_val][x_val] /= ud_tot

# Estimations from Equation (12) (the limit as phi is arbitrarilty narrow) and (13) (the limit as 
# phi samples from the whole landscape)
ud1 = [[]]
ud2 = [[]]
ud1_tot = 0
ud2_tot = 0
for y_val in range(box_height):
    if y_val != 0:
        ud1 += [[]]
        ud2 += [[]]
    for x_val in range(box_width):
        ud1[y_val] += [exp_2layer[y_val][x_val]*math.exp(-2*beta_c*math.sqrt((float(x_val-xc))**2+(float(y_val-yc))**2))]
        ud2[y_val] += [exp_layer[y_val][x_val]*math.exp(-beta_c*math.sqrt((float(x_val-xc))**2+(float(y_val-yc))**2))]
        ud1_tot += ud1[y_val][x_val] 
        ud2_tot += ud2[y_val][x_val]
for y_val in range(box_height):
    for x_val in range(box_width):
        ud1[y_val][x_val] /= ud1_tot
        ud2[y_val][x_val] /= ud2_tot

###############################
# Plot contours
###############################
        
# Plot resource layer in the top-left panel
fig = plt.figure()
fig.set_size_inches(12,12)
fig.add_subplot(2,2,1)
plt.hold('on')
bottomcontour = -4
topcontour = 4
contourres = 1
filled_contours = plt.contourf(r_array, origin='lower', extent=[0,box_width-1,0,box_height-1], levels=[float(x)/contourres for x in range(bottomcontour,topcontour)],cmap=plt.cm.Greens)

# Plot UD Barnett-Moorcroft UD: Equation (11)
fig.add_subplot(2,2,2)
plt.hold('on')
pd_levels = [0.0001,0.001,0.01,0.1]
plt.contour(ud, origin='lower',colors='k',levels=pd_levels, extent=[0,len(ud[0]),0,len(ud)])
plt.text(2,90,'b)',fontsize=26)

# Plot estimated UD in narrow-kernel limit: Equation (12)
fig.add_subplot(2,2,3)
plt.hold('on')
plt.contour(ud1, origin='lower',colors='k',levels=pd_levels, extent=[0,len(ud1[0]),0,len(ud1)])
plt.text(2,90,'c)',fontsize=26)

# Plot estimated UD in wide-kernel limit: Equation (13)
fig.add_subplot(2,2,4)
plt.hold('on')
plt.contour(ud2, origin='lower',colors='k',levels=pd_levels, extent=[0,len(ud2[0]),0,len(ud2)])
plt.text(2,90,'d)',fontsize=26)

# Save and show figure
plt.savefig(savefile)
plt.show()