################################################################################
# Name: htsu_fig3.py
#
# Purpose: Plots Figure 3 of "How to scale up from animal movement decisions to spatio-temporal
#          patterns: an approach via step selection" by JR Potts and L Borger
#
# Usage: python htsu_fig3.py htsu_sim_4paths_full.out htsu_sim_4paths_noavoid.out htsu_sim_4paths_nocp.out htsu_sim_4paths_nores.out random_field_100.inp 0 100 0 100 1 1 0 10 fig3.png
################################################################################

import sys,math,pylab,numpy,random

# Get filename of files containing animal positions and open for reading
curr_arg = 1
pos_file1 = open(sys.argv[curr_arg],'r')
curr_arg += 1
pos_file2 = open(sys.argv[curr_arg],'r')
curr_arg += 1
pos_file3 = open(sys.argv[curr_arg],'r')
curr_arg += 1
pos_file4 = open(sys.argv[curr_arg],'r')

# Get resource layer
curr_arg += 1
resfile = open(sys.argv[curr_arg],'r')

# Get box left and right coords
curr_arg += 1
bl = int(sys.argv[curr_arg])
curr_arg += 1
br = int(sys.argv[curr_arg])

# Get box top and bottom coords
curr_arg += 1
bt = int(sys.argv[curr_arg])
curr_arg += 1
bb = int(sys.argv[curr_arg])

# Get lattice spacing
curr_arg += 1
lat_space = float(sys.argv[curr_arg])

# Get units
curr_arg += 1
unit = float(sys.argv[curr_arg])

# Get start value for calculating KDEs
curr_arg += 1
start_val = int(sys.argv[curr_arg])

# Take points every "subsampling" timesteps 
curr_arg += 1
subsampling = int(sys.argv[curr_arg])

# Get file for saving plot
curr_arg += 1
savefile = sys.argv[curr_arg]

# Get the resource layer
z_array = [[]]
resfile.seek(0)
y_val = -1
for line in resfile:
    y_val += 1
    if y_val != 0:
        z_array += [[]]
    split_line = line.rsplit()
    for x_val in range(len(split_line)):
        z_array[y_val] += [float(split_line[x_val])]
resfile.close()
bottomcontour = -4
topcontour = 4
contourres = 1

######################################################
# First plot
######################################################

# Read in the positions from file
x_array1 = []
y_array1 = []
x_array2 = []
y_array2 = []
x_array3 = []
y_array3 = []
x_array4 = []
y_array4 = []
row = 0
pos_file1.seek(0)
for line in pos_file1:
    row += 1
    split_line = line.rsplit()
    if (row > start_val) and (row % subsampling == 0):
        x_array1 += [float(split_line[0])+random.random()-0.5]
        y_array1 += [float(split_line[1])+random.random()-0.5]
        x_array2 += [float(split_line[2])+random.random()-0.5]
        y_array2 += [float(split_line[3])+random.random()-0.5]
        x_array3 += [float(split_line[4])+random.random()-0.5]
        y_array3 += [float(split_line[5])+random.random()-0.5]
        x_array4 += [float(split_line[6])+random.random()-0.5]
        y_array4 += [float(split_line[7])+random.random()-0.5]
pos_file1.close()

# Calculate the smoothing parameters (h)
h_val1 = math.sqrt((numpy.var(x_array1)+numpy.var(y_array1))/2)*(float(len(x_array1))**(-1.0/6.0))
h_val2 = math.sqrt((numpy.var(x_array2)+numpy.var(y_array2))/2)*(float(len(x_array2))**(-1.0/6.0))
h_val3 = math.sqrt((numpy.var(x_array3)+numpy.var(y_array3))/2)*(float(len(x_array3))**(-1.0/6.0))
h_val4 = math.sqrt((numpy.var(x_array4)+numpy.var(y_array4))/2)*(float(len(x_array4))**(-1.0/6.0))

# Calculate the KDEs
pd1 =[[]]
pd2 =[[]]
pd3 =[[]]
pd4 =[[]]
for y_val in range(0,bb-bt):
    if y_val != 0:
        pd1 += [[]]
        pd2 += [[]]
        pd3 += [[]]
        pd4 += [[]]
    for x_val in range(0,br-bl):
        pd1[y_val] += [0]
        pd2[y_val] += [0]
        pd3[y_val] += [0]
        pd4[y_val] += [0]
        for data in range(len(x_array1)):
            pd1[y_val][x_val] += (1/(unit**2))*math.exp(-(((x_val+bl)*lat_space-x_array1[data])**2+((y_val+bt)*lat_space-y_array1[data])**2)/(2*(h_val1**2)))/(2*math.pi*len(x_array1)*(h_val1**2))
            pd2[y_val][x_val] += (1/(unit**2))*math.exp(-(((x_val+bl)*lat_space-x_array2[data])**2+((y_val+bt)*lat_space-y_array2[data])**2)/(2*(h_val2**2)))/(2*math.pi*len(x_array2)*(h_val2**2))
            pd3[y_val][x_val] += (1/(unit**2))*math.exp(-(((x_val+bl)*lat_space-x_array3[data])**2+((y_val+bt)*lat_space-y_array3[data])**2)/(2*(h_val3**2)))/(2*math.pi*len(x_array3)*(h_val3**2))
            pd4[y_val][x_val] += (1/(unit**2))*math.exp(-(((x_val+bl)*lat_space-x_array4[data])**2+((y_val+bt)*lat_space-y_array4[data])**2)/(2*(h_val4**2)))/(2*math.pi*len(x_array4)*(h_val4**2))

# Plot KDEs
fig = pylab.figure()
fig.set_size_inches(12,12)
fig.add_subplot(2,2,1)
pylab.hold('on')
filled_contours = pylab.contourf(z_array, origin='lower', extent=[0,(br-bl)*lat_space*unit,0,(bb-bt)*lat_space*unit], levels=[float(x)/contourres for x in range(bottomcontour,topcontour)],cmap=pylab.cm.Greens)
pylab.ylabel('Northing',fontsize=20)
pd_levels = [0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.1]
pylab.contour(pd1, origin='lower',colors='k', extent=[0,(br-bl)*lat_space*unit,0,(bb-bt)*lat_space*unit],levels=pd_levels,linewidths=2)
pylab.contour(pd2, origin='lower',colors='b', extent=[0,(br-bl)*lat_space*unit,0,(bb-bt)*lat_space*unit],levels=pd_levels,linewidths=2)
pylab.contour(pd3, origin='lower',colors='y', extent=[0,(br-bl)*lat_space*unit,0,(bb-bt)*lat_space*unit],levels=pd_levels,linewidths=2)
pylab.contour(pd4, origin='lower',colors='m', extent=[0,(br-bl)*lat_space*unit,0,(bb-bt)*lat_space*unit],levels=pd_levels,linewidths=2)
pylab.text(2,93,'a)',fontsize=26)

######################################################
# Second plot
######################################################

# Read in the positions from file
x_array1 = []
y_array1 = []
x_array2 = []
y_array2 = []
x_array3 = []
y_array3 = []
x_array4 = []
y_array4 = []
row = 0
pos_file2.seek(0)
for line in pos_file2:
    row += 1
    split_line = line.rsplit()
    if (row > start_val) and (row % subsampling == 0):
        x_array1 += [float(split_line[0])+random.random()-0.5]
        y_array1 += [float(split_line[1])+random.random()-0.5]
        x_array2 += [float(split_line[2])+random.random()-0.5]
        y_array2 += [float(split_line[3])+random.random()-0.5]
        x_array3 += [float(split_line[4])+random.random()-0.5]
        y_array3 += [float(split_line[5])+random.random()-0.5]
        x_array4 += [float(split_line[6])+random.random()-0.5]
        y_array4 += [float(split_line[7])+random.random()-0.5]
pos_file2.close()

# Calculate the smoothing parameters (h)
h_val1 = math.sqrt((numpy.var(x_array1)+numpy.var(y_array1))/2)*(float(len(x_array1))**(-1.0/6.0))
h_val2 = math.sqrt((numpy.var(x_array2)+numpy.var(y_array2))/2)*(float(len(x_array2))**(-1.0/6.0))
h_val3 = math.sqrt((numpy.var(x_array3)+numpy.var(y_array3))/2)*(float(len(x_array3))**(-1.0/6.0))
h_val4 = math.sqrt((numpy.var(x_array4)+numpy.var(y_array4))/2)*(float(len(x_array4))**(-1.0/6.0))

# Calculate the KDEs
pd1 =[[]]
pd2 =[[]]
pd3 =[[]]
pd4 =[[]]
for y_val in range(0,bb-bt):
    if y_val != 0:
        pd1 += [[]]
        pd2 += [[]]
        pd3 += [[]]
        pd4 += [[]]
    for x_val in range(0,br-bl):
        pd1[y_val] += [0]
        pd2[y_val] += [0]
        pd3[y_val] += [0]
        pd4[y_val] += [0]
        for data in range(len(x_array1)):
            pd1[y_val][x_val] += (1/(unit**2))*math.exp(-(((x_val+bl)*lat_space-x_array1[data])**2+((y_val+bt)*lat_space-y_array1[data])**2)/(2*(h_val1**2)))/(2*math.pi*len(x_array1)*(h_val1**2))
            pd2[y_val][x_val] += (1/(unit**2))*math.exp(-(((x_val+bl)*lat_space-x_array2[data])**2+((y_val+bt)*lat_space-y_array2[data])**2)/(2*(h_val2**2)))/(2*math.pi*len(x_array2)*(h_val2**2))
            pd3[y_val][x_val] += (1/(unit**2))*math.exp(-(((x_val+bl)*lat_space-x_array3[data])**2+((y_val+bt)*lat_space-y_array3[data])**2)/(2*(h_val3**2)))/(2*math.pi*len(x_array3)*(h_val3**2))
            pd4[y_val][x_val] += (1/(unit**2))*math.exp(-(((x_val+bl)*lat_space-x_array4[data])**2+((y_val+bt)*lat_space-y_array4[data])**2)/(2*(h_val4**2)))/(2*math.pi*len(x_array4)*(h_val4**2))

# Plot KDEs
fig.add_subplot(2,2,2)
pylab.hold('on')
pd_levels = [0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.1]
filled_contours = pylab.contourf(z_array, origin='lower', extent=[0,(br-bl)*lat_space*unit,0,(bb-bt)*lat_space*unit], levels=[float(x)/contourres for x in range(bottomcontour,topcontour)],cmap=pylab.cm.Greens)
pylab.contour(pd1, origin='lower',colors='k', extent=[0,(br-bl)*lat_space*unit,0,(bb-bt)*lat_space*unit],levels=pd_levels,linewidths=2)
pylab.contour(pd2, origin='lower',colors='b', extent=[0,(br-bl)*lat_space*unit,0,(bb-bt)*lat_space*unit],levels=pd_levels,linewidths=2)
pylab.contour(pd3, origin='lower',colors='y', extent=[0,(br-bl)*lat_space*unit,0,(bb-bt)*lat_space*unit],levels=pd_levels,linewidths=2)
pylab.contour(pd4, origin='lower',colors='m', extent=[0,(br-bl)*lat_space*unit,0,(bb-bt)*lat_space*unit],levels=pd_levels,linewidths=2)
pylab.text(2,93,'b)',fontsize=26)

######################################################
# Third plot
######################################################

# Read in the positions from file
x_array1 = []
y_array1 = []
x_array2 = []
y_array2 = []
x_array3 = []
y_array3 = []
x_array4 = []
y_array4 = []
row = 0
pos_file3.seek(0)
for line in pos_file3:
    row += 1
    split_line = line.rsplit()
    if (row > start_val) and (row % subsampling == 0):
        x_array1 += [float(split_line[0])+random.random()-0.5]
        y_array1 += [float(split_line[1])+random.random()-0.5]
        x_array2 += [float(split_line[2])+random.random()-0.5]
        y_array2 += [float(split_line[3])+random.random()-0.5]
        x_array3 += [float(split_line[4])+random.random()-0.5]
        y_array3 += [float(split_line[5])+random.random()-0.5]
        x_array4 += [float(split_line[6])+random.random()-0.5]
        y_array4 += [float(split_line[7])+random.random()-0.5]
pos_file3.close()

# Calculate the smoothing parameters (h)
h_val1 = math.sqrt((numpy.var(x_array1)+numpy.var(y_array1))/2)*(float(len(x_array1))**(-1.0/6.0))
h_val2 = math.sqrt((numpy.var(x_array2)+numpy.var(y_array2))/2)*(float(len(x_array2))**(-1.0/6.0))
h_val3 = math.sqrt((numpy.var(x_array3)+numpy.var(y_array3))/2)*(float(len(x_array3))**(-1.0/6.0))
h_val4 = math.sqrt((numpy.var(x_array4)+numpy.var(y_array4))/2)*(float(len(x_array4))**(-1.0/6.0))

# Calculate the KDEs
pd1 =[[]]
pd2 =[[]]
pd3 =[[]]
pd4 =[[]]
for y_val in range(0,bb-bt):
    if y_val != 0:
        pd1 += [[]]
        pd2 += [[]]
        pd3 += [[]]
        pd4 += [[]]
    for x_val in range(0,br-bl):
        pd1[y_val] += [0]
        pd2[y_val] += [0]
        pd3[y_val] += [0]
        pd4[y_val] += [0]
        for data in range(len(x_array1)):
            pd1[y_val][x_val] += (1/(unit**2))*math.exp(-(((x_val+bl)*lat_space-x_array1[data])**2+((y_val+bt)*lat_space-y_array1[data])**2)/(2*(h_val1**2)))/(2*math.pi*len(x_array1)*(h_val1**2))
            pd2[y_val][x_val] += (1/(unit**2))*math.exp(-(((x_val+bl)*lat_space-x_array2[data])**2+((y_val+bt)*lat_space-y_array2[data])**2)/(2*(h_val2**2)))/(2*math.pi*len(x_array2)*(h_val2**2))
            pd3[y_val][x_val] += (1/(unit**2))*math.exp(-(((x_val+bl)*lat_space-x_array3[data])**2+((y_val+bt)*lat_space-y_array3[data])**2)/(2*(h_val3**2)))/(2*math.pi*len(x_array3)*(h_val3**2))
            pd4[y_val][x_val] += (1/(unit**2))*math.exp(-(((x_val+bl)*lat_space-x_array4[data])**2+((y_val+bt)*lat_space-y_array4[data])**2)/(2*(h_val4**2)))/(2*math.pi*len(x_array4)*(h_val4**2))

# Plot KDEs
fig.add_subplot(2,2,3)
pylab.hold('on')
filled_contours = pylab.contourf(z_array, origin='lower', extent=[0,(br-bl)*lat_space*unit,0,(bb-bt)*lat_space*unit], levels=[float(x)/contourres for x in range(bottomcontour,topcontour)],cmap=pylab.cm.Greens)
pylab.ylabel('Northing',fontsize=20)
pylab.xlabel('Easting',fontsize=20)
pd_levels = [0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.1]
pylab.contour(pd1, origin='lower',colors='k', extent=[0,(br-bl)*lat_space*unit,0,(bb-bt)*lat_space*unit],levels=pd_levels,linewidths=2)
pylab.contour(pd2, origin='lower',colors='b', extent=[0,(br-bl)*lat_space*unit,0,(bb-bt)*lat_space*unit],levels=pd_levels,linewidths=2)
pylab.contour(pd3, origin='lower',colors='y', extent=[0,(br-bl)*lat_space*unit,0,(bb-bt)*lat_space*unit],levels=pd_levels,linewidths=2)
pylab.contour(pd4, origin='lower',colors='m', extent=[0,(br-bl)*lat_space*unit,0,(bb-bt)*lat_space*unit],levels=pd_levels,linewidths=2)
pylab.text(2,93,'c)',fontsize=26)

######################################################
# Fourth plot
######################################################

# Read in the positions from file
x_array1 = []
y_array1 = []
x_array2 = []
y_array2 = []
x_array3 = []
y_array3 = []
x_array4 = []
y_array4 = []
row = 0
pos_file4.seek(0)
for line in pos_file4:
    row += 1
    split_line = line.rsplit()
    if (row > start_val) and (row % subsampling == 0):
        x_array1 += [float(split_line[0])+random.random()-0.5]
        y_array1 += [float(split_line[1])+random.random()-0.5]
        x_array2 += [float(split_line[2])+random.random()-0.5]
        y_array2 += [float(split_line[3])+random.random()-0.5]
        x_array3 += [float(split_line[4])+random.random()-0.5]
        y_array3 += [float(split_line[5])+random.random()-0.5]
        x_array4 += [float(split_line[6])+random.random()-0.5]
        y_array4 += [float(split_line[7])+random.random()-0.5]
pos_file4.close()

# Calculate the smoothing parameters (h)
h_val1 = math.sqrt((numpy.var(x_array1)+numpy.var(y_array1))/2)*(float(len(x_array1))**(-1.0/6.0))
h_val2 = math.sqrt((numpy.var(x_array2)+numpy.var(y_array2))/2)*(float(len(x_array2))**(-1.0/6.0))
h_val3 = math.sqrt((numpy.var(x_array3)+numpy.var(y_array3))/2)*(float(len(x_array3))**(-1.0/6.0))
h_val4 = math.sqrt((numpy.var(x_array4)+numpy.var(y_array4))/2)*(float(len(x_array4))**(-1.0/6.0))

# Calculate the KDEs
pd1 =[[]]
pd2 =[[]]
pd3 =[[]]
pd4 =[[]]
for y_val in range(0,bb-bt):
    if y_val != 0:
        pd1 += [[]]
        pd2 += [[]]
        pd3 += [[]]
        pd4 += [[]]
    for x_val in range(0,br-bl):
        pd1[y_val] += [0]
        pd2[y_val] += [0]
        pd3[y_val] += [0]
        pd4[y_val] += [0]
        for data in range(len(x_array1)):
            pd1[y_val][x_val] += (1/(unit**2))*math.exp(-(((x_val+bl)*lat_space-x_array1[data])**2+((y_val+bt)*lat_space-y_array1[data])**2)/(2*(h_val1**2)))/(2*math.pi*len(x_array1)*(h_val1**2))
            pd2[y_val][x_val] += (1/(unit**2))*math.exp(-(((x_val+bl)*lat_space-x_array2[data])**2+((y_val+bt)*lat_space-y_array2[data])**2)/(2*(h_val2**2)))/(2*math.pi*len(x_array2)*(h_val2**2))
            pd3[y_val][x_val] += (1/(unit**2))*math.exp(-(((x_val+bl)*lat_space-x_array3[data])**2+((y_val+bt)*lat_space-y_array3[data])**2)/(2*(h_val3**2)))/(2*math.pi*len(x_array3)*(h_val3**2))
            pd4[y_val][x_val] += (1/(unit**2))*math.exp(-(((x_val+bl)*lat_space-x_array4[data])**2+((y_val+bt)*lat_space-y_array4[data])**2)/(2*(h_val4**2)))/(2*math.pi*len(x_array4)*(h_val4**2))

# Plot KDEs
fig.add_subplot(2,2,4)
pylab.hold('on')
filled_contours = pylab.contourf(z_array, origin='lower', extent=[0,(br-bl)*lat_space*unit,0,(bb-bt)*lat_space*unit], levels=[float(x)/contourres for x in range(bottomcontour,topcontour)],cmap=pylab.cm.Greens)
pylab.xlabel('Easting',fontsize=20)
pd_levels = [0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.1]
pylab.contour(pd1, origin='lower',colors='k', extent=[0,(br-bl)*lat_space*unit,0,(bb-bt)*lat_space*unit],levels=pd_levels,linewidths=2)
pylab.contour(pd2, origin='lower',colors='b', extent=[0,(br-bl)*lat_space*unit,0,(bb-bt)*lat_space*unit],levels=pd_levels,linewidths=2)
pylab.contour(pd3, origin='lower',colors='y', extent=[0,(br-bl)*lat_space*unit,0,(bb-bt)*lat_space*unit],levels=pd_levels,linewidths=2)
pylab.contour(pd4, origin='lower',colors='m', extent=[0,(br-bl)*lat_space*unit,0,(bb-bt)*lat_space*unit],levels=pd_levels,linewidths=2)
pylab.text(2,93,'d)',fontsize=26)

# Save and show figure
pylab.savefig(savefile)
pylab.show()