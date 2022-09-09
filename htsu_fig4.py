################################################################################
# Name: htsu_fig4.py
#
# Purpose: Plot Figure 4
# 
# Usage: python htsu_fig4.py htsu_plot_bif.inp htsu_plot_bif.png
################################################################################

import sys
from matplotlib import pyplot as plt

# Get input file from command line
curr_arg = 1
infile = open(sys.argv[curr_arg],'r')
curr_arg += 1
savefile = sys.argv[curr_arg]

# Average for RW
rwave = 16.50394

# Analytic point of bifurcation
bif_an = 1.875

# Get information from file. 
beta_fwd = []
beta_back = []
seg = []
agg_back = []
agg_fwd = []
for curr_line in infile:
  split_line=curr_line.rsplit()
  beta_fwd += [float(split_line[0])]
  seg += [float(split_line[1])]
  agg_fwd += [float(split_line[2])]
  if split_line[3] != "NA":
    beta_back += [float(split_line[0])]
    agg_back += [float(split_line[3])]
  
# Plot the data
plt.subplot(1,2,1)
plt.hold('on')
fig = plt.gcf()
fig.set_size_inches(12,5)
plt.axis([1.45,3.55,15,35])
plt.xlabel("Strength of avoidance, $a$", fontsize=16)
plt.ylabel(r"Extent of segregation", fontsize=16)
plt.text(1.5,33,'a)', fontsize=24)
plt.scatter(beta_fwd, seg, c='k',edgecolor='face',s=10,label="Simulation output")
plt.plot([1.5,2,2.5,3,3.5],[rwave,rwave,rwave,rwave,rwave],'k--',linewidth=1,label="No segregation")
plt.plot([bif_an,bif_an],[15,20],'b--',linewidth=1,label="Turing bifurcation")
plt.legend(loc=(0.02,0.63), fontsize=14)
plt.subplot(1,2,2)
plt.axis([1.45,3.55,10,110])
plt.text(1.5,100,'b)', fontsize=24)
plt.xlabel("Strength of attraction, $a$", fontsize=16)
plt.ylabel(r"Extent of aggregation", fontsize=16)
plt.scatter(beta_fwd, agg_fwd, c='k',edgecolor='face',s=10,label="Increasing attraction")
plt.scatter(beta_back, agg_back, c='r',edgecolor='face',s=10,label="Decreasing attraction")
plt.plot([1.5,2,2.5,3,3.5],[rwave,rwave,rwave,rwave,rwave],'k--',linewidth=1,label="No aggregation")
beta_min = 1.875
beta_max = 2.025
plt.plot([beta_min,beta_min],[15,100],'b--',linewidth=1,label=r"$a_{\rm min}$")
plt.plot([beta_max,beta_max],[15,100],'b--',linewidth=1,label=r"$a_{\rm max}$")
plt.legend(loc=7, fontsize=14)

# Save and show figure
plt.savefig(savefile)
plt.show()                  