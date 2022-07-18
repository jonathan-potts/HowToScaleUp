#############################################################################
# Name: htsu_bm_from_ssa.R
#
# Purpose: Takes some simulated movement data and a resource layer, uses
#          step selection analysis to parametrise a movement kernel, then
#          uses the Barnett-Moorcroft method to calculate a predicted 
#          utilisation distribution for the simulated animal.
#           --> see Suplementary Appendix F: "Calculating the steady state UD
#          using the Barnett-Moorcroft method"
# 
# Author: Jonathan R. Potts
#
# Feel free to share and adapt
#############################################################################

# Clear the memory
rm(list=ls());gc()

# Get survival package for doing SSA
require('survival')

# Get (simulated) case and control locations
ssfex1<-read.table('sim_ssf_rf1_beta1_1.csv',sep=',',header=TRUE)

# Get the mean step length of observed path
obs<-subset(ssfex1,Observed==1)
sl<-mean(obs$StepLength)

# Get landscape as a matrix
layer<-as.matrix(read.table('random_field_100.inp',sep='\t',header=FALSE))

# Do the SSA
ssf1sim<-clogit(ssfex1$Observed ~ ssfex1$resource + ssfex1$StepLength + strata(ssfex1$strata), data=ssfex1)

# Get the selection-free value of lambda from the observed step length and the correction 
# from SSA (see Forester et al. 2009 for rationale behind this correction).  This is
# the value of lambda that should be used in the movement kernel.
lambda<-(1/sl)-as.numeric(ssf1sim$coefficients[2])

# Grab the beta coefficient for the landscape layer
beta<-as.numeric(ssf1sim$coefficients[1])

box_width = 100
box_height = 100

# Calculate y-values at each point in layer
y_vals<-rep(1:box_width,box_height)     

# Calculate x-values at each point in layer
x_vals<-c()
for(count in 1:box_height)
{
  x_vals<-c(x_vals,rep(count,box_width))  
}

# Calculate the predicted UD using the Barnett-Moorcroft method
ud<-c()
for(x_val in 1:box_height)
{
  for(y_val in 1:box_width)
  {
    # Value of the UD at this point, prior to normalising
    ud<-c(ud,exp(beta*layer[(x_val-1)*box_height+y_val])*sum(exp(-lambda*sqrt((x_vals-x_val)^2+(y_vals-y_val)^2))*exp(beta*layer)))
  }
}
ud<-ud/sum(ud)  # normalise the UD

# Plot the UD as a raster
library(raster)
xy<-cbind(x_vals,y_vals)
ras<-raster(ncols=box_width, nrows=box_height,xmn=0,xmx=box_width,ymn=0,ymx=box_height)
rud<-rasterize(xy,ras,ud)
plot(rud)