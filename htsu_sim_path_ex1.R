###############################################################################
# Name: htsu_sim_path_ex1.R
#
# Purpose: Simulate path through a resource layer
#           --> see Suplementary Appendix A: "Example movement kernels"
#
# Author: Jonathan R. Potts
#
# Feel free to share and adapt, but giving appropriate credit
###############################################################################

# Clear the memory
rm(list=ls());gc()

# Raster library for plot
library(raster)

# Get landscape layer as a matrix
layer<-as.matrix(read.table('random_field_100.inp',sep='\t',header=FALSE))

beta<-1.5       # beta_1 from Equation (S.2) in Supplementary Appendix A
lambda<-0.2     # lambda from Equation (S.2) in Supplementary Appendix A
step_no<-1000       # number of steps to be simulated

box_width<-100        # Width of layer
box_height<-100       # Height of layer
loc_x<-c(box_width/2+1)  
loc_y<-c(box_height/2+1)               # start in centre of box
y_vals<-rep(1:box_width,box_height)    # y-values at each point in layer

# Calculate x-values at each point in layer
x_vals<-c()
for(count in 1:box_height)
{
  x_vals<-c(x_vals,rep(count,box_width))  
}

# Simulate path
for(step in 2:step_no)
{
  # Movement kernel prior to normalising
  unnorm_mk<-exp(-lambda*sqrt((x_vals-loc_x[step-1])^2+(y_vals-loc_y[step-1])^2))*exp(beta*layer)
  # Normalise the movement kernel
  mk<-unnorm_mk/sum(unnorm_mk)
  # Draw sample from the movement kernel to find new location
  newxy<-sample(box_width*box_height,1,prob=mk)
  loc_x<-c(loc_x,((newxy-1) %% box_width)+1)
  loc_y<-c(loc_y,((newxy-1) %/% box_width)+1)
}

loc_x    # Write values of x-locations to screen
loc_y    # Write values of y-locations to screen

# Plot the layer as a raster
xy<-cbind(x_vals,y_vals)
ras<-raster(ncols=box_width, nrows=box_height,xmn=0,xmx=box_width,ymn=0,ymx=box_height)
rlayer<-rasterize(xy,ras,layer)
plot(rlayer)

# Plot the animal locations
points(loc_x,loc_y)