###############################################################################
# Name: htsu_bm_ud_ex1.R
#
# Purpose: An example of calculating the UD from a movement kernel, using the
#          Barnett-Moorcroft formula.  The animal's movement is biased towards 
#          higher quality resources, given by a single resource layer.  Also 
#          calculates the approximations given in Equations (12) and (13).
#           --> see Suplementary Appendix F: "Calculating the steady state UD
#          using the Barnett-Moorcroft method"
#
# Author: Jonathan R. Potts
#
# Feel free to share and adapt, but giving appropriate credit
###############################################################################

# Clear the memory
rm(list=ls());gc()

# Get landscape layer as a matrix
layer<-as.matrix(read.table('random_field_100.inp',sep='\t',header=FALSE))

beta<-1.5       # beta_R from the paper and Supplementary Appendices
lambda<-0.2     # lambda from the paper and Supplementary Appendices

box_width<-100        # Width of layer
box_height<-100       # Height of layer

y_vals<-rep(1:box_width,box_height)      # y-values at each point in layer

# Calculate x-values at each point in layer
x_vals<-c()
for(count in 1:box_height)
{
  x_vals<-c(x_vals,rep(count,box_width))  
}

# Barnett-Moorcroft UD from Equation (11).  This requires an integration.  The integral becomes 
# a sum when translated into computer code.
ud<-c()
for(x_val in 1:box_height)
{
  for(y_val in 1:box_width)
  {
    # Value of the UD at this point, prior to normalising
    ud<-c(ud,exp(beta*layer[(x_val-1)*box_height+y_val])*sum(exp(-lambda*sqrt((x_vals-x_val)^2+(y_vals-y_val)^2))*exp(beta*layer)))    # Add this value to the UD array
  }
}
ud<-ud/sum(ud)  # normalise the UD

# Estimations from Equation (12) (the limit as phi is arbitrarily narrow) and (13) (the limit as 
# phi samples from the whole landscape)
ud1<-exp(2*beta*layer)/sum(exp(2*beta*layer))
ud2<-exp(beta*layer)/sum(exp(beta*layer))

# Plot the UD as a raster
library(raster)
xy<-cbind(x_vals,y_vals)
ras<-raster(ncols=box_width, nrows=box_height,xmn=0,xmx=box_width,ymn=0,ymx=box_height)
rud<-rasterize(xy,ras,ud)
plot(rud)