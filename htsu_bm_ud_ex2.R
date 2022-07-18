###############################################################################
# Name: htsu_bm_ud_ex2.R
#
# Purpose: An example of calculating the UD from a movement kernel, using the
#          Barnett-Moorcroft formula.  The animal's movement is biased towards 
#          higher quality resources, given by a single resource layer, as well as 
#          towards a central point.  Also calculates the approximations given 
#          in Equations (12) and (13).
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

betaR<-1.5       # beta_R from the paper and Supplementary Appendices
betaC<-0.2       # beta_C from the paper and Supplementary Appendices
lambda<-0.2      # lambda from the paper and Supplementary Appendices

box_width<-100        # Width of layer
box_height<-100       # Height of layer

xc<-box_width/2   # x-value of central point
yc<-box_height/2  # y-value of central point

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
    ud<-c(ud,exp(betaR*layer[(x_val-1)*box_height+y_val]-betaC*sqrt((xc-x_val)^2+(yc-y_val)^2))*
             sum(exp(-lambda*sqrt((x_vals-x_val)^2+(y_vals-y_val)^2))*
                   exp(betaR*layer-betaC*sqrt((xc-x_vals)^2+(yc-y_vals)^2))))    # Add this value to the UD array
  }
}
ud<-ud/sum(ud)  # normalise the UD

# Estimations from Equation (12) (the limit as phi is arbitrarily narrow) and (13) (the limit as 
# phi samples from the whole landscape)
ud1<-exp(2*betaR*layer-2*betaC*sqrt((xc-x_vals)^2+(yc-y_vals)^2))/sum(exp(2*betaR*layer-2*betaC*sqrt((xc-x_vals)^2+(yc-y_vals)^2)))
ud2<-exp(betaR*layer-betaC*sqrt((xc-x_vals)^2+(yc-y_vals)^2))/sum(exp(betaR*layer-betaC*sqrt((xc-x_vals)^2+(yc-y_vals)^2)))

# Plot the UD as a raster
library(raster)
xy<-cbind(x_vals,y_vals)
ras<-raster(ncols=box_width, nrows=box_height,xmn=0,xmx=100,ymn=0,ymx=100)
rud<-rasterize(xy,ras,ud)
plot(rud)