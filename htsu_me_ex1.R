###############################################################################
# Name: htsu_me_ex1.R
#
# Purpose: Master equation for calculating UD for an animal moving through
#          a resource layer
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
layer<-as.matrix(read.table('random_field_50.inp',sep='\t',header=FALSE))

beta<-1.5       # beta_R from Equation (S.2) in Supplementary Appendix A
lambda<-0.2     # lambda from Equation (S.2) in Supplementary Appendix A
step_no<-100    # number of steps to be simulated

box_width<-50        # Width of layer
box_height<-50       # Height of layer
ud<-c(rep(0,box_width*box_height)) 
ud[box_height*(box_width/2-1)+box_width/2]=1 # start with animal in centre of box
y_vals<-rep(1:box_width,box_height)          # y-values at each point in layer

# Calculate x-values at each point in layer
x_vals<-c()
for(count in 1:box_height)
{
  x_vals<-c(x_vals,rep(count,box_width))  
}

# Calculate movement kernel
kernel<-c()
for(from_x in 1:box_height)
{
  for(from_y in 1:box_width)
  {
    # Movement kernel prior to normalising
    unnorm_mk<-exp(-lambda*sqrt((x_vals-from_x)^2+(y_vals-from_y)^2))*exp(beta*layer)
    # Normalise the movement kernel
    mk<-unnorm_mk/sum(unnorm_mk)
    # Place into kernel
    kernel<-c(kernel,mk)
  }
}
  
# Solve the master equation over time 
for(step in 2:step_no)
{
  sumval<-c()
  for(to_x in 0:(box_width-1))
  {
    for(to_y in 0:(box_height-1))
    {
      # The sum over P(s|s',t)U(s',t) of the various locations s' in the landscape
      sumval<-c(sumval,sum(ud*kernel[to_x*box_height+to_y+1+
                                       (x_vals-1)*box_height*box_width*box_height+
                                       (y_vals-1)*box_width*box_height]))
    }
  }
      
  # Update utilisation distribution
  ud<-sumval
}

# Plot the UD as a raster
xy<-cbind(x_vals,y_vals)
ras<-raster(ncols=box_width, nrows=box_height,xmn=0,xmx=box_width,ymn=0,ymx=box_height)
rud<-rasterize(xy,ras,ud)
plot(rud)