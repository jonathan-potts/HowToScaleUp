###############################################################################
# Name: htsu_sim_2indivs.R
#
# Purpose: Simulate 2 interacting individuals moving through a resource layer with mutual avoidance
#           --> see Suplementary Appendix G: "Interacting stochastic IBM example"
#
# Author: Jonathan R. Potts
#
# Feel free to share and adapt, but giving appropriate credit
###############################################################################

# Clear the memory
rm(list=ls());gc()

# Get landscape layer as a matrix
layer<-as.matrix(read.table('random_field_100.inp',sep='\t',header=FALSE))

betaR<-0.5              # beta_R from the paper and Supplementary Appendices
betaA<-4000             # Strength of avoidance (beta_{1,2} and beta_{2,1} from Supplementary Appendix G)
step_no = 20000         # Number of steps
burn_in = 19000         # No of steps to burn in
no_steps_for_od = 10000 # Number of steps that OD is calculated over
box_width<-100          # Width of layer
box_height<-100         # Height of layer

dec_inc = 1/no_steps_for_od # Amount to decrease OD by each step

loc_x1<-c(box_width/4)   # Array to store x-locations for individual 1
loc_x2<-c(3*box_width/4) # Array to store x-locations for individual 2
loc_y1<-c(box_height/2)  # Array to store y-locations for individual 1
loc_y2<-c(box_height/2)  # Array to store y-locations for individual 2

occurrence_dist1<-c(rep(0,box_width*box_height))  # Array for occurrence distribution of individual 1
occurrence_dist2<-c(rep(0,box_width*box_height))  # Array for occurrence distribution of individual 2
occurrence_dist1[(loc_x1[1]-1)*box_height+loc_y1]<-1 # Probability of being at start location is 1
occurrence_dist2[(loc_x2[1]-1)*box_height+loc_y2]<-1 # Probability of being at start location is 1

# Simulate paths
for(step in 2:step_no)
{
  ##############################################################################
  # Calculate the probabilities of moving up/down/left/right for individual 1
  ##############################################################################
  if(loc_y1[step-1] == 1)
  {
    # Can't move up if at the top of the domain
    up_prob<-0
  }
  else
  {
    # Include the effect of the resource layer
    up_prob<-exp(betaR*layer[loc_y1[step-1]-1+box_height*(loc_x1[step-1]-1)])
    # Include the effect of the occurrence distribution of the other individual
    up_prob<-up_prob*exp(-betaA*occurrence_dist2[loc_y1[step-1]-1+box_height*(loc_x1[step-1]-1)])
  }
  if(loc_y1[step-1] == box_height)
  {
    # Can't move down if at the bottom of the domain
    down_prob<-0
  }
  else
  {
    # Include the effect of the resource layer
    down_prob<-exp(betaR*layer[loc_y1[step-1]+1+box_height*(loc_x1[step-1]-1)])
    # Include the effect of the occurrence distribution of the other individual
    down_prob<-down_prob*exp(-betaA*occurrence_dist2[loc_y1[step-1]+1+box_height*(loc_x1[step-1]-1)])
  }
  if(loc_x1[step-1] == 1)
  {
    # Can't move left if at the far left of the domain
    left_prob<-0
  }
  else
  {
    # Include the effect of the resource layer
    left_prob<-exp(betaR*layer[loc_y1[step-1]+box_height*(loc_x1[step-1]-2)])
    # Include the effect of the occurrence distribution of the other individual
    left_prob<-left_prob*exp(-betaA*occurrence_dist2[loc_y1[step-1]+box_height*(loc_x1[step-1]-2)])
  }
  if(loc_x1[step-1] == box_width)
  {
    # Can't move right if at the far right of the domain
    right_prob<-0
  }
  else
  {
    # Include the effect of the resource layer
    right_prob<-exp(betaR*layer[loc_y1[step-1]+box_height*(loc_x1[step-1])])
    # Include the effect of the occurrence distribution of the other individual
    right_prob<-right_prob*exp(-betaA*occurrence_dist2[loc_y1[step-1]+box_height*(loc_x1[step-1])])
  }
  tot_prob<-up_prob + down_prob + left_prob + right_prob
  # Normalise probabilities 
  mk<-c(up_prob/tot_prob,down_prob/tot_prob,left_prob/tot_prob,right_prob/tot_prob)
  # Find next location of individual 1
  newxy<-sample(c(-1,1,-2,2),1,prob=mk)
  loc_x1<-c(loc_x1,loc_x1[step-1]+(newxy%/%2)*abs(newxy)%/%2)
  loc_y1<-c(loc_y1,loc_y1[step-1]+newxy*(newxy%%2))
  
  ##############################################################################
  # Calculate the probabilities of moving up/down/left/right for individual 2
  ##############################################################################
  if(loc_y2[step-1] == 1)
  {
    # Can't move up if at the top of the domain
    up_prob<-0
  }
  else
  {
    # Include the effect of the resource layer
    up_prob<-exp(betaR*layer[loc_y2[step-1]-1+box_height*(loc_x2[step-1]-1)])
    # Include the effect of the occurrence distribution of the other individual
    up_prob<-up_prob*exp(-betaA*occurrence_dist1[loc_y2[step-1]-1+box_height*(loc_x2[step-1]-1)])
  }
  if(loc_y2[step-1] == box_height)
  {
    # Can't move down if at the bottom of the domain
    down_prob<-0
  }
  else
  {
    # Include the effect of the resource layer
    down_prob<-exp(betaR*layer[loc_y2[step-1]+1+box_height*(loc_x2[step-1]-1)])
    # Include the effect of the occurrence distribution of the other individual
    down_prob<-down_prob*exp(-betaA*occurrence_dist1[loc_y2[step-1]+1+box_height*(loc_x2[step-1]-1)])
  }
  if(loc_x2[step-1] == 1)
  {
    # Can't move left if at the far left of the domain
    left_prob<-0
  }
  else
  {
    # Include the effect of the resource layer
    left_prob<-exp(betaR*layer[loc_y2[step-1]+box_height*(loc_x2[step-1]-2)])
    # Include the effect of the occurrence distribution of the other individual
    left_prob<-left_prob*exp(-betaA*occurrence_dist1[loc_y2[step-1]+box_height*(loc_x2[step-1]-2)])
  }
  if(loc_x2[step-1] == box_width)
  {
    # Can't move right if at the far right of the domain
    right_prob<-0
  }
  else
  {
    # Include the effect of the resource layer
    right_prob<-exp(betaR*layer[loc_y2[step-1]+box_height*(loc_x2[step-1])])
    # Include the effect of the occurrence distribution of the other individual
    right_prob<-right_prob*exp(-betaA*occurrence_dist1[loc_y2[step-1]+box_height*(loc_x2[step-1])])
  }
  tot_prob<-up_prob + down_prob + left_prob + right_prob
  # Normalise probabilities 
  mk<-c(up_prob/tot_prob,down_prob/tot_prob,left_prob/tot_prob,right_prob/tot_prob)
  # Find next location of individual 2
  newxy<-sample(c(-1,1,-2,2),1,prob=mk)
  loc_x2<-c(loc_x2,loc_x2[step-1]+(newxy%/%2)*abs(newxy)%/%2)
  loc_y2<-c(loc_y2,loc_y2[step-1]+newxy*(newxy%%2))
  
  # Update occurrence distribution
  if(step <= no_steps_for_od)
  {
    # We are within the first no_steps_for_od steps, so decrease the OD at the first step
    occurrence_dist1[loc_y1[1]+box_height*(loc_x1[1]-1)]<-occurrence_dist1[loc_y1[1]+box_height*(loc_x1[1]-1)]-dec_inc
    occurrence_dist2[loc_y2[1]+box_height*(loc_x2[1]-1)]<-occurrence_dist2[loc_y2[1]+box_height*(loc_x2[1]-1)]-dec_inc
  }
  else
  {
    # Deprecate the occurrence distribution at the first position used to calculate this distribution
    occurrence_dist1[loc_y1[step - no_steps_for_od]+box_height*(loc_x1[step - no_steps_for_od]-1)]<-occurrence_dist1[loc_y1[step - no_steps_for_od]+box_height*(loc_x1[step - no_steps_for_od]-1)]-dec_inc
    occurrence_dist2[loc_y2[step - no_steps_for_od]+box_height*(loc_x2[step - no_steps_for_od]-1)]<-occurrence_dist2[loc_y2[step - no_steps_for_od]+box_height*(loc_x2[step - no_steps_for_od]-1)]-dec_inc
  }
  # Increment the occurrence distribution at the current position
  occurrence_dist1[loc_y1[step]+box_height*(loc_x1[step]-1)]<-occurrence_dist1[loc_y1[step]+box_height*(loc_x1[step]-1)]+dec_inc
  occurrence_dist2[loc_y2[step]+box_height*(loc_x2[step]-1)]<-occurrence_dist2[loc_y2[step]+box_height*(loc_x2[step]-1)]+dec_inc
}

# Plot the animal locations
plot(loc_x1,loc_y1)
points(loc_x2,loc_y2,col="red")