# Function to generate idealized (kernel-based) connectivity matrices according to 

# 1) distance matrix (distmat) specifying distance between planning units in m, and  
# 2) estimated larval dispersal distance (ldist) in m (mean +- SD), and/or 
# 3) estimated adult home range movement distance (adist) in m (mean +- SD)

# Output:
# A list termed conmats, containing
# 1) a 2D matrix of larval dispersal probabilities (proportions of all) normalized to 1 for each source location
# 2) a 2D matrix of adult residence probabilities (proportions of all) normalized to 1 for each source location (centre of home range) 
# Note: Both matrices represent sources in rows and destinations in columns
# Note: Probabilities are calculated based on the negative binomial probability density function suitable for over-dispersed count data (no negative values)

# Version 01 by Nils Krueck, November 2020

generate_conmats_func <- function(distmat,mean_ldist,sd_ldist,mean_adist,sd_adist) {

  rm(list = ls()) # remove all parameters from workspace
  
  # calculate larval dispersal probability matrix
  
  if(!is.na(mean_ldist)){ # if larval dispersal distance is specified
    
    if(is.na(sd_ldist)){ # if standard deviation of larval dispersal distance is unspecified
      sd_ldist <- mean_ldist # set equal to mean distance
      print('Standard deviation of larval dispersal distance not specified - setting equal to mean (CV = 1).')
    }
    
    cv_ldist = sd_ldist/mean_ldist # coefficient of variation
    p_ldist = 1/cv_ldist/sd_ldist # probability of sampling success to represent SD
    r_ldist = mean_ldist * p_ldist / (1-p_ldist); # change of mean to sample in order to maintain actual mean 
      
    probs_ldisp <- 1 - pnbinom(distmat,size=r_ldist,prob = p_ldist) # calculate probabilities of larval dispersal
    rowsums <- rowSums(probs_ldisp,na.rm = TRUE) # calculate sums across columns
    nprobs_ldisp <- apply(probs_ldisp,1,function(x) x/rowsums) # normalize dispersal probabilities
    
  } else {
    nprobs_ldisp <- NA # else return NA
  }
  
  # calculate adult movement probability matrix

  if(!is.na(mean_adist)){ # if adult movement distance is specified
  
    mean_adist <- round(mean_adist/2) # divide home range by half, assuming that fish have a centre within their home range and are equally likely to move in any direction within 2D space  
    
    if(is.na(sd_adist)){ # if standard deviation of adult movement distance is unspecified
      sd_adist <- mean_adist # set equal to mean distance
      print('Standard deviation of adult movement distance is not specified - setting equal to mean (CV = 1).')
    } else {(sd_adist <- round(sd_adist/2))} # see above - diving by half to account for movements from origin rather than across entire range
  
  cv_adist = sd_adist/mean_adist # coefficient of variation
  p_adist = 1/cv_adist/sd_adist # probability of sampling success to represent SD
  r_adist = mean_adist * p_adist / (1-p_adist); # change of mean to sample in order to maintain actual mean 
  
  probs_adisp <- 1 - pnbinom(distmat,size=r_adist,prob = p_adist) # calculate probabilities of larval dispersal
  rowsums <- rowSums(probs_adisp,na.rm = TRUE) # calculate sums across columns
  nprobs_adisp <- apply(probs_adisp,1,function(x) x/rowsums) # normalize dispersal probabilities
  
  } else {
  nprobs_adisp <- NA # else return NA
  }

# OUTPUT --

conmats <- list(nprobs_ldisp,nprobs_adisp)
names(conmats) <- c("Probability_Larvae","Probabilities_Adults")
return(conmats)

} # end of function