# Script to generate metadata and variable resolutions of individual movement profiles from tracking data
# Version 1: Nils Krueck, 02/11/21
library(stringr)

rm(list = ls()) # remove all parameters from workspace
data.folder <- paste0(getwd(),'/') #Data/')#'C:/Users/nkrueck/Documents/Github/SA_MPA_model/Data'
#setwd(data.folder)
resolution <- 100 # in m
type <- 'Measured' # 'Measured' or 'Idealized'

# get data on dispersal
probsi <- readRDS(paste0('IndividualMovementProbabilities_',type,'.RDS'))
probsm <- readRDS(paste0('MeanMovementProbabilities_',type,'.RDS'))

# save metadata in loop if desirable
meta.data <- list()
inames <- list() 
inums <- matrix()
meta.data[[1]] <- names(probsi)

for (s in 1:length(probsi)){
  
  #probsi_s <- probsi[[s]]
  #probsm_s <- probsm[[s]]
  meta.data[[2]] <- names(probsi[[s]]) # names of periods
 
  for (p in 1:length(probsi[[s]])){
  
    inames[[s]] <- names(probsi[[s]][[p]]) # names of individuals 
    inums[s] <- length(probsi[[s]][[p]]) # names of individuals 
    
    # calculate revised mean probabilities across individuals
    if(length(probsm[[s]][[p]]$Distance_Metres)>2*resolution){
    
      newmprobs <- matrix()
      
      mdistances <- probsm[[s]][[p]]$Distance_Metres
      mprobs <- probsm[[s]][[p]]$Probability_Occurrence
      
      #recalculate distances
      bdist <- ceiling(max(mdistances)/resolution)*resolution
      newmdistances <- c(-rev(c(seq(resolution,bdist,resolution))),
                         seq(0,bdist,resolution))
      
      newmdistances <- newmdistances[abs(newmdistances)<=max(mdistances)]
      
      # determine startlocs and endlocs
      startlocs <- match(newmdistances - resolution/2,mdistances)
      endlocs <- match(newmdistances + resolution/2,mdistances)
      
      # couple of corrections needed
      if(is.na(startlocs[1])){
        startlocs[1] <- 1}
      #if(startlocs[1] > 1){
      #  startlocs <- c(1,startlocs)
      #  endlocs <- c(startlocs[2],endlocs)
      #  newmdistances <- c(min(mdistances),newmdistances)}
      
      if(is.na(endlocs[length(endlocs)])){
        endlocs[length(endlocs)] <- length(mdistances)}
      #if(endlocs[length(endlocs)] < length(mdistances)){
      #  startlocs <- c(startlocs,endlocs[length(endlocs)]+1)
      #  endlocs <- c(endlocs,length(mdistances))
      #  newmdistances <- c(newmdistances,max(mdistances))}
      
      # round distances to nearest upper value matching the resolution  
      newmdistances <- ceiling(abs(newmdistances)/resolution)*resolution*(newmdistances/abs(newmdistances))
      newmdistances[is.na(newmdistances)] <- 0
      
      # calculate new probs
      for(rows in 1:length(startlocs)){
        newmprobs[rows] <- sum(mprobs[startlocs[rows]:endlocs[rows]],na.rm = TRUE)
      }
      newmprobs <- newmprobs/sum(newmprobsna.rm=TRUE)
    
    } else {
      newmprobs <- 1
      newmdistances <- 0    
    }    
      
    # add new data to list
    probsm[[s]][[p]] <- data.frame(Probability_Occurrence = newmprobs,
                                        Distance_Metres = newmdistances)
    
    # calculate revised probabilities for individuals
    probsi_p <- probsi[[s]][[p]]
    
    for (i in 1:length(probsi_p)){
      
      if(max(probsi_p[[i]]$Distance_Metres)>2*resolution){
      
        newiprobs <- matrix()
        
        idistances <- probsi_p[[i]]$Distance_Metres
        iprobs <- probsi_p[[i]]$Probability_Occurrence
        
        #recalculate distances
        bdist <- ceiling(max(idistances)/resolution)*resolution
        newidistances <- c(-rev(c(seq(resolution,bdist,resolution))),
                           seq(0,bdist,resolution))
        
        newidistances <- newidistances[abs(newidistances)<=max(idistances)]
        
        
        # determine startlocs and endlocs
        startlocs <- match(newidistances - resolution/2,idistances)
        endlocs <- match(newidistances + resolution/2,idistances)
        
        # couple of corrections needed
        if(is.na(startlocs[1])){
          startlocs[1] <- 1}
        
        if(is.na(endlocs[length(endlocs)])){
          endlocs[length(endlocs)] <- length(mdistances)}
        
    
        # round distances to nearest upper value matching the resolution  
        newidistances <- ceiling(abs(newidistances)/resolution)*resolution*(newidistances/abs(newidistances))
        newidistances[is.na(newidistances)] <- 0
        
        # calculate new probs
        for(rows in 1:length(startlocs)){
          newiprobs[rows] <- sum(iprobs[startlocs[rows]:endlocs[rows]],na.rm = TRUE)
        }
        newiprobs <- newiprobs/sum(newiprobs,na.rm = TRUE)
        
      } else {# length condition
        newiprobs <- 1
        newidistances <- 0
      }
      # add new data to list
      probsi[[s]][[p]][[i]] <- data.frame(Probability_Occurrence = newiprobs,
                                            Distance_Metres = newidistances)
    } # individual loop
  } # period loop
} # species loop


# get data on dispersal
saveRDS(probsi,paste0('IndividualMovementProbabilities_',type,'_',resolution,'m.RDS'))
saveRDS(probsm,paste0('MeanMovementProbabilities_',type,'_',resolution,'m.RDS'))

# check distances
#print(probsm[[s]][[p]]$Distance_Metres)
#print(probsi[[s]][[p]][[25]]$Distance_Metres)


meta.data[[3]] <- inums
meta.data[[4]] <- inames 
names(meta.data) <- c('Species_Names','Time_Periods','Number_of_Individuals_Per_Species','IDs_Of_Individuals_Per_Species')
names(meta.data[[4]]) <- meta.data[[1]]
saveRDS(meta.data,'metadata_MovementProfiles.RDS')
