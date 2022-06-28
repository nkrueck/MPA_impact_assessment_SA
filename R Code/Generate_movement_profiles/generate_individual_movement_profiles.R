# Script to generate individual movement profiles from tracking data
# Version 1: Nils Krueck, 31/10/21
library(stringr)

rm(list = ls()) # remove all parameters from workspace
function.folder <- paste0(getwd(),'/R Code/Generate_movement_profiles') #'C:/Users/nkrueck/Dropbox/Projects/MPA SA project with Flinders Uni/SA_MPA_model/new_scripts'
#setwd(function.folder)
source(paste0(function.folder,'/generate_conmats_func.R')) # run function to get idealized movement profiles

# get data on dispersal
data.folder <- paste0(getwd(),'/Data') #'C:/Users/nkrueck/Documents/Github/SA_MPA_model/Data'
#setwd(data.folder)
dispersal.files <- list.files(data.folder,pattern='Dispersal_Timescales_*',all.files=FALSE)
dispersal.files <- sort(dispersal.files) #[2:length(dispersal.files)])
nspecies <- length(dispersal.files)
species.names <- matrix()

#@@@@@@@@@@@@@@@@@@@@@@@@
# Main species loop #####
#@@@@@@@@@@@@@@@@@@@@@@@@

# allocate lists
mprobsm_species <- list() # mean across individual movement probability profiles as measured in the field
mprobsi_species <- list() # individual movement probability profiles as measured in the field
eprobs_species <- list() # movement probability profiles estimated from means and std across individuals
eprobsm_species <- list() # mean across individual movement probability profiles estimated from means and std
eprobsi_species <- list() # individual movement probability profiles estimated from means and std

for (s in 1:nspecies){  
  
  # load data
  sdata <- readRDS(paste0(data.folder,'/',dispersal.files[s]))
  nperiods <- length(sdata)
  species.name <- str_replace_all(str_to_title(unique(sdata[[1]]$species_common_name))," ","_")
  species.names[s] <- species.name[species.name != "Na"] # exclude NAs
  
  # allocate lists
  mprobsi_period <- list()
  mprobsm_period <- list()
  eprobs_period <- list()
  eprobsi_period <- list()
  eprobsm_period <- list()
  
  for (p in 1:nperiods){
    
    # allocate lists
    mprobsi <- list() # empty list for movement probabilities
    eprobsi <- list() # empty list for movement probabilities
    
    # get rid of zeros
    zero_ids <- aggregate(list(maxd=sdata[[p]]$maxDistkm),list(id=sdata[[p]]$transmitter_id),max,na.rm=TRUE)
    zero_ids <- zero_ids$id[which(zero_ids$maxd==0)]
    if(length(zero_ids)>0){
      keep.locs <- which(is.na(match(sdata[[p]]$transmitter_id,zero_ids)))
      sdata[[p]] <- sdata[[p]][keep.locs,]}
    #which(is.na(sdata[[p]]$maxDistkm))
    
    ids <- unique(sdata[[p]]$transmitter_id) # get individual ids
    mean_adist_all <- mean(sdata[[p]]$maxDistkm*1000,na.rm=TRUE)
    sd_adist_all <- sd(sdata[[p]]$maxDistkm*1000,na.rm=TRUE)
    cv_adist_all <- sd_adist_all/mean_adist_all
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    # estimate movement probabilities from means and stds
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
    if(mean_adist_all != 0){
      if(sd_adist_all < 0.25 * mean_adist_all | is.na(sd_adist_all)){sd_adist_all <- cv_adist_all*mean_adist_all}
      max_dist_all <- round(mean_adist_all)*2 # extent of modelling environment
      distmat <- rbind(seq(max_dist_all,0,-1),seq(0,max_dist_all,1)) # generate two-dimensional matrix 
      mean_adist <- mean_adist_all
      sd_adist <- sd_adist_all
      mean_ldist <- NA
      sd_ldist <- NA
      probmats <- generate_conmats_func(distmat,mean_ldist,sd_ldist,mean_adist,sd_adist)
      eprobs_distances <- c(seq(-max_dist_all,-1,1),seq(0,max_dist_all,1)) # generate distances profiles for plotting - deleting one of the zero distance values
      eprobs_probs <- c(probmats$Probabilities_Adults[1:max_dist_all,1],probmats$Probabilities_Adults[,2]) # bind matrices together - deleting one of the zero distance estimates
      # find first and last zero probabilities
      min0loc <- min(which(eprobs_probs!=0))
      max0loc <- max(which(eprobs_probs!=0))
      eprobs_probs <- eprobs_probs[min0loc:max0loc]/sum(eprobs_probs[min0loc:max0loc],na.rm = TRUE)
      eprobs_distances <- eprobs_distances[min0loc:max0loc]
    } else {
        eprobs_probs <- 1
        eprobs_distances <- 0
    }
  
    eprobs_period[[p]] <- data.frame(Probability_Occurrence = eprobs_probs,
                              Distance_Metres = eprobs_distances)
  
    rm(probmats,eprobs_distances,eprobs_probs)
    
    # loop over individuals
    for (id in 1:length(ids)) {
      
      alldists <- sdata[[p]]$maxDistkm[sdata[[p]]$transmitter_id == ids[id]] * 1000
      
      #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      # estimate individual movement probabilities from means and stds
      #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      mean_adist <- mean(alldists,na.rm=TRUE) # mean movement distance 
      sd_adist <- sd(alldists,na.rm=TRUE) # sd
      
      if(mean_adist != 0){
        if(sd_adist < 0.25 * mean_adist | is.na(sd_adist)){sd_adist <- cv_adist_all*mean_adist}
        max_dist <- round(mean_adist)*2 # extent of modelling environment
        distmat <- rbind(seq(max_dist,0,-1),seq(0,max_dist,1)) # generate two-dimensional matrix 
        probmats <- generate_conmats_func(distmat,mean_ldist,sd_ldist,mean_adist,sd_adist)
        eprobs_distances <- c(seq(-max_dist,-1,1),seq(0,max_dist,1)) # generate distances profiles for plotting - deleting one of the zero distance values
        eprobs_probs <- c(probmats$Probabilities_Adults[1:max_dist,1],probmats$Probabilities_Adults[,2]) # bind matrices together - deleting one of the zero distance estimates
        # find first and last zero probabilities
        min0loc <- min(which(eprobs_probs!=0))
        max0loc <- max(which(eprobs_probs!=0))
        eprobs_probs <- eprobs_probs[min0loc:max0loc]/sum(eprobs_probs[min0loc:max0loc],na.rm = TRUE)
        eprobs_distances <- eprobs_distances[min0loc:max0loc]
      } else {
        eprobs_probs <- 1
        eprobs_distances <- 0
      }
      
      eprobsi[[id]] <- data.frame(Probability_Occurrence = eprobs_probs,
                                 Distance_Metres = eprobs_distances)

      #print(c(min(eprobsi[[id]]$Distance_Metres),max(eprobsi[[id]]$Distance_Metres)))
      
      # add individual movement probabilities together
      if(id==1){eprobsm <- eprobsi[[id]] # assign first list 
      } else { # add data to list
        if(length(eprobs_distances) > length(eprobsm$Distance_Metres)){ # id travelled further than previous individuals
          edistances <- eprobs_distances
          eprobs.locs <- which(!is.na(match(edistances,eprobsm$Distance_Metres)))
          eprobs <- eprobs_probs
          eprobs[eprobs.locs] <- eprobs[eprobs.locs] + eprobsm$Probability_Occurrence
        } else { # if id did not travel as widely as previous ones
          eprobs.locs <- which(!is.na(match(eprobsm$Distance_Metres,eprobs_distances)))
          edistances <- eprobsm$Distance_Metres
          eprobs <- eprobsm$Probability_Occurrence
          eprobs[eprobs.locs] <- eprobs[eprobs.locs] + eprobs_probs
        }
        eprobsm <- data.frame(Probability_Occurrence = eprobs/sum(eprobs.na.rm=TRUE),
                              Distance_Metres = edistances)
      }
      
      #print(c(min(eprobsm$Distance_Metres),max(eprobsm$Distance_Metres)))
      
      #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      # estimate individual movement probabilities directly from measurements
      #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      
      maxd.kernel <- max(alldists,na.rm=TRUE) 
      kernel <- as.matrix(rep(1,1,maxd.kernel+(maxd.kernel%%2 == 0)*1)) # generate kernel + 1 for even numbers
      centreloc <- ceiling(max(alldists)/2)
    
      # loop over individual travel distances
      for (idist in 1:length(unique(alldists))) { # for each observed travel distance
          addval <- 1 # set initial number of observed vals to 1
          maxlocs <- which(alldists == max(alldists)) # identify all distance locs
          startloc <- floor(centreloc - alldists[maxlocs[1]]/2) # define kernel start loc
          endloc <- ceiling(startloc + alldists[maxlocs[1]]) # define kernel end loc
          if (length(maxlocs) > 1){addval <- length(maxlocs)} # if observed more than once update value
          kernel[startloc:endloc] <- kernel[startloc:endloc] + addval # add observations along kernel spectrum
          alldists <- alldists[alldists != alldists[maxlocs[1]]] # delete observed distance from list
      }
  
      # get cumulative probability spectrum by dividing by sum
      prob <- kernel/sum(kernel,na.rm=TRUE)
      distances <- (1:length(kernel))-(ceiling(length(kernel)/2))
      
      if(sum(sdata[[p]]$maxDistkm[sdata[[p]]$transmitter_id == ids[id]],na.rm=TRUE) == 0){ # make sure that probs and distances match up there is just zeros
         prob <- 1
         distances <- 0
      }
      
      mprobsi[[id]] <- data.frame(Probability_Occurrence = prob,
                                 Distance_Metres = distances)
  
      #print(c(min(mprobsi[[id]]$Distance_Metres),max(mprobsi[[id]]$Distance_Metres)))
      
      # add individual movement probabilities together
      if(id==1){mprobsm <- mprobsi[[id]]} # assign first list 
      if(id > 1){ # add data to list
        if(length(distances) > length(mprobsm$Distance_Metres)){ # id travelled further than previous individuals
          mdistances <- distances
          mprobs.locs <- which(!is.na(match(distances,mprobsm$Distance_Metres)))
          mprobs <- prob
          mprobs[mprobs.locs] <- mprobs[mprobs.locs] + mprobsm$Probability_Occurrence
        } else { # if id did not travel as widely as previous ones
          mprobs.locs <- which(!is.na(match(mprobsm$Distance_Metres,distances)))
          mdistances <- mprobsm$Distance_Metres
          mprobs <- mprobsm$Probability_Occurrence
          mprobs[mprobs.locs] <- mprobs[mprobs.locs] + prob
        }
        mprobsm <- data.frame(Probability_Occurrence = mprobs/sum(mprobs,na.rm = TRUE),
                              Distance_Metres = mdistances)
      }
      
      #print(c(min(mprobsm$Distance_Metres),max(mprobsm$Distance_Metres)))
      
        
    } # end of individual loop
    
    # save probability profiles per period
    names(mprobsi) <- ids # assign names to individual list objects 
    mprobsi_period[[p]] <- mprobsi

    names(eprobsi) <- ids # assign names to individual list objects 
    eprobsi_period[[p]] <- eprobsi
    
    # normalize summed probs to represent averages
    #eprobsm$Probability_Occurrence <- eprobsm$Probability_Occurrence/sum(eprobsm$Probability_Occurrence,na.rm = TRUE)
    eprobsm_period[[p]] <- eprobsm
    #mprobsm$Probability_Occurrence <- mprobsm$Probability_Occurrence/sum(mprobsm$Probability_Occurrence,na.rm = TRUE)
    mprobsm_period[[p]] <- mprobsm
    
  } # time period loop 
  names(mprobsi_period) <- names(sdata)  
  mprobsi_species[[s]] <- mprobsi_period
  
  names(eprobsi_period) <- names(sdata)  
  eprobsi_species[[s]] <- eprobsi_period

  names(eprobs_period) <- names(sdata)  
  eprobs_species[[s]] <- eprobs_period
  
  names(mprobsm_period) <- names(sdata)  
  mprobsm_species[[s]] <- mprobsm_period
  
  names(eprobsm_period) <- names(sdata)  
  eprobsm_species[[s]] <- eprobsm_period
  
} # species loop

names(mprobsi_species) <- species.names
names(mprobsm_species) <- species.names
names(eprobsi_species) <- species.names
names(eprobs_species) <- species.names
names(eprobsm_species) <- species.names

saveRDS(mprobsi_species,'IndividualMovementProbabilities_Measured.RDS')
saveRDS(mprobsm_species,'MeanMovementProbabilities_Measured.RDS')
saveRDS(eprobsi_species,'IndividualMovementProbabilities_Idealized.RDS')
saveRDS(eprobsm_species,'MeanMovementProbabilities_Idealized.RDS')
saveRDS(eprobs_species,'CombinedMovementProbabilities_Idealized.RDS')
