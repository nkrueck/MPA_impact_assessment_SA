# Script to plot individual movement profiles generated from tracking data
# Version 1: Nils Krueck, 31/10/21
library(stringr)

rm(list = ls()) # remove all parameters from workspace
data.folder <- paste0(getwd(),'/Data/Movement_Profiles/') #'C:/Users/nkrueck/Documents/Github/ParksSA/Data/'
#setwd(data.folder)
resolution <- 100 # in m
type <- 'Measured' # 'Measured' or 'Idealized'

# get data on movement distances
if (resolution > 1){
  probsi <- readRDS(paste0(data.folder,'IndividualMovementProbabilities_',type,'_',resolution,'m.RDS'))
  probsm <- readRDS(paste0(data.folder,'MeanMovementProbabilities_',type,'_',resolution,'m.RDS'))
} else {
  probsi <- readRDS(paste0(data.folder,'IndividualMovementProbabilities_',type,'.RDS'))
  probsm <- readRDS(paste0(data.folder,'MeanMovementProbabilities_',type,'.RDS'))
}
  
# get meta data
readRDS(paste0(getwd(),'/Data/Movement_Profiles/metadata_MovementProfiles.RDS'))

save.folder <- paste0(getwd(),'/Dispersal_kernels/') #'C:/Users/nkrueck/Documents/Github/SA_MPA_model/Figures/'
species.names <- names(probsm)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Plot movement profiles ####
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@

species.to.plot <- "All" #"Bluethroat_Wrasse" #c("Western_Blue_Groper","Snapper")
if (species.to.plot == "All"){species.to.plot <- species.names}

for (s in 1:length(species.to.plot)){

  probsi_s <- probsi[[which(!is.na(match(species.names,species.to.plot[s])))]]
  probsm_s <- probsm[[which(!is.na(match(species.names,species.to.plot[s])))]]

  for (p in 1:length(probsm_s)){

    period <- names(probsm_s)[p]
    dir.create(paste0(save.folder,period,'/'))
      
    # call plot
    png(paste0(save.folder,period,'/',species.to.plot[s],'MovementProbabilities_',type,'_',period,'_',resolution,'m.png'), 
        units="cm", width=14, height=12, res=300)
   
    for (i in 1:length(probsi_s[[p]])){
      
      mdistances <- round(probsm_s[[p]]$Distance_Metres/resolution)*resolution
      idistances <- round(probsi_s[[p]][[i]]$Distance_Metres/resolution)*resolution #probsi_s[[p]][[i]]$Distance_Metres
      iprobs <- probsi_s[[p]][[i]]$Probability_Occurrence
      iprobs.locs <- which(!is.na(match(mdistances,idistances)))
      #uneven.locs <- which(round(idistances)/resolution!=round(idistances/resolution))
      plot.probs <- rep(0,length(probsm_s[[p]]$Distance_Metres))
      plot.probs[iprobs.locs] <- iprobs[]/max(iprobs,na.rm=TRUE)
      
      if(i == 1) {plot(mdistances,plot.probs,type='l',col='grey50',
                       xlab = "Movement Distance (m)",
                       ylab = "Relative Probability",
                       main = paste0(str_replace_all(species.to.plot[s],"_"," ")," (",names(probsi_s)[p],")"))
      } else {lines(probsm_s[[p]]$Distance_Metres,plot.probs,col='grey50')}
        #lines(idistances,iprobs,col='grey50')}
      
    } # individual loop
    lines(mdistances,probsm_s[[p]]$Probability_Occurrence/max(probsm_s[[p]]$Probability_Occurrence,na.rm=TRUE),
          col='black',lwd=2)
    
    # end plot
    dev.off()
    
  } # period loop
} # species loop


# species.name <- species.names[4]
# ind.num <- 4
# mprobs <- mprobsi_species[which(names(mprobsi_species)==species.name)][[1]][[1]][[ind.num]]$Probability_Occurrence
# mdistances <- mprobsi_species[which(names(mprobsi_species)==species.name)][[1]][[1]][[ind.num]]$Distance_Metres
# 
# plot(mdistances,mprobs,type='l',main=species.name,
#      xlab = "Movement Distance (m)",
#      ylab = "Probability")
# 
# eprobsi <- eprobsi_species[which(names(eprobsi_species)==species.name)][[1]][[1]][[ind.num]]$Probability_Occurrence
# edistancesi <- eprobsi_species[which(names(eprobsi_species)==species.name)][[1]][[1]][[ind.num]]$Distance_Metres
# 
# plot(edistancesi,eprobsi,type='l',main=species.name,
#      xlab = "Movement Distance (m)",
#      ylab = "Probability")
# 
# eprobs <- eprobs_species[which(names(eprobs_species)==species.name)][[1]][[1]]$Probability_Occurrence
# edistances <- eprobs_species[which(names(eprobs_species)==species.name)][[1]][[1]]$Distance_Metres
# 
# plot(edistances,eprobs,type='l',main=species.name,
#      xlab = "Movement Distance (m)",
#      ylab = "Probability")

