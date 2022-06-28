# Script to analyse individual movement profiles generated from tracking data
# Version 1: Nils Krueck, 16/03/22
library(stringr)

rm(list = ls()) # remove all parameters from workspace
#data.folder <- paste(getwd(),'/Data/Movement_Profiles/')# 'C:/Users/nkrueck/Documents/Github/ParksSA/Data/'
data.folder <- 'C:/Users/nkrueck/Documents/Github/ParksSA/Data/Movement_Profiles/'
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

#save.folder <- 'C:/Users/nkrueck/Documents/Github/SA_MPA_model/Figures/'
species.names <- names(probsm)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Plot movement profiles ####
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@

species.to.plot <- "All" #"Bluethroat_Wrasse" #c("Western_Blue_Groper","Snapper")
if (species.to.plot == "All"){species.to.plot <- species.names}

srp_sum <- list()
sp_sum <- list()

mrp_sum <- matrix()
sdrp_sum <- matrix()
minrp_sum <- matrix()
maxrp_sum <- matrix()
medrp_sum <- matrix()
max_dist <- matrix()

mp_sum <- matrix()
sdp_sum <- matrix()
minp_sum <- matrix()
maxp_sum <- matrix()
medp_sum <- matrix()
cvp_sum <- matrix()

for (s in 1:length(species.to.plot)){
  
  probsi_s <- probsi[[which(!is.na(match(species.names,species.to.plot[s])))]]
  probsm_s <- probsm[[which(!is.na(match(species.names,species.to.plot[s])))]]
  rp_max <- matrix()
  rp_sum <- matrix()
  p_sum <- matrix()
  
  for (p in 2){  #1:length(probsm_s)){
    
    period <- names(probsm_s)[p]
    #dir.create(paste0(save.folder,period,'/'))
    
    # call plot
    #png(paste0(save.folder,period,'/',species.to.plot[s],'MovementProbabilities_',type,'_',period,'_',resolution,'m.png'), 
    #    units="cm", width=14, height=12, res=300)
    
    for (i in 1:length(probsi_s[[p]])){
      
      mdistances <- round(probsm_s[[p]]$Distance_Metres/resolution)*resolution
      idistances <- round(probsi_s[[p]][[i]]$Distance_Metres/resolution)*resolution #probsi_s[[p]][[i]]$Distance_Metres
      iprobs <- probsi_s[[p]][[i]]$Probability_Occurrence
      iprobs.locs <- which(!is.na(match(mdistances,idistances)))
      plot.probs <- rep(0,length(probsm_s[[p]]$Distance_Metres))
      plot.probs[iprobs.locs] <- iprobs[]/max(iprobs,na.rm=TRUE)
 
      nprobs <- length(iprobs)
      rp_max[i] <- iprobs[1]/(1/nprobs) # relative probablity at maximum movement distances
      p_sum[i] <- sum(plot.probs) # relative probablity at maximum movement distances
      rp_sum[i] <- sum(plot.probs)/nprobs # relative probablity at maximum movement distances
      #uneven.locs <- which(round(idistances)/resolution!=round(idistances/resolution))
      
      # if(i == 1) {plot(mdistances,plot.probs,type='l',col='grey50',
      #                  xlab = "Movement Distance (m)",
      #                  ylab = "Relative Probability",
      #                  main = paste0(str_replace_all(species.to.plot[s],"_"," ")," (",names(probsi_s)[p],")"))
      # } else {lines(probsm_s[[p]]$Distance_Metres,plot.probs,col='grey50')}
      # #lines(idistances,iprobs,col='grey50')}
      # 
      
    } # individual loop
    
    srp_sum[[s]] <- rp_sum;
    mrp_sum[s] <- mean(rp_sum)
    sdrp_sum[s] <- sd(rp_sum)
    minrp_sum[s] <- min(rp_sum)
    maxrp_sum[s] <- max(rp_sum)
    medrp_sum[s] <- median(rp_sum)
 
    sp_sum[[s]] <- p_sum;
    mp_sum[s] <- mean(p_sum)
    sdp_sum[s] <- sd(p_sum)
    minp_sum[s] <- min(p_sum)
    maxp_sum[s] <- max(p_sum)
    medp_sum[s] <- median(p_sum)
    cvp_sum[s] <- mp_sum[s]/sdp_sum[s]
    
    max_dist[s] <- round(abs(mdistances[1])*2/1000)
    
    
    
    #lines(mdistances,probsm_s[[p]]$Probability_Occurrence/max(probsm_s[[p]]$Probability_Occurrence,na.rm=TRUE),
    #      col='black',lwd=2)
    
    # end plot
    #dev.off()
    
  } # period loop
} # species loop


combined.data <- rbind(max_dist,round(rbind(mrp_sum,sdrp_sum,medrp_sum,minrp_sum,maxrp_sum)*100))
colnames(combined.data) <- species.names
rownames(combined.data) <- c('Max_distance','Mean_Psum','SD_Psum','Median_Psum','Min_Psum','Max_Psum')

names(srp_sum) <- species.names
names(mrp_sum) <- species.names
names(sdrp_sum) <- species.names
names(minrp_sum) <- species.names
names(maxrp_sum) <- species.names
names(medrp_sum) <- species.names

