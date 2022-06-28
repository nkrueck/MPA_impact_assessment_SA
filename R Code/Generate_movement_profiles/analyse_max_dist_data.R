# Script to generate individual movement profiles from tracking data
# Version 1: Nils Krueck, 31/10/21
library(stringr)

rm(list = ls()) # remove all parameters from workspace
#function.folder <- paste0(getwd(),'/R Code/Generate_movement_profiles') #'C:/Users/nkrueck/Dropbox/Projects/MPA SA project with Flinders Uni/SA_MPA_model/new_scripts'
#setwd(function.folder)
#source(paste0(function.folder,'/generate_conmats_func.R')) # run function to get idealized movement profiles

# get data on dispersal
data.folder <- paste0(getwd(),'/Data') #'C:/Users/nkrueck/Documents/Github/SA_MPA_model/Data'
#setwd(data.folder)
dispersal.files <- list.files(data.folder,pattern='Dispersal_Timescales_*',all.files=FALSE)
dispersal.files <- sort(dispersal.files) #[2:length(dispersal.files)])
nspecies <- length(dispersal.files)
species.names <- matrix()

period.id <- 1
period <- 'weekly'
if(period == 'weekly'){period.id <- 2}
if(period == 'monthly'){period.id <- 3
} else if (period != 'daily' & period.id ==1 ){stop('Pick a valid period of movement (daily, weekly, monthly)')}

# allocate species matrix
mdist.individual <- list()
maxdist.individual <- list()
sddist.individual <- list()

mdist.species <- matrix()
sddist.species <- matrix()
maxdist.species <- matrix()

mdist.species_m <- matrix()
maxdist.species_m <- matrix()
sddist.species_m <- matrix()

#median.dist.individuals<- matrix()

#@@@@@@@@@@@@@@@@@@@@@@@@
# Main species loop #####
#@@@@@@@@@@@@@@@@@@@@@@@@

# allocate lists

for (s in 1:nspecies){  
  
  # load data
  sdata <- readRDS(paste0(data.folder,'/',dispersal.files[s]))
  nperiods <- length(sdata)
  species.name <- str_replace_all(str_to_title(unique(sdata[[1]]$species_common_name))," ","_")
  species.names[s] <- species.name[species.name != "Na"] # exclude NAs
  
  # get rid of zeros
  zero_ids <- aggregate(list(maxd=sdata[[period.id]]$maxDistkm),list(id=sdata[[period.id]]$transmitter_id),max,na.rm=TRUE)
  zero_ids <- zero_ids$id[which(zero_ids$maxd==0)]
  if(length(zero_ids)>0){
    keep.locs <- which(is.na(match(sdata[[period.id]]$transmitter_id,zero_ids)))
    sdata.period <- sdata[[period.id]][keep.locs,]
    } else {sdata.period <- sdata[[period.id]]}
  
  #sdata.period <- sdata[[period.id]]
  mdist.individual[[s]] <- aggregate(sdata.period$maxDistkm,list(id = sdata.period$transmitter_id),mean,na.rm=TRUE)$x
  maxdist.individual[[s]] <- aggregate(sdata.period$maxDistkm,list(id = sdata.period$transmitter_id),max,na.rm=TRUE)$x
  sddist.individual[[s]] <- aggregate(sdata.period$maxDistkm,list(id = sdata.period$transmitter_id),mean,na.rm=TRUE)$x
  #median.dist.individual[[s]] <- aggregate(sdata.period$maxDistkm,list(id = sdata.period$transmitter_id),median,na.rm=TRUE)$x

  #hist(mdist.individual[[s]],main=paste0(species.name,' n = ',length(mdist.individual[[s]])))
  #abline(v=mean(mdist.individual[[s]]),col='red')
  #abline(v=median(mdist.individual[[s]]),col='blue')
  
  mdist.species[s] <- mean(mdist.individual[[s]])
  sddist.species[s] <- sd(mdist.individual[[s]])
  maxdist.species[s] <- mean(maxdist.individual[[s]])
  
  mdist.species_m[s] <- round(mean(mdist.individual[[s]])*1000)
  sddist.species_m[s] <- round(sd(mdist.individual[[s]])*1000)
  maxdist.species_m[s] <- round(mean(maxdist.individual[[s]])*1000)
  
}

mdist.species_m <- as.data.frame(mdist.species_m)
names(mdist.species_m) <- 'mean.dist'
mdist.species_m$common.name <- species.names

sddist.species_m <- as.data.frame(sddist.species_m)
names(sddist.species_m) <- 'sd.dist'
sddist.species_m$common.name <- species.names

mean.maxdist.species_m <- as.data.frame(maxdist.species_m)
names(mean.maxdist.species_m) <- 'mean.max.dist'
mean.maxdist.species_m$common.name <- species.names

write.csv(mdist.species_m,paste0(getwd(),'/Data/mean_distance_m_',period,'.csv'))
write.csv(sddist.species_m,paste0(getwd(),'/Data/sd_distance_m_',period,'.csv'))
write.csv(mean.maxdist.species_m,paste0(getwd(),'/Data/mean_max_distance_m_',period,'.csv'))
