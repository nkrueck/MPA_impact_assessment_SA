# Script to run loop over species data for SA case study
# Nils Krueck, Jan 2022, version 1.00

# NEED INFORMATION ON:
# - rBRUVcatchment by species.name
# - fmorts by species.name
# - max age by species.name 

rm(list = ls()) # delete parameters in the workspace

modelversion <- "sa_mpasizer100"
modelname <- 'SA_MPASizer_v100.R'
resolution <- 100

scenario.name <- paste0('SA_Res', resolution)
maindir <- 'C:/Users/nkrueck/Dropbox/Projects/MPA SA project with Flinders Uni/SA_MPA_model/'
data.folder <- 'Github Dec 2021/Data/'
if(!dir.exists(paste0(getwd(),'/Results/'))){
  dir.create(paste0(getwd(),'/Results/'))}
results.folder <- paste0('Results/', modelversion,'_',scenario.name,'/') # specify results folder
datadir <- paste0(maindir,data.folder)

library(stringr)
library(dplyr)
library(doParallel)

# load data
num.data.all <- readRDS(paste0(datadir,'MaxNs with absences.RDS'))
mov.data.all <- readRDS(paste0(datadir,'IndividualMovementProbabilities_Measured_',resolution,'m.RDS'))
species.names <- names(mov.data.all)
nspecies <- length(species.names)  

# Prepare input data
#MPAsizes <- c(seq(0,1000,1000)) #,seq(2000,5000,1000),10000)      # vector of MPA sizes to be simulated for conservation impact testing 
MPAsizes <- c(0,seq(500,2000,500),seq(3000,5000,1000),seq(10000,30000,10000),50000,100000)
rBRUVcatchments <- 200

max.ages <- 15
fmort <- seq(0,0.5,.1)  # fishing mortality rate when fully exposed (discrete proportion per year) 
nreplicates <- 1000 # number of replicate simulations to run
mean.extent <- 1 # specify extent of modeling environment based on either mean (1) or max movement distances (0)  

### Optimize performance

nCores<-detectCores() #On Klaas' & Nils' laptops this returns 8, but returns above 4 are diminishing as they are not true separate cores
cl <- makePSOCKcluster(nCores)
registerDoParallel(cl)

#### Run the model for single species at a time


# run MPAsizer function
source(paste0(modelname))

returndat <- list()
#spec.data <- matrix(nrow=nspecies,ncol=8)
location.names <- "SA"
nspecies <- 1 # dummy

for (s in 1:nspecies) {
  
  mov.data <- mov.data.all[[s]]
  mov.data <- mov.data$Weekly # hardwired to weekly analysis
  
  # confirm that resolution matches up
  if(resolution != abs(mov.data[[1]]$Distance_Metres[1]-mov.data[[1]]$Distance_Metres[2])){
    stop("Resolution does not match movement probability data.")
  }
  
  print(paste0('Species ',s,'/',nspecies)) # "Western_Blue_Groper": Achoerodus gouldii
  
  # match BRUV data to movement data 
  species.name <- species.names[s]
  species.num.data.loc <- 5 # which(names(num.dat.all)==species.name) # identify data for species
  num.data <- num.data.all[[species.num.data.loc]] # need to confirm species names!!!!!
  
  max.age <- max.ages[s]
  rBRUVcatchment <- rBRUVcatchments[s]
  
  # FROM OLD SCRIPT TO SAVE METADATA
  #print.data <- data.frame(Species = species.name, Length = max.length, Density_100m2 = mean.density*100,
  #                        Mean_HR = mean.distance,SD_HR = sd.distance, HR_Method = method, Max_age = max.age, Age_method = age.method)
  #spec.data[s,] <- as.matrix(print.data)
  #print(print.data)
  

  returndat[[s]] <- SA_MPASizer_v100(species.name,   # name of species to simulate
                               location.names, # names of locations represented by different density values (optional)
                               MPAsizes,       # vector of MPA sizes to be simulated for conservation impact testing
                               num.data,       # number of individuals / abundance data - define format 
                               mov.data,       # movement data consisting of a list that contains individual movement probability kernels inspecified modelling resolution  
                               rBRUVcatchment, # Radius in m of the plume catchment area around BRUVs (area of attraction)
                               fmort,          # fishing mortality rate(s) when fully exposed (discrete proportion per year)
                               max.age,        # maximum age of species
                               nreplicates = 1000, # number of replicate simulations to run
                               resolution = 100, # resolution of the modelling environment in m
                               mean.extent = 1) # specify extent of modeling environment based on either mean (1) or max movement distances (0)

}

names(returndat) <- species.names[1:nspecies]

# create and assign folders to save results
if(!dir.exists(paste0(getwd(),'/',results.folder))){ #species.name, "_rBRUVcatchment_",rBRUVcatchment,"/", speciesToAnalyse[s]))){
  dir.create(paste0(getwd(),'/',results.folder)) #"/rBRUVcatchment_",rBRUVcatchment,"/", speciesToAnalyse[s]))}
}

#spec.data <- setNames(as.data.frame(spec.data),names(print.data))
save.image(paste0(getwd(),'/',results.folder, modelversion,'_', scenario.name,'.RData'))
saveRDS(returndat,paste0(getwd(),'/',results.folder, modelversion,'_', scenario.name,'_returndat.RDS'))

