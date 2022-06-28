#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# MPA size decision support model for SA case study by Nils Krueck, Jan 2022          @
# Script to run version 1.00 (Nils Krueck), 20 Feb 2022                               @
#   - version run_SA_MPA_sizer_v100_nils - adapted from the first generic run script  @
#   - folders specific to Nils' computer                                              @
#                                                                                     @
# Version: nils02, adding code to include runs based on mean distance                                                                             @
#          nils03, adding code to include runs based on mean max distance
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

rm(list = ls()) # delete parameters in the workspace
library(stringr)
library(dplyr)
library(doParallel)

# Specify model and parameterisation
modelid <- 100 # specify model version / id
resolution <- 100 # modelling resolution in m
rBRUVcatchment <- 100 # radius of BRUV catchment
dispersal.period <- 'Weekly' # specify dispersal period
age <- 'Max' # specify age threshold as either "Max" or "Maturity"
case.study.region <- 'SA' # specify region name in case location names are not specified or pooled
location.names <- case.study.region # set to 'NA' unless data is supposed to be pooled!!!!  
MPAsizes <- c(0,seq(500,2000,500),seq(3000,5000,1000),seq(10000,30000,10000),50000,100000) # specify MPA sizes to be analysed (in m)
fmort <- seq(0,1,0.05)  # fishing mortality rate when fully exposed (discrete proportion per year) 
nreplicates <- 100 # number of replicate simulations to run
mean.extent <- 1 # specify extent of modeling environment based on either mean (1) or max movement distances (0)  
idealized.maxn <- 1 # use uniform maxn data (>0) to get idealized results from movement profiles
mean.max.dist <- 1 # use the grand mean across individual movement distances to 

if(idealized.maxn == 1){case.study.region <- paste0(case.study.region,'_idealized')}
if(mean.max.dist == 1){case.study.region <- paste0(case.study.region,'_meanMaxDist')}

# Load model and data
modelversion <- paste0("SA_MPAsizer",modelid)
modelname <- paste0('SA_MPASizer_v',modelid,'.R')
scenario.name <- paste0(case.study.region,'_Res', resolution,'_',dispersal.period,'_AgeAt',age,'_rBRUV',rBRUVcatchment)
maindir <- paste0(getwd(),'/')
data.folder <- 'Data/'
code.folder <- 'R Code/IBMs/'
if(!dir.exists(paste0(getwd(),'/Results/'))){
  dir.create(paste0(getwd(),'/Results/'))}
results.folder <- paste0('Results/', modelversion,'/') #,'_',scenario.name,'/') # specify results folder
datadir <- paste0(maindir,data.folder)
codedir <- paste0(maindir,code.folder)
source(paste0(codedir,modelname)) # run MPAsizer function

num.data.all <- readRDS(paste0(datadir,'MaxNs with absences.RDS')) # load BURV data
mov.data.all <- readRDS(paste0(datadir,'Movement_Profiles/IndividualMovementProbabilities_Measured_',resolution,'m.RDS')) # load dispersal profiles
species.names <- names(mov.data.all) # get species names
nspecies <- length(species.names)   # specify number of species
age.data <- read.csv(paste0(datadir,'input_parameters_updated.csv')) # load age input data
scientific.names <- age.data$genus.species[match(species.names,age.data$common.name)] # get scientific names
mean.max.distances <- read.csv(paste0(datadir,'mean_max_distance_m_',tolower(dispersal.period),'.csv'))

### Optimize performance in R
nCores<-detectCores() #On Klaas' & Nils' laptops this returns 8, but returns above 4 are diminishing as they are not true separate cores
cl <- makePSOCKcluster(nCores)
registerDoParallel(cl)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#### Run the model for single species at a time @
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

returndat <- list() # aloocate list to store model output
#spec.data <- matrix(nrow=nspecies,ncol=8)

print(scenario.name)

for (s in 1:nspecies) {
  
  mov.data <- mov.data.all[[s]]

  # choose dispersal duration and specify data
  if(dispersal.period=='Daily'){mov.data <- mov.data$Daily}
  if(dispersal.period=='Weekly'){mov.data <- mov.data$Weekly}
  if(dispersal.period=='Monthly'){mov.data <- mov.data$Monthly}
  
  # confirm that resolution matches up
  if(dim(mov.data[[1]])[1] > 1 & resolution != abs(mov.data[[1]]$Distance_Metres[1]-mov.data[[1]]$Distance_Metres[2])){
    stop("Resolution does not match movement probability data.")
  }
  
  species.name <- species.names[s]
  
  print(paste0('Species ',s,'/',nspecies,' : ',species.name)) # "Western_Blue_Groper": Achoerodus gouldii
  
  # match BRUV data to species name 
  species.num.data.loc <- which(names(num.data.all)==species.name) # identify data for species
  num.data <- num.data.all[[species.num.data.loc]] # need to confirm species names!!!!!
  
  # match age data to species name 
  species.age.data.loc <- which(age.data$common.name==species.name) # identify data for species
  if(age=='Max'){max.age <- age.data$max.age[species.age.data.loc]}
  if(age=='Maturity'){max.age <- age.data$age.at.maturity[species.age.data.loc]}
  
  # specify locations to be analysed (Pre-specified names means data are pooled across locations)
  if(is.na(location.names)){location.names <- unique(num.data$location)
  } else {num.data$location <- location.names}
  
  # add idealized number at each location
  if(idealized.maxn > 0){num.data$maxn <- idealized.maxn}
  
  # override data with mean distance profile if species
  if(mean.max.dist == 1){
    species.mean.dist.data.loc <- which(mean.max.distances$common.name==species.name)
    mdata <- data.frame(Distance_Metres = seq(-ceiling(mean.max.distances$mean.max.dist[species.mean.dist.data.loc]/100)*100,
                                    ceiling(mean.max.distances$mean.max.dist[species.mean.dist.data.loc]/100)*100,100))
    mdata$Probability_Occurrence <- 1/dim(mdata)[1]
    mov.data <- list()
    mov.data[[1]] <- mdata
    names(mov.data) <- "Mean_species_dist"
  }
  
  print(paste0(species.name, ', mean.max.dist = ', mean.max.distances$mean.max.dist[species.mean.dist.data.loc]))
  
  #rBRUVcatchment <- rBRUVcatchments[s] # could introduce an internal loop 
  
  # FROM OLD SCRIPT TO SAVE METADATA
  #print.data <- data.frame(Species = species.name, Length = max.length, Density_100m2 = mean.density*100,
  #                        Mean_HR = mean.distance,SD_HR = sd.distance, HR_Method = method, Max_age = max.age, Age_method = age.method)
  #spec.data[s,] <- as.matrix(print.data)
  #print(print.data)
  if(!is.na(max.age)&max(num.data$maxn,na.rm=TRUE)>0){

    output <- SA_MPASizer_v100(species.name,   # name of species to simulate
                                 location.names, # names of locations represented by different density values (optional)
                                 MPAsizes,       # vector of MPA sizes to be simulated for conservation impact testing
                                 num.data,       # number of individuals / abundance data - define format 
                                 mov.data,       # movement data consisting of a list that contains individual movement probability kernels inspecified modelling resolution  
                                 rBRUVcatchment, # Radius in m of the plume catchment area around BRUVs (area of attraction)
                                 fmort,          # fishing mortality rate(s) when fully exposed (discrete proportion per year)
                                 max.age,        # maximum age of species
                                 nreplicates, # number of replicate simulations to run
                                 resolution = 100, # resolution of the modelling environment in m
                                 mean.extent = 1) # specify extent of modeling environment based on either mean (1) or max movement distances (0)
  
    output$scenario_name <- scenario.name
    returndat[[s]] <- output
  } else {
    returndat[[s]] <- NA
  }
}

names(returndat) <- species.names[1:nspecies]
#returndat$scenario_name <- scenario.name

# create and assign folders to save results
if(!dir.exists(paste0(getwd(),'/',results.folder))){ #species.name, "_rBRUVcatchment_",rBRUVcatchment,"/", speciesToAnalyse[s]))){
  dir.create(paste0(getwd(),'/',results.folder)) #"/rBRUVcatchment_",rBRUVcatchment,"/", speciesToAnalyse[s]))}
}

#spec.data <- setNames(as.data.frame(spec.data),names(print.data))
#save.image(paste0(getwd(),'/',results.folder, modelversion,'_', scenario.name,'.RData'))
saveRDS(returndat,paste0(getwd(),'/',results.folder, modelversion,'_', scenario.name,'_returndat.RDS'))

