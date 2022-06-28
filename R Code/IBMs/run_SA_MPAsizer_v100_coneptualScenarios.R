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
max.age <- 1
case.study.region <- 'conceptual' # specify region name in case location names are not specified or pooled
location.names <- case.study.region # set to 'NA' unless data is supposed to be pooled!!!!  
MPAsizes <- c(seq(500,10000,500),seq(20000,50000,10000)) # specify MPA sizes to be analysed (in m)
fmort <- seq(0,1,0.05)  # fishing mortality rate when fully exposed (discrete proportion per year) 
nreplicates <- 1000 # number of replicate simulations to run
mean.extent <- 1 # specify extent of modeling environment based on either mean (1) or max movement distances (0)  
idealized.maxn <- 1 # use uniform maxn data (>0) to get idealized results from movement profiles
#mean.max.dist <- 1 # use the grand mean across individual movement distances to 
if(idealized.maxn == 1){case.study.region <- paste0(case.study.region,'_idealized')}
disp.val.names <- c('uniform','normal','negbin')

# Load model and data
modelversion <- paste0("SA_MPAsizer",modelid)
modelname <- paste0('SA_MPASizer_v',modelid,'.R')
scenario.name <- paste0(case.study.region,'_rBRUV',rBRUVcatchment)
maindir <- paste0(getwd(),'/')
data.folder <- 'Data/'
code.folder <- 'R Code/IBMs/'
if(!dir.exists(paste0(getwd(),'/Results/'))){
  dir.create(paste0(getwd(),'/Results/'))}
results.folder <- paste0('Results/', modelversion,'/') #,'_',scenario.name,'/') # specify results folder
datadir <- paste0(maindir,data.folder)
codedir <- paste0(maindir,code.folder)
source(paste0(codedir,modelname)) # run MPAsizer function
disp.vals <- readRDS(paste0(data.folder,'hist.metrics.RDS'))

### Optimize performance in R
nCores<-detectCores() #On Klaas' & Nils' laptops this returns 8, but returns above 4 are diminishing as they are not true separate cores
cl <- makePSOCKcluster(nCores)
registerDoParallel(cl)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#### Run the model for single species at a time @
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

returndat <- list() # aloocate list to store model output
num.data <- data.frame(maxn=1,location=location.names)
#spec.data <- matrix(nrow=nspecies,ncol=8)

ndisp.vals <- length(disp.val.names)
print(scenario.name)

for (v in 1:ndisp.vals) {
  
  species.name <- 'test_species'
  save.name <- paste0(disp.val.names[v],'_',scenario.name)
  print(scenario.name)
  disp.data <- disp.vals[v,]
  hist(disp.data)
  
  # generate movement profiles
  mov.data <- list()
  for(dv in 1:length(disp.data)){
    mdist <- round(disp.data[dv]*1000)/2
    mdata <- data.frame(Distance_Metres = seq(-ceiling(mdist/100)*100,
                                    ceiling(mdist/100)*100,100))
    mdata$Probability_Occurrence <- 1/dim(mdata)[1]
    mov.data[[dv]] <- mdata
  }
  names(mov.data) <- paste0('I',1:length(disp.data))

  #rBRUVcatchment <- rBRUVcatchments[s] # could introduce an internal loop 
  
  # FROM OLD SCRIPT TO SAVE METADATA
  #print.data <- data.frame(Species = species.name, Length = max.length, Density_100m2 = mean.density*100,
  #                        Mean_HR = mean.distance,SD_HR = sd.distance, HR_Method = method, Max_age = max.age, Age_method = age.method)
  #spec.data[s,] <- as.matrix(print.data)
  #print(print.data)
 
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
  
    output$scenario_name <- save.name
    returndat[[v]] <- output
 
}

names(returndat) <- disp.val.names
#returndat$scenario_name <- scenario.name

# create and assign folders to save results
if(!dir.exists(paste0(getwd(),'/',results.folder))){ #species.name, "_rBRUVcatchment_",rBRUVcatchment,"/", speciesToAnalyse[s]))){
  dir.create(paste0(getwd(),'/',results.folder)) #"/rBRUVcatchment_",rBRUVcatchment,"/", speciesToAnalyse[s]))}
}

#spec.data <- setNames(as.data.frame(spec.data),names(print.data))
#save.image(paste0(getwd(),'/',results.folder, modelversion,'_', scenario.name,'.RData'))
saveRDS(returndat,paste0(getwd(),'/',results.folder, modelversion,'_', scenario.name,'_returndat.RDS'))

