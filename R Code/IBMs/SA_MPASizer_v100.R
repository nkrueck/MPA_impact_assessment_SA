#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# MPA size decision support model for SA case study by Nils Krueck, Jan 2022          @
# Version 1.00 (Nils Krueck), 08 Jan 2022                                             @
#   - adapted from the generic MPAsizer model version 1.02, 31 July 2020              @
# Version 1.01 (Name), Date,                                                          @
#   - added...                                                                        @
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

SA_MPASizer_v100 <- function(species.name,   # name of species to simulate
                             location.names, # names of locations represented by different density values (optional)
                             MPAsizes,       # vector of MPA sizes to be simulated for conservation impact testing
                             num.data,       # number of individuals / abundance data - define format 
                             mov.data,       # movement data consisting of a list that contains individual movement probability kernels inspecified modelling resolution  
                             rBRUVcatchment, # Radius in m of the plume catchment area around BRUVs (area of attraction)
                             fmort,          # fishing mortality rate(s) when fully exposed (discrete proportion per year)
                             max.age,        # maximum age of species
                             nreplicates, #= 1000, # number of replicate simulations to run
                             resolution = 100, # resolution of the modelling environment in m
                             mean.extent = 1) # specify extent of modeling environment based on either mean (1) or max movement distances (0)


{

  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # Uncomment to run as script #
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  # load data
  #rm(list = ls()) # delete parameters in the workspace
  #modelversion <- "sa_mpasizer100"
  #modelname <- 'SA_MPASizer_v100.R'
  #resolution <- 100
  #scenario.name <- paste0('SA_Res', resolution)
  #maindir <- 'C:/Users/nkrueck/Dropbox/Projects/MPA SA project with Flinders Uni/SA_MPA_model/'
  #data.folder <- 'Github Dec 2021/Data/'
  #results.folder <- paste0('Results/Res', resolution,'/') # specify results folder
  #datadir <- paste0(maindir,data.folder)
  #resultsdir <- paste0(maindir,results.folder)
  
  # load data
  #num.data.all <- readRDS(paste0(datadir,'MaxNs with absences.RDS'))
  #mov.data.all <- readRDS(paste0(datadir,'IndividualMovementProbabilities_Measured_',resolution,'m.RDS'))
  #mov.data <- mov.dat.all[[1]]
  #species.name <- names(mov.data)[1]  # name of species to simulate
  #mov.data <- mov.data$Weekly
  #num.data <- num.data[[1]]$maxn
  
  # Uncomment to run as script
  # MPAsizes <- c(200,1000,5000,10000) #c(seq(100,1000,100),seq(2000,5000,1000),10000)      # vector of MPA sizes to be simulated for conservation impact testing 
  # location.names <- 'SA' #NA     # names of locations represented by different density values (optional)  
  # rBRUVcatchment <- 200
  # fmort <- seq(0,0.3,0.1)  # fishing mortality rate when fully exposed (discrete proportion per year) 
  # max.age <- 15       # maximum age of species
  # nreplicates <- 1000 # number of replicate simulations to run
  # resolution <- 100 # resolution of the modelling environment in m
  # mean.extent <- 1 # specify extent of modeling environment based on either mean (1) or max movement distances (0)  

  
  #@@@@@@@@@@@@@@@@@@@@@@@
  # DATA PREP & CHECK ####
  #@@@@@@@@@@@@@@@@@@@@@@@
  
  #location.names <- 'South_Australia'
  if(is.na(location.names)){
  location.names <- as.character(seq(1,length(unique(num.data$location)),1))}
  nlocations <- length(location.names)
  
  # check mortality
  if(length(which(fmort==0))==0){error("first mortality value should be zero.")}
  
  # get distances
  distances <- unlist(lapply(mov.data,dim))
  distances <- distances[match(paste0(names(mov.data),'1'),names(distances))]
  names(distances) <- names(mov.data)
  
  # sort mpa size
  MPAsizes <- sort(MPAsizes)
  
  # allocate output and working matrices
  sprobslist <- vector("list")
  probslist <- vector("list")
  probsdflist <- vector("list")
  
  ioutputHeader <- c("location", "replicate", "mpa_size", "home_range", "protection")
  ioutput = matrix(nrow=0,ncol=length(ioutputHeader)) # output table for individual data 
  statstableHeader <- c("location", "mean_n", "sd_n", "mean_hr","sd_hr", 
                        "mpa_size","coastline_m","ncells_simulated", "resolution", 
                        "mean_nfull", "mean_npart","mean_prot", "sd_prot", "median_prot")
  statstable = matrix(nrow=nlocations*length(MPAsizes),ncol=length(statstableHeader)) # allocate output table for individual data 
  statsmatsHeader = c("location", "mpa_size",as.character(fmort))
  
  meanmeanp_survival = matrix(data=0,nrow=nlocations*length(MPAsizes),ncol=length(statsmatsHeader))
  
  minminp_survival = meanmeanp_survival
  maxmaxp_survival = meanmeanp_survival
  meanminp_survival = meanmeanp_survival
  meanmaxp_survival = meanmeanp_survival
  p05meanp_survival = meanmeanp_survival
  p95meanp_survival = meanmeanp_survival
  minmeanp_survival = meanmeanp_survival
  maxmeanp_survival = meanmeanp_survival
  meansdp_survival = meanmeanp_survival
  meanp05p_survival = meanmeanp_survival
  meanp95p_survival = meanmeanp_survival
  
  meansumn_survival = meanmeanp_survival
  mediansumn_survival = meanmeanp_survival
  sdsumn_survival = meanmeanp_survival
  p05sumn_survival = meanmeanp_survival
  p95sumn_survival = meanmeanp_survival
  
  meanpropn_mortality = meanmeanp_survival
  medianpropn_mortality = meanmeanp_survival
  sdpropn_mortality = meanmeanp_survival
  p05propn_mortality = meanmeanp_survival
  p95propn_mortality = meanmeanp_survival
  
  mean_mortality_offset = matrix(nrow=nlocations*length(MPAsizes),ncol=length(statsmatsHeader))
  sd_mortality_offset = mean_mortality_offset
  mean_mortalityInst_offset = sd_mortality_offset
  sd_mortalityInst_offset = sd_mortality_offset
  
  
  #@@@@@@@@@@@@@@@@@@@@
  # START OF MODEL ####
  #@@@@@@@@@@@@@@@@@@@@
  
  # count for indexing purposes
  count = 0
  
  starttime <- Sys.time() # save time
  
  
  #@@@@@@@@@@@@@@@
  # location loop (# coded as 1, 2, 3 in vector location_code for low, medium, high abundance)
  #@@@@@@@@@@@@@@@
  
  for (loc in 1:nlocations){
    
    locname = location.names[loc]
    loc.num.data <- num.data[which(num.data$location == locname),]
    loc.num.data <- as.numeric(loc.num.data$maxn)  #[which(specdata$location_code == loc),] # | specdata$location_code == nlocations+1),]
   
    #@@@@@@@@@@@@@@@@
    # MPA size loop #
    #@@@@@@@@@@@@@@@@
    
    for (mpas in 1:length(MPAsizes)){
      
      # count for indexing 
      count = count + 1
      
      # allocate data matrix
      mpastats = matrix(nrow=length(MPAsizes),ncol=5)
      
      # determine mpa size
      mpasize = MPAsizes[mpas]/resolution
      #print(mpasize)
      
      # prepare modelling environment
      maxdist <- max(distances) # maximum travel distance by any individual 
      meandist <- mean(distances) # maximum travel distance by any individual 
      sddist <- sd(distances) # maximum travel distance by any individual 
      clength <- round((mpasize + 2 * maxdist)) # actual extent of simulated coastline
      if (mean.extent == 1) {clength <- round(mpasize + 2 * mean(distances))} # actual length of coastline
      
      coastline <- seq(1,clength+1,1) # index for hypothetical coastline in defined resolution
      mpalocs <- floor(seq(length(coastline)/2-round((mpasize)/2),
                          length(coastline)/2+round((mpasize)/2),1)) # index of mpa locations
      if (mpasize == 0) {mpalocs <- NA}
      
      # specify BRUV locations to determine abundance
      #ndrops = 2*ceiling((clength/(rBRUVcatchment/resolution))/2)+1 # ensure sufficent hypothetical BRUV samples - odd number
      fishlocs = seq(1,length(coastline),rBRUVcatchment/resolution) # BRUV locations
      nfishlocs <- length(fishlocs)
      
      # allocate data matrices
      snfish <- matrix(nrow = nreplicates,ncol=1)
      nfull <- snfish
      npart <- snfish
      meanp_survival <- matrix(ncol=length(fmort),nrow=nreplicates)
      sdp_survival <- meanp_survival
      p95p_survival <- meanp_survival
      p05p_survival <- meanp_survival
      sumn_survival <- meanp_survival
      propn_mortality <- meanp_survival 
      minp_survival <- meanp_survival
      maxp_survival <- meanp_survival
      repdata <- matrix(nrow = nreplicates, ncol= 3)
      pidata <- matrix(nrow=0,ncol=2) # allocate matrix for individual protection across all replicates
    
      #@@@@@@@@@@@@@@@@@
      # Replicate loop #
      #@@@@@@@@@@@@@@@@@
      
      for (rep in 1:nreplicates){
  
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        # Hypothetical fish surveys ####
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        
        nfish <- 0 # set to zero before each sampling
        
        # allocate matrix for individual data (idata) by location
        hpidata = matrix(nrow=0,ncol=2) # allocate matrix for individual home range and protection
        sfishlocs = matrix(nrow=0,ncol=2)
        
        #@@@@@@@@@@@@@@@@@@
        # fish group loop #
        #@@@@@@@@@@@@@@@@@@
        
        for (floc in 1:nfishlocs) {
          
          fishloc <- fishlocs[floc] # index of fish/group location
          #fishAtLoc = nfish.groups[floc] # fish abundance at location
          fishAtLoc <- loc.num.data[round(runif(1,1,length(loc.num.data)))] # randomly selects observations
          nfish <- nfish + fishAtLoc # keep track of total fish number per simulation
          
            if(fishAtLoc > 0){
            
              # determine home range extent for each individual in group at location
              ilocs <- sample(length(distances),fishAtLoc,replace=TRUE)
              idists <- distances[ilocs]
                
              #@@@@@@@@@@@@@@@@@@
              # individual loop #
              #@@@@@@@@@@@@@@@@@@
              
              idata = matrix(nrow=fishAtLoc,ncol=2) # allocate matrix for individual home range and protection
              
              # determine exposure of each individual
              for (ifloc in 1:fishAtLoc){
                
                  hr = idists[ifloc] # extent of home range of individual
                  
                  hrprobs <- mov.data[[ilocs[ifloc]]]$Probability_Occurrence 
                  hrlocs <- seq(1,hr,1) # assign home range vector of cell ids
                  randloc <- sample(hrlocs,1,replace=TRUE,prob=hrprobs) # pick a cell where fish/group is according to probability 
                  
                  idistribution <- seq(fishloc-randloc,fishloc+hr-randloc,1) # assign individual distribution / activity space (overlapping visual census site) along modelling environment
                  
                  # save home range and activity outside of MPA
                  idata[ifloc,1] = hr*resolution
                  idata[ifloc,2] = sum(hrprobs[hrlocs[is.element(idistribution,mpalocs)]],na.rm = TRUE)
                  
              } # end of individual fish loop (ifloc)
              
              hpidata = rbind(hpidata,idata) # store individual data
              
            } # end of abundance > 0 condition 
          
          } # end of fish location loop (floc)
        
          #@@@@@@@@@@@@@@@@@@@@@@@@@@
          # Save replicate data #####
          #@@@@@@@@@@@@@@@@@@@@@@@@@@
    
          snfish[rep] <- sum(nfish)
        
          # save all individual results in big table: c("location", "mpa_size", "replicate", "home range", "protection")
          pidata <- rbind(pidata,hpidata)
          #itable <- cbind(matrix(c(loc,mpasize,rep),nrow = snfish[rep], ncol=3,byrow=TRUE), hpidata)
          
          # save stats per replicate
          #print(hpidata)
          #print(sum(as.numeric(hpidata[,2])==1))
          #print(sum(as.numeric(hpidata[,2])>0))
          
          nfull[rep] <- sum(as.numeric(hpidata[,2])==1) # number of fully protected individuals
          npart[rep] <- sum(as.numeric(hpidata[,2])>0) # number of partially protected individuals
          repdata[rep,] <- c(snfish[rep],nfull[rep],npart[rep]) 
          
          # quantify chances of survival and mortality of fully and partially protected individuals
          exposure <- (1-hpidata[hpidata[,2]>0,2]) # of all at least partially protected individuals
          exposure[exposure<0] <- 0 # exclude numerical diffusion errors
          
          exposuremat <- do.call(cbind, replicate(length(fmort), exposure, simplify=FALSE)) 
          pmortmat <- do.call(rbind, replicate(length(exposure), fmort, simplify=FALSE))
          pstable <- (1 - exposuremat * pmortmat) # probability of survival 
          # creates matrix with individual ages, representing survivors according to mortality risk
          stables <- matrix(rbinom(1:length(pstable),
                                  size = 1, prob = pstable),dim(pstable)[1],dim(pstable)[2]) # survival drawn from probability
          if(max.age>1){ # need to enable option for max age == 1
          for (addyears in 2:max.age) {
            stables <- stables+matrix(rbinom(1:length(pstable),
                                            size = 1, prob = pstable),dim(pstable)[1],dim(pstable)[2]) # survival drawn from probability
            }    
          }
          
          stable <- matrix(as.integer(stables==max.age),dim(pstable)[1],dim(pstable)[2]) 
          
          meanp_survival[rep,] <- apply(pstable,2,mean)
          sdp_survival[rep,] <- apply(pstable,2,sd)
          p95p_survival[rep,] <- apply(pstable,2,quantile,0.95)
          p05p_survival[rep,] <- apply(pstable,2,quantile,0.05)
          minp_survival[rep,] <- apply(pstable,2,min)
          maxp_survival[rep,] <- apply(pstable,2,max)
          sumn_survival[rep,] <- apply(stable,2,sum)
          propn_mortality[rep,] <- 1-sumn_survival[rep,]/sumn_survival[rep,which(fmort==0)] 
          #print(sumn_survival)
          
        } # end of replicate loop (rep)
        
      #reptable = rbind(reptable,itable)
      
      # calculate stats across replicates
      mnfish <- mean(snfish)
      sdnfish <- sd(snfish)
      mnfull <- mean(nfull[npart>0]) # only count fish overlapping the MPA
      mnpart <- mean(npart[npart>0]) # only count fish overlapping the MPA
      mprot <- mean(pidata[pidata[,2]>0,2])
      sdprot <- sd(pidata[pidata[,2]>0,2])
      medprot <- median(pidata[pidata[,2]>0,2])
      mpastats[mpas,] <- c(mnfull,mnpart,mprot,sdprot,medprot)
      
      meanmeanp_survival[count,] = c(loc,mpasize*resolution,apply(meanp_survival,2,mean,na.rm=TRUE))
      minminp_survival[count,] = c(loc,mpasize*resolution,apply(minp_survival,2,min,na.rm=TRUE))
      maxmaxp_survival[count,] = c(loc,mpasize*resolution,apply(maxp_survival,2,max,na.rm=TRUE))
      meanminp_survival[count,] = c(loc,mpasize*resolution,apply(meanp_survival,2,min,na.rm=TRUE))
      meanmaxp_survival[count,] = c(loc,mpasize*resolution,apply(meanp_survival,2,max,na.rm=TRUE))
      p05meanp_survival[count,] = c(loc,mpasize*resolution,apply(meanp_survival,2,quantile,0.05,na.rm=TRUE))
      p95meanp_survival[count,] = c(loc,mpasize*resolution,apply(meanp_survival,2,quantile,0.95,na.rm=TRUE))
      minmeanp_survival[count,] = c(loc,mpasize*resolution,apply(meanp_survival,2,min,na.rm=TRUE))
      maxmeanp_survival[count,] = c(loc,mpasize*resolution,apply(meanp_survival,2,max,na.rm=TRUE))
      meansdp_survival[count,] = c(loc,mpasize*resolution,apply(sdp_survival,2,mean,na.rm=TRUE))
      meanp05p_survival[count,] = c(loc,mpasize*resolution,apply(p05p_survival,2,mean,na.rm=TRUE))
      meanp95p_survival[count,] = c(loc,mpasize*resolution,apply(p95p_survival,2,mean,na.rm=TRUE))
      
      meansumn_survival[count,] = c(loc,mpasize*resolution,apply(sumn_survival,2,mean,na.rm=TRUE))
      mediansumn_survival[count,] = c(loc,mpasize*resolution,apply(sumn_survival,2,median,na.rm=TRUE))
      sdsumn_survival[count,] = c(loc,mpasize*resolution,apply(sumn_survival,2,sd,na.rm=TRUE))
      p05sumn_survival[count,] = c(loc,mpasize*resolution,apply(sumn_survival,2,quantile,0.05,na.rm=TRUE))
      p95sumn_survival[count,] = c(loc,mpasize*resolution,apply(sumn_survival,2,quantile,0.95,na.rm=TRUE))
      
      meanpropn_mortality[count,] = c(loc,mpasize*resolution,apply(propn_mortality,2,mean,na.rm=TRUE))
      medianpropn_mortality[count,] = c(loc,mpasize*resolution,apply(propn_mortality,2,median,na.rm=TRUE))
      sdpropn_mortality[count,] = c(loc,mpasize*resolution,apply(propn_mortality,2,sd,na.rm=TRUE))
      p05propn_mortality[count,] = c(loc,mpasize*resolution,apply(propn_mortality,2,quantile,0.05,na.rm=TRUE))
      p95propn_mortality[count,] = c(loc,mpasize*resolution,apply(propn_mortality,2,quantile,0.95,na.rm=TRUE))
      
      statstable[count,] = c(loc,round(mnfish),round(sdnfish),round(meandist)*resolution,round(sddist)*resolution,
                             mpasize*resolution,clength*resolution,clength,resolution,mpastats[mpas,])
      
      moffset = (1-mprot)*fmort
      sdoffset = (1-sdprot)*fmort
      mean_mortality_offset[count,] = c(loc,mpasize*resolution,moffset)
      sd_mortality_offset[count,] = c(loc,mpasize*resolution,sdoffset)
      mean_mortalityInst_offset[count,] = c(loc,mpasize*resolution,-log(1-moffset))
      sd_mortalityInst_offset[count,] = c(loc,mpasize*resolution,-log(1-sdoffset))
      
      # computing time status report
      mpatime <- Sys.time()
      print(paste0("Location ", loc, "/" , nlocations, ": ", "MPA size ", mpas, "/" , length(MPAsizes) , " completed."))
      print(mpatime - starttime)  
      
    } # mpa loop   
    
  } # location loop
  

  endtime <- Sys.time()
  print(endtime - starttime)
  
  # convert results to data frames
    #ioutput <- as.data.frame(ioutput)
    #colnames(ioutput) = ioutputHeader
    
    statstable <- as.data.frame(statstable)
    colnames(statstable) = statstableHeader
    
    mean_survival_probability = as.data.frame(meanmeanp_survival, row.names = FALSE) 
    colnames(mean_survival_probability) = statsmatsHeader
    
    sd_survival_probability = as.data.frame(meansdp_survival, row.names = FALSE) 
    colnames(sd_survival_probability) = statsmatsHeader
    
    mean_protected_individuals = as.data.frame(meansumn_survival, row.names = FALSE)
    colnames(mean_protected_individuals) = statsmatsHeader
    #mean_protected_density <- mean_protected_individuals[,c(1,3:dim(mean_protected_individuals)[2])]
    #mean_protected_density[,2:dim(mean_protected_density)[2]] <- mean_protected_density[,2:dim(mean_protected_density)[2]]/statstable$coastline_m*100
        
    median_protected_individuals = as.data.frame(mediansumn_survival, row.names = FALSE)
    colnames(median_protected_individuals) = statsmatsHeader
    
    sd_protected_individuals = as.data.frame(sdsumn_survival, row.names = FALSE)
    colnames(sd_protected_individuals) = statsmatsHeader
    
    mean_Fmortality = as.data.frame(meanpropn_mortality, row.names = FALSE)
    colnames(mean_Fmortality) = statsmatsHeader
    
    median_Fmortality = as.data.frame(medianpropn_mortality, row.names = FALSE)
    colnames(median_Fmortality) = statsmatsHeader
    
    sd_Fmortality = as.data.frame(sdpropn_mortality, row.names = FALSE)
    colnames(sd_Fmortality) = statsmatsHeader
    
    mean_Foffset = as.data.frame(mean_mortality_offset, row.names = FALSE)
    colnames(mean_Foffset) = statsmatsHeader
    
    sd_Foffset = as.data.frame(sd_mortality_offset, row.names = FALSE)
    colnames(sd_Foffset) = statsmatsHeader
    
    other_parameters <-  as.data.frame(t(c(nreplicates,max.age,rBRUVcatchment,mean.extent,resolution)))
    colnames(other_parameters) <- c("num_replicates", "lifetime_years","radius_BRUV_catchment", "mean_extent_chosen","resolution")
 
    #return results
    returndat <- list(statstable = statstable,
                    mean_survival_probability = mean_survival_probability,
                    sd_survival_probability = sd_survival_probability,
                    mean_protected_individuals = mean_protected_individuals,
                    median_protected_individuals = median_protected_individuals,
                    sd_protected_individuals = sd_protected_individuals,
                    mean_Fmortality = mean_Fmortality,
                    median_Fmortality = median_Fmortality,
                    sd_Fmortality = sd_Fmortality,
                    mean_Foffset = mean_Foffset,
                    sd_Foffset = sd_Foffset,
                    other_parameters = other_parameters)
  
  return(returndat)
  

}# Function end


