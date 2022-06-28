#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# MPA size decision support model for SA case study by Nils Krueck, Jan 2022          @
# Script to analyse & plot results from model v1.00 (Nils Krueck), 20 Feb 2022        @
#   - version analyse_SA_MPA_sizer_v100_nils - adapted from the first generic srcipt  @
#   - folders specific to Nils' computer                                              @
#                                                                                     @
# Version:                                                                            @
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

rm(list = ls()) # delete parameters in the workspace
library(stringr)
library(dplyr)
library(matrixStats)

# DUSKY WHALER MAXN DATA - JUST ZEROS!!!!!!!!!!

# Specify model and parameterisation
modelid <- 100 # specify model version / id
resolution <- 100 # modelling resolution in m
rBRUVcatchment <- 100 # radius of BRUV catchment
dispersal.period <- 'Weekly' # specify dispersal period
age <- "Max" # specify age threshold as either "Max" or "Maturity"
case.study.region <- "SA" # specify region name in case location names are not specified or pooled
location.names <- case.study.region # set to 'NA' unless data is supposed to be pooled!!!!  
MPAsizes <- c(0,seq(500,2000,500),seq(3000,5000,1000),seq(10000,30000,10000),50000,100000) # specify MPA sizes to be analysed (in m)
#fmort <- seq(0,0.5,.1)  # fishing mortality rate when fully exposed (discrete proportion per year) 
#nreplicates <- 1000 # number of replicate simulations to run
#mean.extent <- 1 # specify extent of modeling environment based on either mean (1) or max movement distances (0)  
fmortplot <- c(0.05,0.2,0.5) # lifetime mortality plots
foffplot <- c(1) # annual mortality offset/risk/survival plots (all the same so need to pick one level) 

idealized.maxn <- 1
if(idealized.maxn == 1){case.study.region <- paste0(case.study.region,'_idealized')}
mean.dist <- 1
if(mean.dist == 1){case.study.region <- paste0(case.study.region,'_meanDist')}


# plot specs
lwd <- 2 # line width
res <- 300 # resolution of plot (pixels)
xlabel.size <- "Protected area size (km)"
xlabel.area <- "Protected area size (km^2)"

# Load model and data
modelversion <- paste0("SA_MPAsizer",modelid)
modelname <- paste0('SA_MPASizer_v',modelid,'.R')
scenario.name <- paste0(case.study.region,'_Res', resolution,'_',dispersal.period,'_AgeAt',age,'_rBRUV',rBRUVcatchment)
maindir <- paste0(getwd(),'/')

#setwd(datadir)

# load data
datadir <- paste0(maindir, 'Data/')
num.data.all <- readRDS(paste0(datadir,'MaxNs with absences.RDS'))
mov.data.all <- readRDS(paste0(datadir,'Movement_Profiles/IndividualMovementProbabilities_Measured_',resolution,'m.RDS'))
# species.names <- names(mov.dat.all)
# nspecies <- length(species.names)
age.data <- read.csv(paste0(datadir,'input_parameters_updated.csv'))

# specify data folders and load data
results.folder <- 'Results/'
resultsdir <- paste0(maindir,results.folder,modelversion,'/')
returndat <- readRDS(paste0(resultsdir,modelversion,'_',scenario.name,'_returndat.RDS'))

scenario.results.folder <- paste0(resultsdir,scenario.name,'/')
if(!dir.exists(paste0(scenario.results.folder))){ #species.name, "_rBRUVcatchment_",rBRUVcatchment,"/", speciesToAnalyse[s]))){
  dir.create(paste0(scenario.results.folder)) #"/rBRUVcatchment_",rBRUVcatchment,"/", speciesToAnalyse[s]))}
}

nona.locs <- which(!is.na(returndat)&!is.null(returndat))
returndat <- returndat[nona.locs]
nspecies <- length(returndat)
species.names <- names(returndat) 
species.groups <- age.data$species.group[match(species.names,age.data$common.name)]
species.order <- order(species.groups,species.names) #species.groups,species.names) # determine order of species to plot
species.names.ordered <- species.names[species.order] # get names of species in order
ngroups <- xtabs(~species.groups)

returndat <- returndat[species.names.ordered]

# subselect data by matching MPA sizes (if needed) 
for(rd in 1:length(returndat)){
  rdat <- returndat[[rd]]
  listids <- which(is.na(match(names(rdat),c("other_parameters","scenario_name"))))
  sizeids <- which(!is.na(match(rdat$statstable$mpa_size,MPAsizes)))
  for(rd2 in listids){rdat[[rd2]] <- rdat[[rd2]][sizeids,]} # replace content that matches MPA sizes
  returndat[[rd]] <- rdat # replace list
}

# specify colours and line types to use for different species/groups
ucols <- c("blue","grey50","black"); ltys <- NA; cols <- NA
for(g in 1:length(ngroups)){
  ltys <- c(na.omit(ltys),1:ngroups[g])
  cols <- c(na.omit(cols),rep(ucols[g],ngroups[g]))} # specify colours
#cols <- cols[species.order]
#ltys <- ltys[species.order]


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# extract MPA sizes needed to achieve certain protection levels #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

mean_iprot <- matrix(nrow=nspecies,ncol=length(MPAsizes))
mean_prot <- mean_iprot
median_prot <- mean_iprot
mean_hr <- matrix()
mean_nfish <- mean_iprot
mean_nfull <- mean_iprot
mean_npart <- mean_iprot
#floc <- 2 #which(fmort==1)+2

for (s in 1:nspecies) {
  mean_iprot[s,] <- returndat[[s]]$statstable$mean_nfull/returndat[[s]]$statstable$mean_npart #mean_protected_individuals[,3]#returndat[[s]]$mean_protected_individuals[,floc]/
                    #returndat[[s]]$mean_protected_individuals[,3] #/returndat[[s]]$statstable$mean_nfish
  mean_nfull[s,] <- returndat[[s]]$statstable$mean_nfull
  mean_npart[s,] <- returndat[[s]]$statstable$mean_npart
  mean_prot[s,] <- returndat[[s]]$statstable$mean_prot
  median_prot[s,] <- returndat[[s]]$statstable$median_prot
  mean_hr[s] <- unique(returndat[[s]]$statstable$mean_hr)
  #mean_nfish[s,] <- returndat[[s]]$statstable$mean_nfish#/(returndat[[s]]$statstable$mpa_size+2*mean_hr[s])
}

mean_prot[is.na(mean_prot)]<-1
mean_protf0 <- mean_prot[,1]
mean_prot[,1] <- rep(0,nspecies)
#plot(rep(MPAsizes,nspecies),mean_prot)

if(idealized.maxn == 1) { # if abudance increased
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@
  # Plot mean protection ####
  #@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  # specify effective range
  min.effective.size <- NA #MPAsizes[min(which(colMeans(mean_prot)>0.1))]/1000
  max.effective.size <- NA #MPAsizes[min(which(colMeans(mean_prot)>=0.5))]/1000
  
  xlims <- c(10,100)
  
  for(xl in 1:length(xlims)){
    
    xlim <- xlims[xl]
    
    tiff(paste0(scenario.results.folder, scenario.name,'_Protection',xlim,'km.tiff'), units="cm", width=18, height=12, res=res)
  
    par(mar=c(5.1, 4.1, 4.1, 12.1)) # this is usually the default
  
    plot(t(MPAsizes/1000),mean_prot[1,],type='l',xlim = c(0,xlim),ylim = c(0,1),
         xlab = xlabel.size, ylab = "Time spent in protected area",
         lty=ltys[s],col=cols[s],lwd=lwd)
    axis(side = 1, at = seq(0,xlim,xlim/10))
    
    # add effective area
    polygon(c(min.effective.size,min.effective.size,max.effective.size,max.effective.size),c(0,1,1,0),col='light green',border="light green")
    
    for (s in 1:nspecies) {
      lines(t(MPAsizes/1000),mean_prot[s,],lty=ltys[s],col=cols[s],lwd=lwd)
    }
    #abline(v=min.effective.size,lty=1,lw=2,col='red')
    #abline(v=max.effective.size,lty=1,lw=2,col='red')
    
    #crit.species <- species.names[which(mean_prot[,4]<0.5)]
    #text(8,0.4,crit.species[1])
    #text(8,0.3,crit.species[2])
    
    # add legend
    par(xpd=TRUE)
    legend("topright", inset=c(-0.6,0), legend=str_replace_all(species.names.ordered,'_',' '), 
           lty=ltys[1:nspecies],col=cols[1:nspecies],lwd=2) #, title="Group")
    
    dev.off()
    par(xpd=FALSE)
  }
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # Calculate MPA size to achieve certain protection levels ##### COULD do by group!!!!!!!!
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  prot.levels <- c(0.25,0.5,0.75,0.95)
  sizes <- matrix(ncol = length(prot.levels),nrow = nspecies)
  cnspeciesprotected <- matrix(nrow = length(prot.levels), ncol = length(MPAsizes))
  for (p in 1:length(prot.levels)){
    sizes[,p] <- MPAsizes[apply(mean_iprot,1,function(x) min(which(x>prot.levels[p])))]
    minmpasizemat <- matrix(rep(sizes[,p],length(MPAsizes)),nrow=nspecies,ncol=length(MPAsizes))
    mpasizemat <- t(matrix(rep(MPAsizes,nspecies),ncol=nspecies,nrow=length(MPAsizes)))
    speciesprotmat <- minmpasizemat <= mpasizemat
    cnspeciesprotected[p,] <- colSums(speciesprotmat,na.rm=TRUE)
  }
  
  pspeciesprotected <- cnspeciesprotected/nspecies*100
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # Plot species protection (four protection levels) ####
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # two different resolutions
  
  xlim <- 10
  
  legendlabels <- paste0(prot.levels * 100, '%')
  
  tiff(paste0(scenario.results.folder, scenario.name,'_Species_protection_',xlim,'km.tiff'), units="cm", width=18, height=12, res=res)
  
    par(mar=c(5.1, 4.1, 4.1, 12.1)) # this is usually the default
  
    plot(MPAsizes/1000, pspeciesprotected[1,],type='l',lty=ltys[1],lwd=lwd,xlim = c(0,xlim),
         xlab = xlabel.size, ylab = "Species protected (%)")
    axis(side = 1, at = seq(0,xlim,xlim/10))
    
    for (p in 2:length(prot.levels)) {
      lines(MPAsizes/1000, pspeciesprotected[p,],type='l',lty=ltys[p],lwd=lwd)
    }
    #abline(v=.5,lty=1,lw=2,col='red')
    #abline(v=4,lty=1,lw=2,col='red')
    par(xpd=TRUE)
    legend('topright',legendlabels,inset=c(-0.6,0),lty=ltys,lwd=lwd,title="Individual protection level")
    dev.off()
    par(xpd=FALSE)
    
    
    
    xlim <- 100
  
    tiff(paste0(scenario.results.folder, scenario.name,'_Species_protection_',xlim,'km.tiff'), units="cm", width=18, height=12, res=res)
    
    par(mar=c(5.1, 4.1, 4.1, 12.1)) # this is usually the default
    
    
    plot(MPAsizes/1000, pspeciesprotected[1,],type='l',lty=ltys[1],lwd=lwd,xlim = c(0,xlim),
         xlab = xlabel.size, ylab = "Species protected (%)")
    axis(side = 1, at = seq(0,xlim,xlim/10))
    for (p in 2:length(prot.levels)) {
      lines(MPAsizes/1000, pspeciesprotected[p,],type='l',lty=ltys[p],lwd=lwd)
    }
    
    par(xpd=TRUE)
    legend('topright',legendlabels,inset=c(-0.6,0),lty=ltys,lwd=lwd,title="Individual protection level")
    dev.off()
    par(xpd=FALSE)
    
    #abline(v=.5,lty=1,lw=2,col='red')
    #abline(v=4,lty=1,lw=2,col='red')
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # Plot lifetime fishing mortality and fishing mortality offsets ####
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  #xlim <- 10
  ylim <- 100
  #ltys <- 1
  #lwds <- 1
  
  if(age == 'Max'){s.ylabel <- 'Lifetime survival (%)'}
  if(age == 'Maturity'){s.ylabel <- 'Survival to reproduction (%)'} # Survival to reproduction

  fmortcols <- match(fmortplot,names(returndat[[1]]$mean_Fmortality))
  lfmorts <- array(dim=c(nspecies,length(MPAsizes),length(fmortplot)))
  foffs <- lfmorts
  fsurvs <- lfmorts
  
  xlims <- c(10,100)
  
  for(xl in 1:length(xlims)){
    
    xlim <- xlims[xl]
    
    # lifetime fishing mortality
    for (plotnum in 1:length(fmortcols)){
      lfmorts[1,,plotnum] <- 1-returndat[[1]]$mean_Fmortality[,fmortcols[plotnum]]
    
      tiff(paste0(scenario.results.folder, scenario.name,'_Lifetime_survival_F',
                  fmortplot[plotnum]*100,'_',xlim,'km.tiff'), units="cm", width=18, height=12, res=res)
      
      par(mar=c(5.1, 4.1, 4.1, 12.1)) # this is usually the default
      
      plot(MPAsizes/1000,lfmorts[1,,plotnum] * 100,
           type='l',col=cols[species.order[1]],lty=ltys[species.order[1]],lwd=lwd,xlim = c(0,xlim),ylim = c(0,ylim),
           xlab = xlabel.size, ylab = s.ylabel,
           main = paste0('Age at ',age,': ',fmortplot[plotnum] * 100, '% annual fishing mortality'))
      axis(side = 1, at = seq(0,xlim,xlim/10))
      if(length(returndat)>1){
        for (s in 1:nspecies) {
         lfmorts[s,,plotnum] <- 1-returndat[[s]]$mean_Fmortality[,fmortcols[plotnum]]
         lines(MPAsizes/1000, lfmorts[s,,plotnum] * 100,col=cols[s],lty=ltys[s],lwd=lwd)
        }
      }
      par(xpd=TRUE)
      legend("topright", inset=c(-0.6,0), legend=str_replace_all(species.names.ordered,'_',' '), 
             lty=ltys,col=cols,lwd=2) #, title="Group")
      dev.off()
      par(xpd=FALSE)
    }
  
  }
    
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # fishing mortality offsets @
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  foffcols <- match(foffplot,names(returndat[[1]]$mean_Fmortality))
  
  xlims <- c(10,100)
  
  for(xl in 1:length(xlims)){
    
    xlim <- xlims[xl]
    
    for (plotnum in 1:length(foffcols)){
    
      ylim <- ceiling(foffplot[plotnum]*100)/100 #ceiling(-log(fmortplot[plotnum])*10)/10 # get maximum according to specs (turned into instantaneous mortality)
      #print(ylim) 
      
      foffs[1,,plotnum] <- returndat[[1]]$mean_Foffset[,foffcols[plotnum]]
      foffs[1,1,plotnum] <- foffplot[plotnum] # fix first value at specified mortality
      foffs[1,foffs[1,,plotnum]>1,plotnum] <- 1 # cap to mean risk of 1
      tiff(paste0(scenario.results.folder, scenario.name,'_Fmortality_risk_F',foffplot[plotnum]*100,'_',xlim,'km.tiff'), 
           units="cm", width=18, height=12, res=res)
      
      par(mar=c(5.1, 4.1, 4.1, 12.1)) # this is usually the default
      
      
      plot(MPAsizes/1000,foffs[1,,plotnum],
           type='l',lwd=lwd,xlim = c(0,xlim),ylim = c(0,ylim),
           xlab = xlabel.size, ylab = "Fishing mortality offset",
           main = 'Annual fishing mortality', #paste0(foffplot[plotnum]*100, '% Annual fishing mortality'),
           col=cols[1],lty=ltys[1])
      axis(side = 1, at = seq(0,xlim,xlim/10))
    
      if(length(returndat)>1){
       for (s in 1:nspecies) {
        foffs[s,,plotnum] <- returndat[[s]]$mean_Foffset[,foffcols[plotnum]]
        foffs[s,1,plotnum] <- foffplot[plotnum] # fix first value at specified mortality
        foffs[s,foffs[s,,plotnum]>1,plotnum] <- 1 # cap at 1
        lines(MPAsizes/1000, foffs[s,,plotnum],col=cols[s],lty=ltys[s],lwd=lwd)
       }
      }
      par(xpd=TRUE)
      legend("topright", inset=c(-0.6,0), legend=str_replace_all(species.names.ordered,'_',' '), 
             lty=ltys,col=cols,lwd=2) #, title="Group")
      dev.off()
      par(xpd=FALSE)
    }
  
  } # xlim loop
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@
  # Probability of survival @
  #@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  ylim <- 1
  
  xlims <- c(10,100)
  
  for(xl in 1:length(xlims)){
    
    xlim <- xlims[xl]
    
    for (plotnum in 1:length(foffcols)){
      fsurvs[1,,plotnum] <- returndat[[1]]$mean_survival_probability[,foffcols[plotnum]]
      fsurvs[1,1,plotnum] <- 1-foffplot[plotnum]
      
      tiff(paste0(scenario.results.folder, scenario.name,'_Psurvival_F',foffplot[plotnum]*100,'_',xlim,'km.tiff'), 
         units="cm", width=18, height=12, res=200)
    
      par(mar=c(5.1, 4.1, 4.1, 12.1)) # this is usually the default
      
      plot(MPAsizes/1000,fsurvs[1,,plotnum],
           type='l',col=cols[1],lty=ltys[1],lwd=lwd,xlim = c(0,xlim),ylim = c(0,ylim),
           xlab = xlabel.size, ylab = "Annual probability of survival",
           main = paste0(foffplot[plotnum]*100, '% annual fishing mortality'))
      axis(side = 1, at = seq(0,xlim,xlim/10))
      
      if(length(returndat)>1){
        for (s in 1:nspecies) {
          fsurvs[s,,plotnum] <- returndat[[s]]$mean_survival_probability[,foffcols[plotnum]]
          fsurvs[s,1,plotnum] <- 1-foffplot[plotnum]
          lines(MPAsizes/1000, fsurvs[s,,plotnum],col=cols[s],lty=ltys[s],lwd=lwd)
        }
      }
      par(xpd=TRUE)
      legend("topright", inset=c(-0.6,0), legend=str_replace_all(species.names.ordered,'_',' '), 
             lty=ltys,col=cols,lwd=lwd) #, title="Group")
      dev.off()
      par(xpd=FALSE)
    }
  
  } # xlim loop

} # idealized abundance condition  
    
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Plot relative density of protected individuals ####
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

if(idealized.maxn == 0){ # only plot if densities are not idealized

  
  # # specify effective range
  # #min.effective.size <- MPAsizes[min(which(colMeans(mean_prot)>0.1))]/1000
  # #max.effective.size <- MPAsizes[min(which(colMeans(mean_prot)>=0.5))]/1000
  # 
  # protection.type = 'Fully' # or partially
  # response = 'Density'
  # response2 <- "relative"
  # pMPAsizes <- MPAsizes[MPAsizes!=1]
  # 
  # if(protection.type=='Fully'){mean_nprot <- mean_nfull
  # }else{mean_nprot <- mean_nprot}
  # 
  # # individuals per 100 m2
  # if(response == 'Density'){
  # mean_nprot <- mean_nprot/t(matrix(rep(MPAsizes,nspecies),length(MPAsizes),nspecies))*100
  # ylabel <- paste0(protection.type, ' protected individuals (/100 m^2)')
  # } else {ylabel <- paste0(protection.type,' protected individuals')}
  # mean_nprot <- mean_nprot[,match(pMPAsizes,MPAsizes)]
  # 
  # if(response2 == 'relative'){
  #   mean_nprot <- mean_nprot/rowMaxs(mean_nprot,na.rm = TRUE)*100
  #   ylabel <- paste0(protection.type, ' protected individuals (% max density)')}
  # 
  # xlims <- c(10,100)
  # 
  # for(xl in 1:length(xlims)){
  # 
  #   xlim <- xlims[xl]
  #   ylim <- max(mean_nprot[,3:max(which(pMPAsizes<=xlim*1000))],na.rm = TRUE)
  # 
  #   tiff(paste0(scenario.results.folder, scenario.name,'_',protection.type,'_protected_individuals_',xlim,'km.tiff'), units="cm", width=18, height=12, res=res)
  #   
  #     par(mar=c(5.1, 4.1, 4.1, 12.1)) # this is usually the default
  #     
  #     plot(t(pMPAsizes/1000),mean_nprot[1,],type='l',xlim = c(0,xlim),ylim = c(0,ylim),
  #          xlab = xlabel.size, ylab = ylabel,
  #          lty=ltys[s],col=cols[s],lwd=lwd)
  #     axis(side = 1, at = seq(0,xlim,xlim/10))
  #     
  #     # add effective area
  #     #polygon(c(min.effective.size,min.effective.size,max.effective.size,max.effective.size),c(0,1,1,0),col='light green',border="light green")
  #     
  #     for (s in 1:nspecies) {
  #       lines(t(pMPAsizes/1000),mean_nprot[s,],lty=ltys[s],col=cols[s],lwd=lwd)
  #     }
  #     #abline(v=min.effective.size,lty=1,lw=2,col='red')
  #     #abline(v=max.effective.size,lty=1,lw=2,col='red')
  #     
  #     #crit.species <- species.names[which(mean_prot[,4]<0.5)]
  #     #text(8,0.4,crit.species[1])
  #     #text(8,0.3,crit.species[2])
  #     
  #     # add legend
  #     par(xpd=TRUE)
  #     legend("topright", inset=c(-0.6,0), legend=str_replace_all(species.names.ordered,'_',' '), 
  #            lty=ltys[1:nspecies],col=cols[1:nspecies],lwd=2) #, title="Group")
  #   
  #   dev.off()
  #   par(xpd=FALSE)
  # 
  # }


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Plot number of protected individuals ####
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  # calculate area multipliers
  MPAmult <- MPAsizes/rBRUVcatchment

  # specify effective range
  #min.effective.size <- MPAsizes[min(which(colMeans(mean_prot)>0.1))]/1000
  #max.effective.size <- MPAsizes[min(which(colMeans(mean_prot)>=0.5))]/1000
  #fmortplot <- c(0.05,0.10,0.20,0.50)
  responses = c('Density') # 'Density' or 'Number'
  responses2 <- c('Absolute','Relative') # "Relative" or 'Absolute'
  pMPAsizes <- MPAsizes[MPAsizes!=1] # exclude MPA size of 1 m
  xlims <- c(10,100)
  fmortcols <- match(fmortplot,(names(returndat[[1]]$mean_protected_individuals)))
  
  
  for(r1 in 1:length(responses)){
  
    response <- responses[r1]
    
    for(r2 in 1:length(responses2)){

      response2 <- responses2[r2]
      
      for(fp in 1:length(fmortplot)){
        
        mean_nprot <- matrix(rep(NA,length(mean_nfull)),dim(mean_nfull)) #array(NA,c(dim(mean_nfull),length(fmortcols)))
        mean_nprot_area <- mean_nprot
        mean_dens <-mean_nprot
        median_nprot <- mean_nprot
        sd_nprot <- mean_nprot
        logy <- 'y'
        
        for (s in 1:nspecies) {
          mean_nprot[s,] <- returndat[[s]]$mean_protected_individuals[,fmortcols[fp]]
          mean_nprot_area[s,] <- returndat[[s]]$mean_protected_individuals[,fmortcols[fp]]*MPAmult # same as / MPAmult * MPAmult^2 (representing fish protected per )
          mean_dens[s,] <- returndat[[s]]$mean_protected_individuals[,fmortcols[fp]]/MPAsizes*100
          median_nprot[s,] <- returndat[[s]]$median_protected_individuals[,fmortcols[fp]]
          sd_nprot[s,] <- returndat[[s]]$sd_protected_individuals[,fmortcols[fp]]
        }
        #mean_nport <- median_nprot
      
        #if(protection.type=='Fully'){mean_nprot <- mean_nfull
        #}else{mean_nprot <- mean_npart}
        
        xlabel <- xlabel.size
        # individuals per 100 m2
        if(response == 'Numbers_p100m'){
          mean_nprot <- mean_dens#mean_nprot/t(matrix(rep(MPAsizes,nspecies),length(MPAsizes),nspecies))*100
          ylabel <- paste0('Protected individuals (/100 m^2)')
        } else {ylabel <- paste0('Protected individuals')}
        mean_nprot <- mean_nprot[,match(pMPAsizes,MPAsizes)]
        
        if(response == 'Density'){
          mean_nprot <- mean_nprot_area
          ylabel <- paste0('Protected individuals')
          xlabel <- xlabel.area}
        
        if(response2 == 'Relative'){
          mean_nprot <- mean_dens/rowMaxs(mean_dens,na.rm = TRUE)*100
          ylabel <- paste0('Protected individuals (% max)')}
     
        for(xl in 1:length(xlims)){
          
          xlim <- xlims[xl]
          ylim <- max(mean_nprot[,3:max(which(pMPAsizes<=xlim*1000))],na.rm = TRUE)
          
          tiff(paste0(scenario.results.folder, scenario.name,'_',
                      'Protected_individuals_',response,'_',response2,'_F',fmortplot[fp]*100,'_',xlim,'km.tiff'), units="cm", width=18, height=12, res=res)
          
          par(mar=c(5.1, 4.1, 4.1, 12.1)) # this is usually the default
          
          if(response2=='Relative'){
          plot(t(pMPAsizes/1000),mean_nprot[1,],type='l',xlim = c(0,xlim),ylim = c(0,ylim),
               xlab = xlabel, ylab = ylabel,
               lty=ltys[s],col=cols[s],lwd=lwd,
               main = paste0('Age at ',age,': ',fmortplot[fp] * 100, '% annual fishing mortality'))
          }else{
          plot(t(pMPAsizes/1000),mean_nprot[1,],type='l',xlim = c(0,xlim),ylim = c(1,ylim),
               xlab = xlabel, ylab = ylabel,
               lty=ltys[s],col=cols[s],lwd=lwd, log='y',
               main = paste0('Age at ',age,': ',fmortplot[fp] * 100, '% annual fishing mortality'))}
          axis(side = 1, at = seq(0,xlim,xlim/10))
          
          # add effective area
          #polygon(c(min.effective.size,min.effective.size,max.effective.size,max.effective.size),c(0,1,1,0),col='light green',border="light green")
          
          for (s in 1:nspecies) {
            lines(t(pMPAsizes/1000),mean_nprot[s,],lty=ltys[s],col=cols[s],lwd=lwd)
          }
          #abline(v=min.effective.size,lty=1,lw=2,col='red')
          #abline(v=max.effective.size,lty=1,lw=2,col='red')
          
          #crit.species <- species.names[which(mean_prot[,4]<0.5)]
          #text(8,0.4,crit.species[1])
          #text(8,0.3,crit.species[2])
          
          # add legend
          par(xpd=TRUE)
          legend("topright", inset=c(-0.6,0), legend=str_replace_all(species.names.ordered,'_',' '), 
                 lty=ltys[1:nspecies],col=cols[1:nspecies],lwd=2) #, t itle="Group")
          
          dev.off()
          par(xpd=FALSE)
          
        } # MPA size limits
      
      } # mortality levels

    } # response 2
  
  } # response 1
  
    
} # only plot of densities are not idealized