
rm(list = ls()) # delete parameters in the workspace

# Specify model and parameterisation
modelid <- 100 # specify model version / id
resolution <- 100 # modelling resolution in m
rBRUVcatchment <- 200 # radius of BRUV catchment
dispersal.period <- 'Weekly' # specify dispersal period
age <- "Maturity" # specify age threshold as either "Max" or "Maturity"
case.study.region <- "SA" # specify region name in case location names are not specified or pooled
location.names <- case.study.region # set to 'NA' unless data is supposed to be pooled!!!!  
MPAsizes <- c(0,seq(500,2000,500),seq(3000,5000,1000),seq(10000,30000,10000),50000,100000) # specify MPA sizes to be analysed (in m)
fmort <- seq(0,0.5,.1)  # fishing mortality rate when fully exposed (discrete proportion per year) 
nreplicates <- 1000 # number of replicate simulations to run
mean.extent <- 1 # specify extent of modeling environment based on either mean (1) or max movement distances (0)  

# Load model and data
modelversion <- paste0("SA_MPAsizer",modelid)
modelname <- paste0('SA_MPASizer_v',modelid,'.R')
scenario.name <- paste0(case.study.region,'_Res', resolution,'_',dispersal.period,'_AgeAt',age,'_rBRUV',rBRUVcatchment)
maindir <- paste0(getwd(),'/')
data.folder <- 'Data/'

library(stringr)
library(dplyr)
library(doParallel)

# specify data folders and load data
results.folder <- 'Results'
datadir <- paste0(maindir,data.folder)
scenario.results.folder <- paste0(maindir,results.folder,'/',modelversion,'_',scenario.name,'/')

#setwd(datadir)

# load data
# num.data.all <- readRDS(paste0(datadir,'MaxNs with absences.RDS'))
# mov.data.all <- readRDS(paste0(datadir,'IndividualMovementProbabilities_Measured_',resolution,'m.RDS'))
# species.names <- names(mov.dat.all)
# nspecies <- length(species.names)

returndat <- readRDS(paste0(scenario.results.folder,modelversion,'_',scenario.name,'_returndat.RDS'))
nspecies <- length(returndat)
species.names <- names(returndat) 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# extract MPA sizes needed to achieve certain protection levels #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

MPAsizes <- c(0,seq(500,2000,500),seq(3000,5000,1000),seq(10000,30000,10000),50000,100000)
mean_iprot <- matrix(nrow=nspecies,ncol=length(MPAsizes))
mean_prot <- mean_iprot
median_prot <- mean_iprot
mean_hr <- matrix()
mean_nfish <- mean_iprot
mean_nfull <- mean_iprot
floc <- 2 #which(fmort==1)+2

for (s in 1:length(returndat)) {
  mean_iprot[s,] <- returndat[[s]]$statstable$mean_nfull/returndat[[s]]$statstable$mean_npart #mean_protected_individuals[,3]#returndat[[s]]$mean_protected_individuals[,floc]/
                    #returndat[[s]]$mean_protected_individuals[,3] #/returndat[[s]]$statstable$mean_nfish
  mean_nfull[s,] <- returndat[[s]]$statstable$mean_nfull
  mean_prot[s,] <- returndat[[s]]$statstable$mean_prot
  median_prot[s,] <- returndat[[s]]$statstable$median_prot
  mean_hr[s] <- unique(returndat[[s]]$statstable$mean_hr)
  #mean_nfish[s,] <- returndat[[s]]$statstable$mean_nfish#/(returndat[[s]]$statstable$mpa_size+2*mean_hr[s])
}

mean_prot[is.na(mean_prot)]<-1
mean_protf0 <- mean_prot[,1]
mean_prot[,1] <- rep(0,nspecies)
#plot(rep(MPAsizes,nspecies),mean_prot)


# Plot mean protection ####
tiff(paste0(scenario.results.folder, scenario.name,'_Protection.tiff'), units="cm", width=16, height=12, res=200)

plot(t(MPAsizes/1000),mean_prot[1,],type='l',xlim = c(0,10),
     xlab = "MPA size (km)", ylab = "Time spent in MPA")
axis(side = 1, at = seq(0,10,1))
for (s in 1:length(returndat)) {
  lines(t(MPAsizes/1000),mean_prot[s,])
}
abline(v=.5,lty=1,lw=2,col='red')
abline(v=4,lty=1,lw=2,col='red')
crit.species <- species.names[which(mean_prot[,4]<0.5)]
text(8,0.4,crit.species[1])
text(8,0.3,crit.species[2])
dev.off()


# Calculate MPA size to achieve certain protection levels ##### 
prot.levels <- c(0.25,0.5,0.75,0.95)
sizes <- matrix(ncol = length(prot.levels),nrow = s)
cnspeciesprotected <- matrix(nrow = length(prot.levels), ncol = length(MPAsizes))
for (p in 1:length(prot.levels)){
  sizes[,p] <- MPAsizes[apply(mean_iprot,1,function(x) min(which(x>prot.levels[p])))]
  minmpasizemat <- matrix(rep(sizes[,p],length(MPAsizes)),nrow=nspecies,ncol=length(MPAsizes))
  mpasizemat <- t(matrix(rep(MPAsizes,nspecies),ncol=nspecies,nrow=length(MPAsizes)))
  speciesprotmat <- minmpasizemat <= mpasizemat
  cnspeciesprotected[p,] <- colSums(speciesprotmat,na.rm=TRUE)
}

pspeciesprotected <- cnspeciesprotected/nspecies*100

# Plot species protection (four protection levels) ####
xlim <- 10
ltys <- c(1,2,3,4)
lwds <- 2
legendlabels <- paste0(prot.levels * 100, '% individual protection')

tiff(paste0(scenario.results.folder, scenario.name,'_Species_protection_',xlim,'km.tiff'), units="cm", width=16, height=12, res=200)

plot(MPAsizes/1000, pspeciesprotected[1,],type='l',lty=ltys[1],lwd=lwds,xlim = c(0,xlim),
     xlab = "MPA size (km)", ylab = "Species protected (%)")
axis(side = 1, at = seq(0,xlim,xlim/10))

for (p in 2:length(prot.levels)) {
  lines(MPAsizes/1000, pspeciesprotected[p,],type='l',lty=ltys[p],lwd=lwds)
}
legend('bottomright',legendlabels,lty=ltys,lwd=lwds)
abline(v=.5,lty=1,lw=2,col='red')
abline(v=4,lty=1,lw=2,col='red')
dev.off()

xlim <- 100
tiff(paste0(scenario.results.folder, scenario.name,'_Species_protection_',xlim,'km.tiff'), units="cm", width=16, height=12, res=200)
plot(MPAsizes/1000, pspeciesprotected[1,],type='l',lty=ltys[1],lwd=lwds,xlim = c(0,xlim),
     xlab = "MPA size (km)", ylab = "Species protected (%)")
axis(side = 1, at = seq(0,xlim,xlim/10))
for (p in 2:length(prot.levels)) {
  lines(MPAsizes/1000, pspeciesprotected[p,],type='l',lty=ltys[p],lwd=lwds)
}
legend('bottomright',legendlabels,lty=ltys,lwd=lwds)
abline(v=.5,lty=1,lw=2,col='red')
abline(v=4,lty=1,lw=2,col='red')
dev.off()


# Plot lifetime fishing mortality and fishing mortality offsets ####

xlim <- 10
ylim <- 100
ltys <- 1
lwds <- 1
fmortplot <- c(0.2,0.5)
fmortcols <- match(fmortplot,names(returndat[[1]]$mean_Fmortality))
lfmorts <- array(dim=c(nspecies,length(MPAsizes),length(fmortplot)))
foffs <- lfmorts
fsurvs <- lfmorts

# lifetime fishing mortality
for (plot in 1:length(fmortcols)){
  lfmorts[1,,plot] <- returndat[[s]]$mean_Fmortality[,fmortcols[plot]]
  tiff(paste0(scenario.results.folder, scenario.name,'_Lifetime_mortality_F',fmortplot[plot]*100,'_',xlim,'km.tiff'), units="cm", width=16, height=12, res=200)
  plot(MPAsizes/1000,lfmorts[1,,plot] * 100,
       type='l',lty=ltys[1],lwd=lwds,xlim = c(0,xlim),ylim = c(0,ylim),
       xlab = "MPA size (km)", ylab = "Lifetime fishing mortality (%)",
       main = paste0(fmortplot[plot] * 100, '% annual fishing mortality'))
  axis(side = 1, at = seq(0,xlim,xlim/10))
  if(length(returndat)>1){
    for (s in 2:length(returndat)) {
     lfmorts[s,,plot] <- returndat[[s]]$mean_Fmortality[,fmortcols[plot]]
     lines(MPAsizes/1000, lfmorts[s,,plot] * 100,lty=ltys,lwd=lwds)
    }
  }
  dev.off()
}

# fishing mortality offsets
ylim <- 1
for (plot in 1:length(fmortcols)){
  foffs[1,,plot] <- returndat[[s]]$mean_Foffset[,fmortcols[plot]]
  foffs[1,foffs[1,,plot]>1,plot] <- 1 # cap to mean risk of 1
  tiff(paste0(scenario.results.folder, scenario.name,'_Fmortality_risk_F',fmortplot[plot]*100,'_',xlim,'km.tiff'), units="cm", width=16, height=12, res=200)
  plot(MPAsizes/1000,foffs[1,,plot],
       type='l',lty=ltys[1],lwd=lwds,xlim = c(0,xlim),ylim = c(0,ylim),
       xlab = "MPA size (km)", ylab = "Fishing mortality risk (%)",
       main = paste0(fmortplot[plot]*100, '% annual fishing mortality'))
  axis(side = 1, at = seq(0,xlim,xlim/10))

  if(length(returndat)>1){
   for (s in 2:length(returndat)) {
    foffs[s,,plot] <- returndat[[s]]$mean_Foffset[,fmortcols[plot]]
    foffs[s,foffs[s,,plot]>1,plot] <- 1
    lines(MPAsizes/1000, foffs[s,,plot],lty=ltys,lwd=lwds)
   }
  }
  dev.off()
}


# Probability of survival
ylim <- 1
for (plot in 1:length(fmortcols)){
  fsurvs[1,,plot] <- returndat[[s]]$mean_survival_probability[,fmortcols[plot]]
  fsurvs[1,1,plot] <- 1-fmortplot[plot]
  tiff(paste0(scenario.results.folder, scenario.name,'_Psurvival_F',fmortplot[plot]*100,'_',xlim,'km.tiff'), units="cm", width=16, height=12, res=200)
  plot(MPAsizes/1000,fsurvs[1,,plot],
       type='l',lty=ltys[1],lwd=lwds,xlim = c(0,xlim),ylim = c(0,ylim),
       xlab = "MPA size (km)", ylab = "Annual probability of survival",
       main = paste0(fmortplot[plot]*100, '% annual fishing mortality'))
  axis(side = 1, at = seq(0,xlim,xlim/10))
  
  if(length(returndat)>1){
   for (s in 2:length(returndat)) {
    fsurvs[s,,plot] <- returndat[[s]]$mean_survival_probability[,fmortcols[plot]]
    fsurvs[s,1,plot] <- 1-fmortplot[plot]
    lines(MPAsizes/1000, fsurvs[s,,plot],lty=ltys,lwd=lwds)
   }
  }
  dev.off()
}
