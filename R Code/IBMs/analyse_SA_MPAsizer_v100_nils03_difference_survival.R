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
fmortplot <- c(0.05,0.2,0.5,1) # lifetime mortality plots
foffplot <- c(1) # annual mortality offset/risk/survival plots (all the same so need to pick one level) 
idealized.maxn <- 1
para.scenarios <- data.frame(age=c("Max","Max","Maturity"),mean.max.dist=c(1,0,0))

# plot specs
lwd <- 2 # line width
res <- 300 # resolution of plot (pixels)
xlabel.size <- "Protected area size (km)"
xlabel.area <- "Protected area size (km^2)"

# Load model and data
modelversion <- paste0("SA_MPAsizer",modelid)
modelname <- paste0('SA_MPASizer_v',modelid,'.R')
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

# allocate mortality matrix
lfmorts <- array(dim=c(nspecies,length(MPAsizes),length(fmortplot),dim(para.scenarios)[1]))


for (scenario in 1:dim(para.scenarios)[1]){ # loop over scenarios
  
  scenario.name <-case.study.region
  # determine scenario name
  age <- para.scenarios$age[scenario]
  mean.max.dist <- para.scenarios$mean.max.dist[scenario]
  if(idealized.maxn == 1){scenario.name <- paste0(case.study.region,'_idealized')}
  if(mean.max.dist == 1){scenario.name <- paste0(scenario.name,'_meanMaxDist')}
  scenario.name <- paste0(scenario.name,'_Res', resolution,'_',dispersal.period,'_AgeAt',age,'_rBRUV',rBRUVcatchment)

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
 
    fmortcols <- match(fmortplot,names(returndat[[1]]$mean_Fmortality))
    
    for (plotnum in 1:length(fmortcols)){
      lfmorts[1,,plotnum,scenario] <- 1-returndat[[1]]$mean_Fmortality[,fmortcols[plotnum]]
      
      for (s in 1:nspecies) {
        lfmorts[s,,plotnum,scenario] <- 1-returndat[[s]]$mean_Fmortality[,fmortcols[plotnum]]
      } # species loop
    } # mortality loop
    
    if(scenario==2){
      lfmorts.diff21 <- lfmorts[,,,2]-lfmorts[,,,1] # compare to mean scenario
    }
    if(scenario==3){
      lfmorts.diff32 <- lfmorts[,,,3]-lfmorts[,,,2] # compare to mean scenario
    }
    if(scenario==3){
      lfmorts.diff31 <- lfmorts[,,,3]-lfmorts[,,,1] # compare to mean scenario
    }
    
} # scenario loop
    


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Plot lifetime fishing mortality and fishing mortality offsets ####
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

library(RColorBrewer)

# specify colours and line types to use for different species/groups
ucols <- c("Blues","Greens","Reds"); 
ltys <- rep(1,sum(ngroups)); 
#cols <- NA
for(g in 1:length(ngroups)){
#ltys <- c(na.omit(ltys),1:ngroups[g])
cols <- c(na.omit(cols),brewer.pal(ngroups[g],ucols[g]))} # specify colours


lfmorts.diff <- list()
lfmorts.diff[[1]] <- lfmorts.diff21
lfmorts.diff[[2]] <- lfmorts.diff32

s.ylabel2 <- 'Difference in survival (%)'
#xlim <- 10
ylim.left <- 100
#ltys <- 1
#lwds <- 1
xlim <- 100
cex.labels <- 0.7

# lifetime fishing mortality
for (plotnum in 1:length(fmortcols)){
  
    tiff(paste0('Scenario_differences_Lifetime_survival_F',
              fmortplot[plotnum]*100,'_',xlim,'km.tiff'), units="cm", 
              width=14, height=16, res=res)
  
    par(mfrow = c(3,2))
  
    for(scenario in 1:dim(para.scenarios)[1]){ # loop over scenarios
      
      age <- para.scenarios$age[scenario]
      mean.dist <- para.scenarios$mean.dist[scenario]
      
      if(age == 'Max'){s.ylabel <- 'Lifetime survival (%)'}
      if(age == 'Maturity'){s.ylabel <- 'Survival to reproduction (%)'} # Survival to reproduction
    
      #xlim <- 10
      #ylim <- 100
      #ltys <- 1
      #lwds <- 1
      #xlims <- c(10,100)
        
        #for(xl in 1:length(xlims)){
          
          #xlim <- xlims[xl]
          
            
            #tiff(paste0(scenario.results.folder, scenario.name,'_Lifetime_survival_F',
            #            fmortplot[plotnum]*100,'_',xlim,'km.tiff'), units="cm", width=18, height=12, res=res)
            
            # plot left panel
            
            par(mar=c(4.1, 4.1, 0.1, 2.1)) # this is usually the default
            
            plot(MPAsizes/1000,lfmorts[1,,plotnum,scenario] * 100,
                 type='l',col=cols[species.order[1]],lty=ltys[species.order[1]],lwd=lwd,xlim = c(0,xlim),ylim = c(0,ylim.left),
                 xlab = NA, ylab = s.ylabel,cex=cex.labels) # xlab = xlabel.size, 
                 #main = paste0('Age at ',age,': ',fmortplot[plotnum] * 100, '% annual fishing mortality'))
            axis(side = 1, at = seq(0,xlim,xlim/10))
            if(length(returndat)>1){
              for (s in 1:nspecies) {
               lines(MPAsizes/1000, lfmorts[s,,plotnum,scenario] * 100,col=cols[s],lty=ltys[s],lwd=lwd)
              }
            }
    
            if(scenario==3){mtext(side=1, line=3,xlabel.size,cex=cex.labels)}
            
            
            # plot right panel 
            
            if (scenario==1){
              plot.new()
              par(xpd=TRUE)
              legend("topleft", inset=c(0,0), 
                     legend=str_replace_all(species.names.ordered,'_',' '), 
                   lty=ltys,col=cols,lwd=2)
              par(xpd=FALSE)
              } #, title="Group")
            #dev.off()
            
            
            if (scenario>1){
              
              slfmorts.diff <- lfmorts.diff[[scenario-1]] # define scenario
              xlim.loc <- max(which(MPAsizes/1000<=xlim))
              ylim.min <- floor(min(slfmorts.diff[,1:xlim.loc,plotnum],na.rm=TRUE)*10)*10
              ylim.max <- ceiling(max(slfmorts.diff[,1:xlim.loc,plotnum],na.rm=TRUE)*10)*10
              
              plot(MPAsizes/1000,slfmorts.diff[1,,plotnum] * 100,
                   type='l',col=cols[species.order[1]],lty=ltys[species.order[1]],lwd=lwd,
                   xlim = c(0,xlim),ylim = c(ylim.min,ylim.max),
                   xlab = NA,ylab = s.ylabel2) # xlab = xlabel.size, 
              #main = paste0('Age at ',age,': ',fmortplot[plotnum] * 100, '% annual fishing mortality'))
              
              axis(side = 1, at = seq(0,xlim,xlim/10))
              
              #if(scenario==3){axis(side=1,xlab = xlabel.size)
              if(length(returndat)>1){
                for (s in 1:nspecies) {
                  lines(MPAsizes/1000, slfmorts.diff[s,,plotnum] * 100,col=cols[s],lty=ltys[s],lwd=lwd)
                }
              }
               
            } # scenario condition
            
            if(scenario==3){mtext(side=1, line=3,xlabel.size,cex=cex.labels)}
            
          } # scenario loop
  
          dev.off()
    
#    } # scale loop

    
      
} # plot loop

    
