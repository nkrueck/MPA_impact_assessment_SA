##########################################################
# Shark MPA model plots of movement profiles by Nils, 1 Nov 2018 #
# Also plot of modelling environment

rm(list = ls()) # delete parameters in the workspace

#versionfolder <- "shark_mpa_model/Nils"
versionfolder <- ""
#install.packages('stringr')
library('stringr')

# load data
paste('Need to specifiy datadir - wherever you place the folder')
#datadir <- paste0('C:/Users/uqnkruec_local/Dropbox/Papers/Shark MPAs/', versionfolder ,'/')
datadir <- ""
#setwd(datadir)

#1: daily; 2: weekly; 3: monthly; 4: yearly
timescale <- 2 

#maxddata <- read.csv(paste0(datadir, "maxddata.csv"))

#Dispersal_Timescales <- readRDS("Dispersal_Timescales.RDS")
Dispersal_Timescales <- readRDS("2018-11-01_Dispersal_Timescales.RDS")

maxddata <- data.frame(Dispersal_Timescales[[timescale]]) # Load the monthly data
maxddata$species_name <- maxddata$species
maxddata <- maxddata[!is.na(maxddata$common_name),] # for some odd reason there are NAs in the dispersal file
maxddata$species <- str_extract(maxddata$species, '[^ ]+$')
rBRUVcatchment = 400
maxd.per.ind = 1
probs.per.ind = 1

if(timescale == 1) {timescalestring = "daily"}
if(timescale == 2) {timescalestring = "weekly"} 
if(timescale == 3) {timescalestring = "monthly"}
if(timescale == 4) {timescalestring = "annual"}
savefolder = paste0('results_', timescalestring)

movedat <- maxddata
speciesToSimulate <- c(5,2,3,1,4)

library(data.table)

# explore nmax data
lat.spec.names <- unique(maxddata$species)
comm.spec.names <- unique(maxddata$common_name)
#speciesToAnalyse = c(as.character(lat.spec.names[speciesToSimulate]),"") # has to include empty cells melanopterus","") #c("obesus","melanopterus","amblyrhynchos","perezi","cirratum") 
speciesToAnalyse = c(as.character(lat.spec.names)) # has to include empty cells melanopterus","") #c("obesus","melanopterus","amblyrhynchos","perezi","cirratum")  
nspecs <- length(speciesToSimulate) #updated when turning into function

# allocate output and working matrices
sprobslist <- vector("list")
probslist <- vector("list")
probsdflist <- vector("list")
distslist <- vector("list")
distsdflist <- vector("list")
sallprobs <- vector("list")
salldists <- vector("list")

##################
# START OF MODEL #
##################

## Only have to run this once at the beginning

# count for indexing purposes
count <- 0

#creates a series of list object which has the 5 species in the order  
# Whitetip, Blacktip, Grey, Caribbean, Nurse

# species loop
for (s in 1:nspecs){  # Removed species loop to run as part of a function

  count <- count + 1  
  # select data
  spec = speciesToSimulate[s] 
  specname = comm.spec.names[spec]
  speclatname = speciesToAnalyse[spec]
  
  # speciesToSimulate[s] # This has been standardised to a 1 as a result of turning the loop into a function which runs independently on each species
  #print(speclatname)
  #print(spec)
  
  if (maxd.per.ind == 1){
    ids = unique(maxddata$tag_id[maxddata$species == speclatname]) # get shark ids
    dists = matrix(nrow=length(ids),ncol=1)
    #print(speclatname) # latin name for species
    #print(spec) # number in vector 5 2 3 1 4
    #print(length(ids)) # how many individuals
    for (id in 1:length(ids)) { # Max distance moved per individual
      dists[id] = max(maxddata$max.dispersal_km[maxddata$tag_id == ids[id]]) * 1000
      }
  }else{
    dists = maxddata$max.dispersal_km[maxddata$species == speclatname] * 1000
  } # travel distances in m
  #hist(dists,breaks=length(dists)/4)
  #print(length(dists))
  #print(length(ids))
  
  
  # # RD This is the bit which numbers the sample size incorrectly
  if (maxd.per.ind == 1 & probs.per.ind == 1){
    # generate probability of occurrence along the distance spectrum
    for (id in 1:length(ids)) {
      alldists = maxddata$max.dispersal_km[maxddata$tag_id == ids[id]] * 1000
      # hist(alldists, breaks = 20, xlab = "Travel distance (m)", ylab = "Frequency")
      kernel = as.matrix(rep(1,1,max(alldists)))
      centreloc = ceiling(max(alldists)/2)
      for (dshark in 1:length(unique(alldists))) { # for each observed travel distance (RD WHY dshark???)
        addval = 1 # set initial number of observed vals to 1
        maxlocs = which(alldists == max(alldists)) # identify all distance locs
        startloc = round(centreloc - alldists[maxlocs[1]]/2) # define kernel start loc
        endloc = floor(startloc + alldists[maxlocs[1]]) # define kernel end loc
        if (length(maxlocs) > 1){addval = length(maxlocs)} # if observed more than once update value
        kernel[startloc:endloc] = kernel[startloc:endloc] + addval # add observations along kernel spectrum
        alldists = alldists[alldists != alldists[maxlocs[1]]] # delete observed distance from list
      }
      probsid = kernel/sum(kernel,na.rm=TRUE)
      probsid[which(is.na(probsid))] = 0

      probslist[[id]] = probsid # get cumulative probability spectrum by dividing by sum
      distslist[[id]] = seq(-centreloc+1,centreloc-1,1)
      
      # generate movement probability dataframes for plotting
      #probsdf <- data.frame(ID=as.character(ids[id]),
      #                      DISTANCE=1:length(probslist[[id]]),
      #                      P=probslist[[id]])
      #probsdflist[[id]] <- probsdf

    } # end of individual loop
    salldists[[s]] <- distslist
    sallprobs[[s]] <- probslist
    
    ### IMPORTANT FIND
    distslist <- NULL #RD Reset these values
    probslist <- NULL #RD Reset these values
  } # end of species loop

  #print(lengths(sallprobs))
  
  #probsdflistdf <- do.call(rbind,probsdflist)
  # probsdflistdf <- data.table(probsdflistdf)
  #return(probsdflistdf)
  #return(list(probslist,distslist))
}


##########################
# Plot movement profiles

# Emulate ggplot colour hues seen in the other figures
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols = gg_color_hue(5)

xlabel = "Distance from origin (km)"
ylabel = "Relative likelihood of occurrence"
linewidth1 = 0.7 # width of lines for individuals
linewidth2 = 1.5 # width of line for average
distlim = 5 # cut off plots at certain distance (km) or not (distlim = 0) 
distticks = distlim/5 # specify ticks in plot
standardize = 1 # standardize probabilities based on maximum

meanlist <- list(0) # list object where mean values are saved

# Plot all the species on a single figure
speciesToPlot = c(1:5)
nspecsplotting = length(speciesToPlot)

# Names for the top left of plots
figi <- c("A",
          "B",
          "C",
          "D",
          "E")
resultsdir <- paste0('Ross/results/week/rBRUVcatchment_400/') # Better if working with projects


## No need to run this as it takes ages and only rewrites the figure

# png(filename = paste0(resultsdir,"Figure 1_Dispersal kernels_species.png"),
#     width = 220, height = 150, units = "mm",
#     res=600)

pdf(file = paste0(resultsdir,"Figure 1.pdf"),
    width = 8.74, height = 5.91)
    #width = 220, height = 150, units = "mm",
    #res=600)


par(mfrow=c(2,3),las=1,
    oma = c(3, 3, 0, 0), # make room (i.e. the 4's) for the overall x and y axis titles
    mar = c(2, 2, 2, 2), # make the plots be closer together
    xpd = TRUE)  # Allow legend plotting in margins

for (s in 1:5){
  
  #s=2
  #count <- count + 1
  
  spec <- speciesToPlot[s]
  probs <- sallprobs[[spec]]
  dists <- salldists[[spec]]
  
  maxprob = 1 
  if (standardize == 0) {maxprob <- max(sapply(probs,max,na.rm=T),na.rm=T)}
  if (distlim == 0){distlim = max(sapply(dists,max))}
  
  id_dists = seq(-distlim*1000,distlim*1000,1)+distlim*1000+1 # vector of cell positions
  centreloc = median(id_dists) # centre location of vector
  plot_dists = seq(-distlim,distlim,1/1000) # distance along vector
  
  sumprobs = rep(NA,length(id_dists)) # initial probaility vector of zeros 
  
  plot(plot_dists,sumprobs,type="l",
       ylim=c(0,maxprob),
       xlim=c(-distlim,distlim),
       xaxt = "n", 
       yaxt = "n", 
       lwd=linewidth1,
       col="grey",
       xlab = "", 
       ylab = "",
       bty='l',
       cex.main=1, cex.sub=1)
  axis(1, seq(-distlim,distlim,distticks),cex.axis=1.5)
  axis(2, seq(0,maxprob,maxprob/5),cex.axis=1.5)

  for (id in 1:length(probs)) {
    # get data
    probsdat = probs[[id]]
    if (standardize == 1){probsdat <- probsdat/max(probsdat)} # standardize if specified
    # get indices for plotting
    id_probs = seq(centreloc-(length(probs[[id]])-1)/2,centreloc+(length(probs[[id]])-1)/2,1) # get vecor of cell IDs
    if (length(id_probs) > length(id_dists)){ # if probs extent spefied distance range to plot
      idx_sprobs <- id_dists; idx_probs <- id_dists + (length(id_probs) - length(id_dists))/2
    } else if (length(id_probs) <= length(id_dists)){
      idx_sprobs <- id_probs; idx_probs <- seq(1,length(probsdat),1)}
    plotprobs = rep(0,length(id_dists))
    plotprobs[idx_sprobs] = probsdat[idx_probs]
    lines(plot_dists,plotprobs,lwd=linewidth1,col="grey")
    sumprobs = apply(rbind(sumprobs,plotprobs),2,sum,na.rm=T)
  }
  meanprobs = sumprobs/length(probs)
  lines(plot_dists,meanprobs,lwd=linewidth2,col= cols[s])

  meanlist[[s]] <- data.frame(Species=s,plot_dists,meanprobs)

  text(-4.5, 1.0, figi[s], cex=2,font=2)
  text(4, 1.0, paste0("n = ",length(probs)), cex=1.5)
  
  }

mtext(xlabel, side = 1, outer = TRUE, line = 1, cex=1)
par(las=0)
mtext(ylabel, side = 2, outer = TRUE, line = 1, cex=1)

plot.new() 

legend("bottomright", title="",
       c("Whitetip reef shark","Blacktip reef shark","Grey reef shark","Caribbean reef shark","Nurse shark"),
       lty= 1,
       lwd=3,
       bty="n",
       col = cols, 
       inset=c(0.1,0.3),
       cex=1.5)

dev.off()





###############

## how many fixes do we have per animal

maxddata_days <- data.frame(Dispersal_Timescales[[1]])

library(dplyr)

days_sum <- maxddata_days %>%
  select(date,tag_id,common_name) %>%
  na.omit() %>%
  group_by(tag_id,common_name) %>%
  summarize(number = n()) %>%
  arrange(number)

days_sum %>%
  group_by(common_name) %>%
  summarize(number = n()) %>%
  arrange(number)

# How many fixes with >n days
days_sum7 <- days_sum %>%
  filter(number > 7)

days_sum7 %>%
  group_by(common_name) %>%
  summarize(number = n()) %>%
  arrange(number)

###############
## Generate an animation of our plots




#count <- count + 1
s=3
par(mfrow=c(1,1),las=1,
    oma = c(3, 3, 0, 0), # make room (i.e. the 4's) for the overall x and y axis titles
    mar = c(2, 2, 2, 2), # make the plots be closer together
    xpd = TRUE)  # Allow legend plotting in margins
  
  spec <- speciesToPlot[s]
  probs <- sallprobs[[spec]]
  dists <- salldists[[spec]]
  
  maxprob = 1 
  if (standardize == 0) {maxprob <- max(sapply(probs,max,na.rm=T),na.rm=T)}
  if (distlim == 0){distlim = max(sapply(dists,max))}
  
  id_dists = seq(-distlim*1000,distlim*1000,1)+distlim*1000+1 # vector of cell positions
  centreloc = median(id_dists) # centre location of vector
  plot_dists = seq(-distlim,distlim,1/1000) # distance along vector
  
  saveHTML({
 
  for (pl in 1:length(probs)){
    #for (pl in 1:3){
    sumprobs = rep(NA,length(id_dists)) # initial probaility vector of zeros 

    
     # png(filename = paste0("Ross/animation/Figure 1_",pl,".png"),
     #     width = 220, height = 150, units = "mm",
     #     res=600)
     # 
     
    plot(plot_dists,sumprobs,
         type="l",
         ylim=c(0,maxprob),
         xlim=c(-distlim,distlim),
         xaxt = "n", 
         yaxt = "n", 
         lwd=linewidth1,
         col="white",
         xlab = "", 
         ylab = "",
         bty='l',
         cex.main=1, cex.sub=1)
    axis(1, seq(-distlim,distlim,distticks),cex.axis=1.5)
    axis(2, seq(0,maxprob,maxprob/5),cex.axis=1.5)
    
    for (id in 1:pl) {
      # get data
      probsdat = probs[[id]]
      if (standardize == 1){probsdat <- probsdat/max(probsdat)} # standardize if specified
      # get indices for plotting
      id_probs = seq(centreloc-(length(probs[[id]])-1)/2,
                     centreloc+(length(probs[[id]])-1)/2,1) # get vecor of cell IDs
      if (length(id_probs) > length(id_dists)){ # if probs extent spefied distance range to plot
        idx_sprobs <- id_dists; idx_probs <- id_dists + (length(id_probs) - length(id_dists))/2
      } else if (length(id_probs) <= length(id_dists)){
        idx_sprobs <- id_probs; idx_probs <- seq(1,length(probsdat),1)}
      plotprobs = rep(0,length(id_dists))
      plotprobs[idx_sprobs] = probsdat[idx_probs]
      
      #print(max(plotprobs))
      
      lines(plot_dists,plotprobs,lwd=linewidth1,col="grey")
      sumprobs = apply(rbind(sumprobs,plotprobs),2,sum,na.rm=T)
    }
    meanprobs = sumprobs/max(sumprobs)
    lines(plot_dists,meanprobs,lwd=linewidth2,col= cols[s])
    
    text(4, 1.0, paste0("n = ",pl), cex=1.5)
  }
  
    
    #dev.off()
 # }



###############

# Now the species comparison plot (i.e. all on the same panel)
library(ggplot2)
meanlinesdf <- do.call(rbind,meanlist)

meanlinesdf$Species <- as.factor(meanlinesdf$Species)
levels(meanlinesdf$Species) <- c("Whitetip reef shark",
                                 "Blacktip reef shark",
                                 "Grey reef shark",
                                 "Caribbean reef shark",
                                 "Nurse shark")


# Compare- proportion protected with increasing MPA size
## plot with species as each line 
ggplot(meanlinesdf, 
       aes(plot_dists, y = meanprobs)) + 
  geom_line(aes(colour = Species)) + 
  labs(linetype = "F")+
  theme_classic(base_size = 13) + 
  theme(#legend.position = c(0.5, 0.2),
        legend.position = c(0.2, 0.8),
        #legend.title=element_text(size=10),
        legend.title = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        strip.background = element_blank(),
        #strip.placement = "inside",
        strip.text=element_text(vjust=-1),
        panel.background = element_blank())+ 
  ylim(0,1)+
  xlab(xlabel) + 
  ylab(ylabel)+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)

###############

# Supp figure showing model structure

# Which species to visualise
s= 3

spec <- speciesToPlot[s] # which species to plot. #3 is Grey reef shark
probs <- sallprobs[[spec]] # extract probabilities
dists <- salldists[[spec]] # extract distances

maxprob = 1 
if (standardize == 0) {maxprob <- max(sapply(probs,max,na.rm=T),na.rm=T)} # standardise probabilities so max = 1
if (distlim == 0){distlim = max(sapply(dists,max))}

id_dists = seq(-distlim*1000,distlim*1000,1)+distlim*1000+1 # vector of cell positions
centreloc = median(id_dists) # centre location of vector
plot_dists = seq(-distlim,distlim,1/1000) # distance along vector

sumprobs = rep(NA,length(id_dists)) # initial probaility vector of zeros 

Maxdistshark <- max(do.call(rbind,lapply(dists,max)))/1000 # Max distance covered by a shark
distlim= (Maxdistshark+50)/2 # max area where we are working (i.e. the modelling area) that takes into consideration the MPA size and the max shark movement
#distlim=55 # Max and min limits (55 = 110 km modelling env)
distlim=30 # Max and min limits (30 = 60 km modelling env)
distticks=10 # plotted tick on x axis
maxprob = 20 # number of animals to be plotted
yid <- c(1:maxprob)

png(filename = paste0(resultsdir,"Supp Figure 4_model environment zoom_3_30.png"),
    width = 700, height = 480, units = "px")


pdf(file = paste0(resultsdir,"Figure S3.pdf"),
    width = 8.74, height = 5.91)
#width = 220, height = 150, units = "mm",
#res=600)


# Plot empty modelling environment
par(mfrow=c(1,1),
    bty="o",
    mar = c(4,1,4,10),
    xpd=TRUE
)
plot(plot_dists,sumprobs,type="n",
     ylim=c(0,maxprob),
     xlim=c(-distlim,distlim),
     xaxt = "n", 
     yaxt = "n", 
     xlab = "Coastline (km)", 
     ylab = "",
     cex.main=1.5,
     cex.lab=1.5,
     cex.sub=1.5)
axis(side = 1, 
     at = seq(-distlim,distlim,distticks),
     labels = seq(0,distlim*2,distticks),
     cex.axis=1.5)


# MPA polygon
#xpol2 <- c(-1,1,1,-1,-1) # 2 km MPA
xpol10 <- c(-5,5,5,-5,-5) # 10 km MPA
xpol20 <- c(-10,10,10,-10,-10) # 20 km MPA
xpol30 <- c(-15,15,15,-15,-15) # 30 km MPA
#xpol100 <- c(-50,50,50,-50,-50) # 100 km MPA
ypol <- c(-.4,-.4,maxprob+.4,maxprob+.4,-.4)
#plot the 20 km MPA boundary
polygon(xpol30, ypol,
        col="lightgrey",border="darkgrey",
        lty=1)

# polygon(xpol100, ypol,
#         border="lightgrey",lty=2)

# Add text at the header of the plot
text(-20, maxprob+2.5, "Fished", cex=1.5)
text(-20, maxprob+1.5, "F > 0", cex=1, col="darkgrey")
text(0, maxprob+2.5, "MPA", cex=1.5)
text(0, maxprob+1.5, "F = 0", cex=1, col="darkgrey")
text(20, maxprob+2.5, "Fished", cex=1.5)
text(20, maxprob+1.5, "F > 0", cex=1, col="darkgrey")

# Positioning of BRUVS
bruvspace <- 2  # suggests bruvs are 2 km apart
xpoi <- seq(-distlim,distlim,by=bruvspace) 

all0s <- data.frame(xid1=xpoi,Freq1=rep(0,length(xpoi)))

# Pick how many sharks were detected at each BRUV
set.seed(3)
xid1 <- base::sample(xpoi,size=(maxprob-1),replace = TRUE) # recordings at BRUVS
xid1_df <- data.frame(table(xid1))
xid1_df[,1] <- as.numeric(as.character(xid1_df[,1]))
xid1_df[,2] <- as.character(xid1_df[,2])

library(dplyr)
xid1_df1 <- left_join(all0s,xid1_df)
xid1_df1$Freq[is.na(xid1_df1$Freq)] <- 0

# Plot BRUV locations along foot of image
#points(xpoi,rep(0,length(xpoi)),pch="0",cex=1,col="grey")

# plot counts
points(xid1_df1[,1],rep(0,nrow(xid1_df1)),pch=xid1_df1[,3],cex=1,col="black")

# Randomise where on the BRUVS spectrum a shark sits
set.seed(5)
xid <- xid1 + base::sample(-.2:.2) #

# Hard wire % overlap with MPA
perover <- c(0,0,100,45,100,100,0,0,100,100,100,100,100,100,0,0,0,10,0)
mean(perover)
sd(perover)/sqrt(length(perover))

#perover <- 0

##


# Plot the shark movement profiles
for (id in 1:(maxprob-1)) {
  
  probsdat = probs[[id]]
  if (standardize == 1){probsdat <- probsdat/max(probsdat)} # standardize if specified
  
  # get indices for plotting
  id_probs = seq(centreloc-(length(probs[[id]])-1)/2,centreloc+(length(probs[[id]])-1)/2,1) # get vecor of cell IDs 
  if (length(id_probs) > length(id_dists)){ # if probs extent spefied distance range to plot
    idx_sprobs <- id_dists; idx_probs <- id_dists + (length(id_probs) - length(id_dists))/2 
  } else if (length(id_probs) <= length(id_dists)){
    idx_sprobs <- id_probs; idx_probs <- seq(1,length(probsdat),1)}
  plotprobs = rep(0,length(id_dists))
  plotprobs[idx_sprobs] = probsdat[idx_probs]
  
  # Makes the lines look pretty
  first0 <-  head(which(plotprobs!=0),1) # first value which isn't zero    
  last0 <-  tail(which(plotprobs!=0),1) # last value which isn't zero
  plotprobs[plotprobs == 0] <- NA #replace all zeros wih NAs
  plotprobs[(first0-500):first0] <- 0 # Make the lines fix to the zero probability
  plotprobs[last0:(last0+500)] <- 0
  
  #Add the probability lines on the figure
  lines(xid[id]+ plot_dists,
        yid[id]+ plotprobs,
        lwd=2,
        col=cols[3])
  
   # Text containing % overlap with MPA for each individual (FOR 50 km mobelling env)
   text(37, maxprob+0.5, "Protection", cex=1.0,col="black")
   #text(36, id+0.5, paste0(perover[id],"%"), cex=1.0,col="black")
  
  sumprobs = apply(rbind(sumprobs,plotprobs),2,sum,na.rm=T)
}

# Add legend to top of plot
# legend("topright", inset=c(-0.28,-0.1),
legend("topright",
       legend=c("MaxN at BRUV",
                "Shark movement profile"),
       col = c("#808080",cols[3]),
       bty="o",
       text.col = "black",
       lty = c(-1, 1),
       lwd = c(-1, 2),
       pch = c("0", NA),
       merge = TRUE, bg = "white")

dev.off()

###############

