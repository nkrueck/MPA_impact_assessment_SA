##########################################################
# Exploratory look at data by Ross Dwyer, 11 Sep 2018 #
##########################################################

# in contrast to v1 this version can calculate exposure also based on movement probabilities
# which is initiated by setting probs.per.ind to 1

rm(list = ls()) # delete parameters in the workspace

versionfolder <- "shark_mpa_model_v301"

# load data
paste('Need to specifiy datadir - wherever you place the folder')
#datadir <-  paste0('C:/Users/uqnkruec_local/Dropbox/Papers/Shark MPAs/', version,'/')
#datadir <- #paste0('C:/Users/weizenkeim/Dropbox/Papers/Shark MPAs/', version,'/')
#datadir <- paste0('E:/OneDrive/Australia/SharkRay MPAs/Agent based models/', versionfolder,'/')
datadir <- ""

maxndata <- read.csv(paste0(datadir, "maxndata.csv")) # RD added perezi and set location_code to 4
maxddata <- read.csv(paste0(datadir, "maxddata.csv"))

## Explore the data ##

library(ggplot2)
library(scales)
library(dplyr)
library(broom)

sharkspecies <- c("Whitetip Reef Shark" , "Blacktip Reef Shark", "Grey Reef Shark", "Caribbean Reef Shark")

##  Generate mean dispersal distances table for all 5 species combined
source(paste0(datadir,'summarySE.R'))

# MEan max per id
sumddata_species_all <- maxddata %>%
  summarySE_disp(measurevar="max.dispersal_km", groupvars=c("species","tag_id"))
# MEan max per species
sumddata_species_all %>%
summarySE_disp(measurevar="max", groupvars=c("species"))

                                                                     



## Dispersal distance plots

# Generate histogram of frequency of dispersal distances 
ggplot(data=maxddata, aes(max.dispersal_km, fill=common_name)) + 
  stat_density(aes(y=stat(count))) +
  scale_x_continuous(breaks=c(0,1,2,3,4,5,10,30,100,300,1000), trans="log1p", expand=c(0,0)) +
  theme_bw()+
  facet_wrap( ~ common_name, ncol=2)

# Generate histogram of proportions of dispersal distances 
ggplot(data=maxddata, aes(max.dispersal_km, fill=common_name)) + 
  stat_density(aes(y=..density..)) +
  scale_x_continuous(breaks=c(0,1,2,3,4,5,10,30,100,300,1000), trans="log1p", expand=c(0,0)) +
  theme_bw()+
  facet_wrap( ~ common_name, ncol=2)

source(paste0(datadir,'summarySE.R'))


#####################################

'%!in%' <- function(x,y)!('%in%'(x,y))

### Now for the BRUV observation frequencies

# Function to generate species specific bruv data with correct 'absences
fextract_sp <- function(data,sharksp){
  
  data <- data %>% select(-starts_with("location_code"))
  
  #sharksp = "Whitetip reef shark"
  
  ## BRUV frequency of shark occurrence of species named in our sharksp vector 
  maxndata_spec_shark <- data[data$common_name %in% sharksp, ]  
  
  shark_sites <- unique(as.character(maxndata_spec_shark$site_name))   # Which sites have logged at least one of the species named in our sharksp vector?
  maxndata_sites <- data[data$site_name %in% shark_sites, ] # Filter those bruv sites which have logged at least one of the species named in our sharksp vector
  
  # this step consolidates multiple species into a single category (stops duplicate lines for multiple species)
  if(length(sharksp)>1){ 
    # Sum max n values for multiple species 
    unique_maxn <- with(maxndata_spec_shark,aggregate(maxn, list(set_id = set_id), sum))
    names(unique_maxn)[2] <- "maxn"  
    
    # Generate a unique row in BRUV dataset accordng to set code
    subsetteddf <- maxndata_spec_shark %>%
      select("region_name","location_name","site_name","reef_name","trip_year",
             "set_code","set_lat","set_long","depth","reef_id","set_id") %>%
      distinct() 
    subsetteddf$id <- 1:nrow(subsetteddf) 
    subsetteddf$common_name <- "multiple"
    subsetteddf$genus <- "multiple"
    subsetteddf$species <- "multiple"  
    subsetteddf$animal_id <- 9998
    subsetteddf$location_code <- 2
    subsetteddf$location_code_high <- 2
    
    # Join the 2 dataframes and rearrange to correct order and re-save as maxndata_spec_shark
    subsetteddf_maxn <- left_join(subsetteddf,unique_maxn)
    maxndata_spec_shark <- subsetteddf_maxn %>%
      select(names(maxndata_spec_shark)) 
  }
  
  # Which drops were in a site with one of our shark species but did not see ANYTHING (fish, turtles or sharks)? 
  maxndata_spec_noanimals <- maxndata_sites[maxndata_sites$animal_id %in% 9999, ]
  
  # Which drops were in a site with one of our named shark species, but didnt see any of these species on a survey? 
  sharksp_id <-  unique(maxndata_spec_shark$animal_id) # ids for each of our 5 species
  maxndata_noshark <- maxndata_sites[maxndata_sites$animal_id %!in% c(sharksp_id,9999), ] # removes these ids and the one with no animals
  noshark_setids <- sort(unique(as.character(maxndata_noshark$set_id))) # what are the set ids where species other than our 5 shark speces were spotted  
  nosharksp_setids <- noshark_setids[noshark_setids %!in% maxndata_spec_shark$set_id] # Extract set_ids not present in the subsetted shark data
  maxndata_spec_nosharksdups <- maxndata_sites[maxndata_sites$set_id %in% nosharksp_setids, ] # filter dataset for these set codes
  maxndata_spec_nosharks <- maxndata_spec_nosharksdups[!duplicated(maxndata_spec_nosharksdups$set_id),] # remove the duplicated set ids where more than one species was spotted
  maxndata_spec_nosharks$animal_id <- 9999 # reset animal_id to 9999 for sites with fish but none of our 5 species
  maxndata_spec_nosharks$genus <- ""  # reset genus to blank for sites with fish but none of our 5 species
  maxndata_spec_nosharks$species <- "" # reset species to blank for sites with fish but none of our 5 species
  maxndata_spec_nosharks$maxn <- 0 # reset to 0 for sites with fish but none of our 5 species
  maxndata_spec_nosharks$common_name <- ""
  
  #Assign the correct zeros to the shark bruv data dataframe (bind our shark species, other animals and no animals dataframes)
  maxndata_spec <- rbind(maxndata_spec_shark,maxndata_spec_noanimals,maxndata_spec_nosharks)
  return(maxndata_spec)
}

sharksp <- c("Whitetip reef shark","Blacktip reef shark","Grey reef shark","Caribbean reef shark","Nurse shark")

# Run function to generate species specific bruv data with correct 'absences'
maxndata_spec_all <- fextract_sp(maxndata,sharksp)
maxndata_spec1 <- fextract_sp(maxndata,sharksp[1]) # Whitetip reef shark
maxndata_spec2 <- fextract_sp(maxndata,sharksp[2]) # Blacktip reef shark
maxndata_spec3 <- fextract_sp(maxndata,sharksp[3]) # Grey reef shark
maxndata_spec4 <- fextract_sp(maxndata,sharksp[4]) # Caribbean reef shark
maxndata_spec5 <- fextract_sp(maxndata,sharksp[5]) # Nurse shark
# Grey reef shark - remove 2 sites with very large GRS
maxndata_spec3.1 <- maxndata_spec3[maxndata_spec3$site_name %!in% c("Jarvis Island", "Beveridge Reef"),]
  
library(nlme)
# Run stats to compare atlantic vs pacific species
m1 <- glm(maxn ~ region_name,data= maxndata_spec_all)
anova(m1,test="F")  

library(lme4)
m0 <- glmer(maxn ~ 1 +
             (1|site_name),family=poisson(link = "log"),
           REML=FALSE,
           data= maxndata_spec_all) 

m1 <- lmer(maxn ~ region_name +
             (1|site_name),family=poisson(link = "log"),
           REML=FALSE,
           data= maxndata_spec_all) 

anova(m1,m0)

# Generate summarySE table for all 5 species combined
sumndata_site_all <- maxndata_spec_all %>%
  summarySE(measurevar="maxn", groupvars=c("site_name","region_name")) %>%
  select('site_name','region_name','Ndrops','maxn','se','ci') %>%
  rename(c('maxn' = 'meanmaxn')) 

# Generate summarySE table for whitetip reef sharks
sumndata_site_WRS <- maxndata_spec1 %>%
  summarySE(measurevar="maxn", groupvars=c("site_name","region_name")) %>%
  select('site_name','region_name','Ndrops','maxn','se','ci') %>%
  rename(c('maxn' = 'meanmaxn')) 

# Generate summarySE tables for blacktip reef sharks
sumndata_site_BRS <- maxndata_spec2 %>%
  summarySE(measurevar="maxn", groupvars=c("site_name","region_name")) %>%
  select('site_name','region_name','Ndrops','maxn','se','ci') %>%
  rename(c('maxn' = 'meanmaxn'))

# Generate summarySE tables for grey reef sharks
sumndata_site_GRS <- maxndata_spec3.1 %>%
  summarySE(measurevar="maxn", groupvars=c("site_name","region_name")) %>%
  select('site_name','region_name','Ndrops','maxn','se','ci') %>%
  rename(c('maxn' = 'meanmaxn'))

# Generate summarySE tables for caribbean reef sharks
sumndata_site_CRS <- maxndata_spec4 %>%
  summarySE(measurevar="maxn", groupvars=c("site_name","region_name")) %>%
  select('site_name','region_name','Ndrops','maxn','se','ci') %>%
  rename(c('maxn' = 'meanmaxn'))

# Generate summarySE tables for nurse sharks
sumndata_site_NS <- maxndata_spec5 %>%
  summarySE(measurevar="maxn", groupvars=c("site_name","region_name")) %>%
  select('site_name','region_name','Ndrops','maxn','se','ci') %>%
  rename(c('maxn' = 'meanmaxn'))

# Now draw the plot of mean max n coloured by Region
plotmaxn_region <- function(sumndata_site,plottitle){
  # set the levels in order we want
  sumndata_site$site_name <- factor(sumndata_site$site_name, 
                                    levels = sumndata_site$site_name[order(sumndata_site$meanmaxn,
                                                                           decreasing = TRUE)])
  newtitle <- paste0(plottitle,"\n")
  
  ## plot
  ggplot(sumndata_site, 
         aes(x = site_name, y = meanmaxn, fill = region_name, alpha = 0.5,legend=FALSE)) + 
    geom_errorbar(width=.1, aes(ymin=meanmaxn-se, ymax=meanmaxn+se), colour="black")+
    labs(title = newtitle, x = "BRUV site", y = "maxn", fill = "Region\n") +
    theme_bw() + 
    geom_bar(stat = "identity")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = c(0.8, 0.8))+
    scale_alpha(guide = 'none')
}

plotmaxn_region(sumndata_site_all,plottitle="All species")
plotmaxn_region(sumndata_site_WRS,plottitle="Whitetip reef shark") # 
plotmaxn_region(sumndata_site_BRS,plottitle="Blacktip reef shark") # Rangiroa - Maupiti
plotmaxn_region(sumndata_site_GRS,plottitle="Grey reef shark") # Jarvis - Beveridge
plotmaxn_region(sumndata_site_CRS,plottitle="Caribbean reef shark") #Exumas to Fernando de Noronha
plotmaxn_region(sumndata_site_NS,plottitle="Nurse shark") #Glovers to barbuda

####

## Now generate clusters for each species according to site max in
##https://uc-r.github.io/kmeans_clustering
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

siteclusters <- function(data,icenters){
  #data = sumndata_site_WRS
  df <- data[,c(1,4)]
  row.names(df) <- df[,1]
  df[,1] <- 0
  df <- scale(df)
  df[,1] <- 0
  distance <- get_dist(df)
  #fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
  k2 <- kmeans(df, centers = icenters, nstart = 25)
  #str(k2)
  return(k2$cluster)
}
## 2 clusters  - high - low
WRSclust.2 <- data.frame(siteclusters(sumndata_site_WRS,2))
BRSclust.2 <- data.frame(siteclusters(sumndata_site_BRS,2))
GRSclust.2 <- data.frame(siteclusters(sumndata_site_GRS,2))
CRSclust.2 <- data.frame(siteclusters(sumndata_site_CRS,2))
NSclust.2 <- data.frame(siteclusters(sumndata_site_NS,2))
## 3 clusters  - high - medium -  low
WRSclust.3 <- data.frame(siteclusters(sumndata_site_WRS,3))
BRSclust.3 <- data.frame(siteclusters(sumndata_site_BRS,3))
GRSclust.3 <- data.frame(siteclusters(sumndata_site_GRS,3))
CRSclust.3 <- data.frame(siteclusters(sumndata_site_CRS,3))
NSclust.3 <- data.frame(siteclusters(sumndata_site_NS,3))

### Now combine into lookup table
library(tibble)
library(dplyr)
allclust.2 <- full_join(rownames_to_column(WRSclust.2), rownames_to_column(BRSclust.2), by = ("rowname" = "rowname"))
allclust.2 <- full_join(allclust.2, rownames_to_column(GRSclust.2), by = ("rowname" = "rowname"))
allclust.2 <- full_join(allclust.2, rownames_to_column(CRSclust.2), by = ("rowname" = "rowname"))
allclust.2 <- full_join(allclust.2, rownames_to_column(NSclust.2), by = ("rowname" = "rowname"))
names(allclust.2) <- c("site_name","Whitetip.reef.shark","Blacktip.reef.shark",
                      "Grey.reef.shark","Caribbean.reef.shark","Nurse.shark")

allclust.3 <- full_join(rownames_to_column(WRSclust.3), rownames_to_column(BRSclust.3), by = ("rowname" = "rowname"))
allclust.3 <- full_join(allclust.3, rownames_to_column(GRSclust.3), by = ("rowname" = "rowname"))
allclust.3 <- full_join(allclust.3, rownames_to_column(CRSclust.3), by = ("rowname" = "rowname"))
allclust.3 <- full_join(allclust.3, rownames_to_column(NSclust.3), by = ("rowname" = "rowname"))
names(allclust.3) <- c("site_name","Whitetip.reef.shark","Blacktip.reef.shark",
                      "Grey.reef.shark","Caribbean.reef.shark","Nurse.shark")

write.csv(allclust.2,"kmeans_clusters_lowhigh.csv")
write.csv(allclust.3,"kmeans_clusters_lowmediumhigh.csv")
####

# Create pivot table for each species containing site groupings with means and ses
common_name <- "Whitetip reef shark"
WRSclust.23 <- full_join(rownames_to_column(WRSclust.2), rownames_to_column(WRSclust.3), by = ("rowname" = "rowname"))
names(WRSclust.23)<- c('site_name','location_code','location_code_high')
WRSclust.count23 <- full_join(sumndata_site_WRS,WRSclust.23)
WRSclust.count23 <- cbind(common_name,WRSclust.count23)
WRSclust.count23$region_name <- ifelse(WRSclust.count23$region_name=="Western Atlantic","Western Atlantic","Western Pacific")
WRSclust.count23 <- WRSclust.count23 %>%
  group_by(region_name) %>%
  arrange(site_name) %>%
  mutate(location_code = factor(location_code),
         location_code_high = factor(location_code_high))
names(WRSclust.count23)<- c('common_name','site_name','region_name','Ndrops','Mean','SE','CI','location_code','location_code_high')
medorder <- order(with(WRSclust.count23,tapply(Mean,location_code,mean))) # Ensures level ordering is correct
highorder <- order(with(WRSclust.count23,tapply(Mean,location_code_high,mean))) # Ensures level ordering is correct
levels(WRSclust.count23$location_code) <- list("1"=medorder[1], "2"=medorder[2])
levels(WRSclust.count23$location_code_high) <- list("1"=highorder[1], "2"=highorder[2], "3"=highorder[3])
with(WRSclust.count23,tapply(Mean,location_code,mean))
with(WRSclust.count23,tapply(Mean,location_code_high,mean))



common_name <- "Blacktip reef shark"
BRSclust.23 <- full_join(rownames_to_column(BRSclust.2), rownames_to_column(BRSclust.3), by = ("rowname" = "rowname"))
names(BRSclust.23)<- c('site_name','location_code','location_code_high')
BRSclust.count23 <- full_join(sumndata_site_BRS,BRSclust.23)
BRSclust.count23 <- cbind(common_name,BRSclust.count23)
BRSclust.count23$region_name <- ifelse(BRSclust.count23$region_name=="Western Atlantic","Western Atlantic","Western Pacific")
BRSclust.count23 <- BRSclust.count23 %>%
  group_by(region_name) %>%
  arrange(site_name) %>%
  mutate(location_code = factor(location_code),
         location_code_high = factor(location_code_high))
names(BRSclust.count23)<- c('common_name','site_name','region_name','Ndrops','Mean','SE','CI','location_code','location_code_high')
medorder <- order(with(BRSclust.count23,tapply(Mean,location_code,mean))) # Ensures level ordering is correct
highorder <- order(with(BRSclust.count23,tapply(Mean,location_code_high,mean))) # Ensures level ordering is correct
levels(BRSclust.count23$location_code) <- list("1"=medorder[1], "2"=medorder[2])
levels(BRSclust.count23$location_code_high) <- list("1"=highorder[1], "2"=highorder[2], "3"=highorder[3])

with(BRSclust.count23,tapply(Mean,location_code,mean))
with(BRSclust.count23,tapply(Mean,location_code_high,mean))

common_name <- "Grey reef shark"
GRSclust.23 <- full_join(rownames_to_column(GRSclust.2), rownames_to_column(GRSclust.3), by = ("rowname" = "rowname"))
names(GRSclust.23)<- c('site_name','location_code','location_code_high')
GRSclust.count23 <- full_join(sumndata_site_GRS,GRSclust.23)
GRSclust.count23 <- cbind(common_name,GRSclust.count23)
GRSclust.count23$region_name <- ifelse(GRSclust.count23$region_name=="Western Atlantic","Western Atlantic","Western Pacific")
GRSclust.count23 <- GRSclust.count23 %>%
  group_by(region_name) %>%
  arrange(site_name) %>%
  mutate(location_code = factor(location_code),
         location_code_high = factor(location_code_high))
names(GRSclust.count23)<- c('common_name','site_name','region_name','Ndrops','Mean','SE','CI','location_code','location_code_high')
medorder <- order(with(GRSclust.count23,tapply(Mean,location_code,mean))) # Ensures level ordering is correct
highorder <- order(with(GRSclust.count23,tapply(Mean,location_code_high,mean))) # Ensures level ordering is correct
levels(GRSclust.count23$location_code) <- list("1"=medorder[1], "2"=medorder[2])
levels(GRSclust.count23$location_code_high) <- list("1"=highorder[1], "2"=highorder[2], "3"=highorder[3])

with(GRSclust.count23,tapply(Mean,location_code,mean))
with(GRSclust.count23,tapply(Mean,location_code_high,mean))


common_name <- "Caribbean reef shark"
CRSclust.23 <- full_join(rownames_to_column(CRSclust.2), rownames_to_column(CRSclust.3), by = ("rowname" = "rowname"))
names(CRSclust.23)<- c('site_name','location_code','location_code_high')
CRSclust.count23 <- full_join(sumndata_site_CRS,CRSclust.23)
CRSclust.count23 <- cbind(common_name,CRSclust.count23)
CRSclust.count23$region_name <- ifelse(CRSclust.count23$region_name=="Western Atlantic","Western Atlantic","Western Pacific")
CRSclust.count23 <- CRSclust.count23 %>%
  group_by(region_name) %>%
  arrange(site_name) %>%
  mutate(location_code = factor(location_code),
         location_code_high = factor(location_code_high))
names(CRSclust.count23)<- c('common_name','site_name','region_name','Ndrops','Mean','SE','CI','location_code','location_code_high')
medorder <- order(with(CRSclust.count23,tapply(Mean,location_code,mean))) # Ensures level ordering is correct
highorder <- order(with(CRSclust.count23,tapply(Mean,location_code_high,mean))) # Ensures level ordering is correct
levels(CRSclust.count23$location_code) <- list("1"=medorder[1], "2"=medorder[2])
levels(CRSclust.count23$location_code_high) <- list("1"=highorder[1], "2"=highorder[2], "3"=highorder[3])

with(CRSclust.count23,tapply(Mean,location_code,mean))
with(CRSclust.count23,tapply(Mean,location_code_high,mean))

common_name <- "Nurse shark"
NSclust.23 <- full_join(rownames_to_column(NSclust.2), rownames_to_column(NSclust.3), by = ("rowname" = "rowname"))
names(NSclust.23)<- c('site_name','location_code','location_code_high')
NSclust.count23 <- full_join(sumndata_site_NS,NSclust.23)
NSclust.count23 <- cbind(common_name,NSclust.count23)
NSclust.count23$region_name <- ifelse(NSclust.count23$region_name=="Western Atlantic","Western Atlantic","Western Pacific")
NSclust.count23 <- NSclust.count23 %>%
  group_by(region_name) %>%
  arrange(site_name) %>%
  mutate(location_code = factor(location_code),
         location_code_high = factor(location_code_high))
names(NSclust.count23)<- c('common_name','site_name','region_name','Ndrops','Mean','SE','CI','location_code','location_code_high')
medorder <- order(with(NSclust.count23,tapply(Mean,location_code,mean))) # Ensures level ordering is correct
highorder <- order(with(NSclust.count23,tapply(Mean,location_code_high,mean))) # Ensures level ordering is correct
levels(NSclust.count23$location_code) <- list("1"=medorder[1], "2"=medorder[2])
levels(NSclust.count23$location_code_high) <- list("1"=highorder[1], "2"=highorder[2], "3"=highorder[3])

with(NSclust.count23,tapply(Mean,location_code,mean))
with(NSclust.count23,tapply(Mean,location_code_high,mean))


####

clust.count23 <- rbind(WRSclust.count23,
                       BRSclust.count23,
                       GRSclust.count23,
                       CRSclust.count23,
                       NSclust.count23)
clust.count23$common_name2 <- factor(clust.count23$common_name)
# Relevel how we would like to plot the data
clust.count23$common_name2 <- factor(clust.count23$common_name2, levels = c("Whitetip reef shark","Blacktip reef shark",
                                                                            "Grey reef shark","Caribbean reef shark",
                                                                            "Nurse shark"))
levels(clust.count23$common_name2) <- c("Whitetip", "Blacktip", "Grey",
                                       "Caribbean","Nurse")

# This is the order in which we want the sites plotted
site_name_order <- clust.count23 %>%
  #filter(region_name=='Western Atlantic') %>% 
  group_by(region_name)  %>% 
  distinct(site_name)

clust.count23$site_name <- factor(clust.count23$site_name, levels = as.data.frame(site_name_order)[,1])

plot.hl <- ggplot(clust.count23, 
                 aes(x = site_name, y = Mean, 
                     fill = factor(location_code, labels = c("Low","High")), 
                     legend=FALSE)) + 
  #coord_flip() +
  geom_bar(stat = "identity")+
  geom_errorbar(width=.1, aes(ymin=Mean-SE, ymax=Mean+SE), colour="black")+
  labs(title = "", x = "", y = "Max N", fill = "Abundance") +
  theme_classic()+
  facet_grid(common_name2 ~ .)+
  theme(axis.text.x = element_text(size=8,angle = 270, vjust = 0.5, hjust=0),
        legend.position=c(0,1), 
        legend.justification=c(0, 0),
        legend.direction="horizontal",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        #legend.position= "none",
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color="black", size = 0.3),
        axis.line.y = element_line(color="black", size = 0.3),
        axis.ticks = element_line(color="black", size = 0.3),       
        axis.ticks.length = unit(0.2, "cm"),
        axis.text=element_text(size=10,color="black"),  
        axis.title.y = element_text(size=12,margin=margin(0,16,0,0)),
        axis.title.x = element_text(size=12,margin=margin(10,0,0,0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        #text = element_text(size=10),
        strip.text.x = element_text(size = 12, colour = "orange"))
  #scale_alpha(guide = 'none')

## FIGURE for manuscript

tiff(filename = paste0("Images/Figure S2_Dispersal kernels_species.tiff"),
    #width = 1308, height = 877, units = "px")
    width = 300, height = 180, units = "mm",
    res=600)
plot.hl 
dev.off()

## Values probvided in manuscript
clust.count23

uniquesites <- clust.count23 %>%
  select(site_name,Ndrops)%>%
  distinct()
## How many drops?
sum(uniquesites$Ndrops)
##How many sites
nrow(uniquesites)



detach(package:plyr) # Need this to get the aggregating command working properly

# First in general for all sites
clust.count23 %>%
  group_by(common_name2) %>% 
  summarise(mean = mean(Mean),
            min = min(Mean),
            max = max(Mean))

# Second in split into our high vs low categories 
clust.count23 %>%
  group_by(location_code,common_name2) %>% 
  summarise(mean = mean(Mean),
            min = min(Mean),
            max = max(Mean))



##

## Redundant plot for high medium low abundance categories
plot.hml <- ggplot(clust.count23, 
                  aes(x = site_name, y = Mean, 
                      fill = factor(location_code_high, labels = c("Low","Medium","High")), 
                      legend=FALSE)) + 
  geom_bar(stat = "identity")+
  geom_errorbar(width=.1, aes(ymin=Mean-SE, ymax=Mean+SE), colour="black")+
  labs(title = "", x = "BRUV site", y = "Max N", fill = "Abundance") +
  theme_classic()+
  facet_grid(common_name2 ~ .)+
  theme(axis.text.x = element_text(size=8,angle = 90, vjust = 0.5, hjust=1),
        legend.position=c(0,1), 
        legend.justification=c(0, 0),
        legend.direction="horizontal",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        #legend.position= "none",
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color="black", size = 0.3),
        axis.line.y = element_line(color="black", size = 0.3),
        axis.ticks = element_line(color="black", size = 0.3),       
        axis.ticks.length = unit(0.2, "cm"),
        axis.text=element_text(size=10,color="black"),  
        axis.title.y = element_text(size=12,margin=margin(0,16,0,0)),
        axis.title.x = element_text(size=12,margin=margin(10,0,0,0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        #text = element_text(size=10),
        strip.text.x = element_text(size = 12, colour = "orange"))
#scale_alpha(guide = 'none')

plot.hml

################

# Recalibrate BRUV spreadsheet with sites sorted into new categories - high-low ; high-med-low
maxndata_spec1 <- fextract_sp(maxndata,sharksp[1]) # Whitetip reef shark
maxndata_spec2 <- fextract_sp(maxndata,sharksp[2]) # Blacktip reef shark
maxndata_spec3 <- fextract_sp(maxndata,sharksp[3]) # Grey reef shark
maxndata_spec4 <- fextract_sp(maxndata,sharksp[4]) # Caribbean reef shark
maxndata_spec5 <- fextract_sp(maxndata,sharksp[5]) # Nurse shark

clust.count23 <- rbind(WRSclust.count23, # Whitetip reef shark
                       BRSclust.count23, # Blacktip reef shark
                       GRSclust.count23, # Grey reef shark
                       CRSclust.count23, # Caribbean reef shark
                       NSclust.count23)  # Nurse shark

# Join the Species specific MaxN tables (with no duplicate counts) with the new HCA groupings
maxndata_WRS <- data.frame(clust.count23) %>%
  filter(common_name == "Whitetip reef shark") %>%
  select(site_name,location_code,location_code_high) %>%
  full_join(maxndata_spec1) %>%
  select(c(names(maxndata_spec1),"location_code","location_code_high"))
nrow(maxndata_WRS) # 3510 drops

maxndata_BRS <- data.frame(clust.count23) %>%
  filter(common_name == "Blacktip reef shark") %>%
  select(site_name,location_code,location_code_high) %>%
  full_join(maxndata_spec2) %>%
  select(c(names(maxndata_spec2),"location_code","location_code_high"))
nrow(maxndata_BRS) # 3136 drops

maxndata_GRS <- data.frame(clust.count23) %>%
  filter(common_name == "Grey reef shark") %>%
  select(site_name,location_code,location_code_high) %>%
  full_join(maxndata_spec3) %>%
  select(c(names(maxndata_spec3),"location_code","location_code_high"))
nrow(maxndata_GRS) # 3495 drops

maxndata_CRS <- data.frame(clust.count23) %>%
  filter(common_name == "Caribbean reef shark") %>%
  select(site_name,location_code,location_code_high) %>%
  full_join(maxndata_spec4) %>%
  select(c(names(maxndata_spec4),"location_code","location_code_high"))
nrow(maxndata_CRS) # 4295 drops

maxndata_NS <- data.frame(clust.count23) %>%
  filter(common_name == "Nurse shark") %>%
  select(site_name,location_code,location_code_high) %>%
  full_join(maxndata_spec5) %>%
  select(c(names(maxndata_spec5),"location_code","location_code_high"))
nrow(maxndata_NS) # 4231 drops

write.csv(maxndata_WRS,paste0(datadir, "maxndata_WRS.csv"))
write.csv(maxndata_BRS,paste0(datadir, "maxndata_BRS.csv"))
write.csv(maxndata_GRS,paste0(datadir, "maxndata_GRS.csv"))
write.csv(maxndata_CRS,paste0(datadir, "maxndata_CRS.csv"))
write.csv(maxndata_NS,paste0(datadir, "maxndata_NS.csv"))

############################################

