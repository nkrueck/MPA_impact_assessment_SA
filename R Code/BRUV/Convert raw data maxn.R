## MaxN standardization for BRUV observation frequencies
## Code will take Max Ns and plug in correct zeros for each site for our species of interest

library(ggplot2)
library(scales)
library(dplyr)
library(broom)

source('R Code/summarySE.R') # Loads bespoke functions for subsetting and plotting data

# Location where the raw acoustic telemetry and Bruv data are stored
datafolder <- "/Users/uqrdwye2/Dropbox/shark_mpa_model_v400/SA/DEW Marine Parks project/"

# load the BRUVS data
maxndata <- read.csv(paste0(datafolder,"BRUVS data/SA.Synthesis.MaxN.csv"))

# load the BRUVS metadata
maxnmeta <- read.csv(paste0(datafolder,"BRUVS data/SA.Synthesis.Metadata.csv"))

# Combine the data and metadata files
maxndata1 <- maxnmeta %>%
  dplyr::select(sample,latitude,longitude,location,status,site,depth) %>%
  right_join(maxndata) %>%
  dplyr::select(-family,-genus,-species,-id)

# Function to generate species specific bruv data with correct absences
fextract_sp <- function(data,singlesp){
  
  maxndata_spec <- data %>% 
    filter(genus.species ==singlesp)  # Extract data for single species
  
  all_drops <- maxnmeta %>% 
    dplyr::select(sample,latitude,longitude,location,status,site,depth,campaignid) # Only the columns we are interested in
  
  all_drops_sample <- unique(as.character(all_drops$sample))   # What are the names of the BRUV drops from all drops?
  
  # Which drops were in a site with one of our named species, but didn't see any of this species on a survey? 
  # I.e. generates the 'absences'
  nosinglesp_setids <- all_drops[all_drops_sample %!in% maxndata_spec$sample,'sample'] # Extract set_ids not present in the single species data (i.e. absences)
  maxndata_nospec <- all_drops %>% 
    filter(sample %in% nosinglesp_setids) %>% # filter absences dataset for these set codes
    mutate(maxn = 0, # Add the missing columns and the Zeros for these BRUV surveys
           genus.species = singlesp) %>%
    distinct()
    
  #Assign the correct zeros to the shark bruv data dataframe (bind our shark species, other animals and no animals dataframes)
  maxndata_spec <- rbind(maxndata_spec,maxndata_nospec) %>% arrange(sample) # Arrange by BRUV deployment id
  
  return(maxndata_spec)
}

# Run each species in a function be calling the species name
genus.sp <- c("Myliobatis tenuicaudatus", #Southern eagle ray (Age unknown)
             "Othos dentex", #Harlequin Fish (Age unknown)
             "Seriola lalandi", #Yellowtail Kingfish 
             "Chrysophrys auratus", #Snapper 
             "Achoerodus gouldii", #Western Blue Groper 
             "Notolabrus tetricus", #Bluethroat Wrasse 
             "Carcharhinus brachyurus", #Bronze Whaler 
             "Carcharhinus obscurus", #Dusky Whaler  (Missing from BRUV dataset)
             "Pseudocaranx spp", #Silver Trevally (Coded as this in the BRUVS dataset)
             "Carcharodon carcharias") # White Shark

# Run function to generate species specific bruv data with correct 'absences'
# Also patch on common names to align with the Dispersal kernels dataset
maxndata_spec1 <- fextract_sp(maxndata1,genus.sp[1]) # Southern Eagle Ray
maxndata_spec1 <- maxndata_spec1 %>% mutate(species_common_name="Southern Eagle Ray")
maxndata_spec2 <- fextract_sp(maxndata1,genus.sp[2]) # Harlequin Fish
maxndata_spec2 <- maxndata_spec2 %>% mutate(species_common_name="Harlequin")
maxndata_spec3 <- fextract_sp(maxndata1,genus.sp[3]) # Yellowtail Kingfish
maxndata_spec3 <- maxndata_spec3 %>% mutate(species_common_name="Yellowtail Kingfish")
maxndata_spec4 <- fextract_sp(maxndata1,genus.sp[4]) # Snapper
maxndata_spec4 <- maxndata_spec4 %>% mutate(species_common_name="Snapper")
maxndata_spec5 <- fextract_sp(maxndata1,genus.sp[5]) # Western Blue Groper 
maxndata_spec5 <- maxndata_spec5 %>% mutate(species_common_name="Western Blue Groper")
maxndata_spec6 <- fextract_sp(maxndata1,genus.sp[6]) # Bluethroat Wrasse 
maxndata_spec6 <- maxndata_spec6 %>% mutate(species_common_name="Bluethroat Wrasse")
maxndata_spec7 <- fextract_sp(maxndata1,genus.sp[7]) # Bronze Whaler
maxndata_spec7 <- maxndata_spec7 %>% mutate(species_common_name="Bronze Whaler")
maxndata_spec8 <- fextract_sp(maxndata1,genus.sp[8]) # Dusky Whaler 
maxndata_spec8 <- maxndata_spec8 %>% mutate(species_common_name="Dusky Whaler")
maxndata_spec9 <- fextract_sp(maxndata1,genus.sp[9]) # Silver Trevally
maxndata_spec9 <- maxndata_spec9 %>% mutate(species_common_name="Silver Trevally")
maxndata_spec10 <- fextract_sp(maxndata1,genus.sp[10]) # White Shark
maxndata_spec10 <- maxndata_spec10 %>% mutate(species_common_name="White Shark")

# Check raw data
check0 <- maxndata1 %>% filter(genus.species=="Chrysophrys auratus") #records with at least one of species seen 
dim(check0) # number of times species 'present'
sum(check0$maxn) # Total no animals seen

# Check against my data
check1 <- maxndata_spec4 %>% filter(maxn>0) #records with at least one of species seen 
head(check1,1)
dim(check1) # number of times species 'present'
sum(check1$maxn) # Total no animals seen

## Looks good!

# Combine into single dataset
maxndata_spec_all <- list(maxndata_spec1,
                           maxndata_spec2,
                           maxndata_spec3,
                           maxndata_spec4,
                           maxndata_spec5,
                           maxndata_spec6,
                           maxndata_spec7,
                           maxndata_spec8,
                           maxndata_spec9,
                           maxndata_spec10)

names(maxndata_spec_all) <- c("Southern_Eagle_Ray", 
                              "Harlequin",
                              "Yellowtail_Kingfish",
                              "Snapper",
                              "Western_Blue_Groper", 
                              "Bluethroat_Wrasse",
                              "Bronze_Whaler",
                              "Dusky_Whaler",
                              "Silver_Trevally",
                              "White_Shark")

saveRDS(maxndata_spec_all, file = "Data/MaxNs with absences.RDS") # Save to github

# Generate summarySE table for all 5 species combined
maxndata_spec_all.df <- do.call(rbind,maxndata_spec_all) #Also prepare a dataframe

maxndata_spec_all.df %>%
  summarySE(measurevar="maxn", groupvars=c("genus.species")) %>%
  dplyr::rename('meanmaxn'='maxn') %>%
ggplot(aes(x = genus.species, y = meanmaxn, alpha = 0.5,legend=FALSE)) + 
  geom_errorbar(width=.1, aes(ymin=meanmaxn-se, ymax=meanmaxn+se), colour="black")+
  labs(title = "Max N with absences", x = "Species", y = "Max N") +
  theme_bw() + 
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        legend.position = c(0.8, 0.8))+
  scale_alpha(guide = 'none')
ggsave("Images/MaxN_1.jpg")


