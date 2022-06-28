## Dispersal Kernel comparisons
## Code will take various detection datasets for each species and combine with the tag release metadata
## Updated function to only work with SA fish

#install_github("RossDwyer/VTrack", configure.args = "--with-proj-lib=/usr/local/lib/")

library(tidyverse)
library(lubridate)
library(VTrack)
library(remora)
library(sf)
library(sp)
library(ggpubr)
#source("R/2018-10-15_displot.R")

## Get the tag files ready and run the QC -------------

# Location where the raw acoustic telemetry and Bruv data are stored
datafolder <- "/Users/uqrdwye2/Dropbox/shark_mpa_model_v400/SA/DEW Marine Parks project/"

#Snapper
sp_det <- paste0(datafolder,"Snapper/IMOS_detections.csv")
sp_receivermet <- paste0(datafolder,"Snapper/IMOS_receiver_deployment_metadata.csv")
sp_tagmet <- paste0(datafolder,"Snapper/IMOS_transmitter_deployment_metadata.csv") # duplicate sensor tags (some) need to combine. Filter to Gulf stvincent - ignore NSW DPI and Vic
sp_meas <- paste0(datafolder,"Snapper/IMOS_animal_measurements.csv")

## specify files to QC - use supplied example .csv data
sp_files <- list(det = sp_det,
                 rmeta = sp_receivermet,
                 tmeta = sp_tagmet,
                 meas = sp_meas)
qc.out <- runQC(sp_files) # Run remora::runQC() to combine fields 
#saveRDS(qc.out, file = "Data/snapper_detQC.RDS") # Save to github

## Get the data ready for generating temporal dispersal metrics --------------

# Read in the detection dataset
#qc.out <- readRDS("Data/snapper_detQC.RDS") # Load detection data
d.qc <- grabQC(qc.out, what = "dQC") # Grab QC detection-only data

# As Snapper have multiple locations where tagged, lets make sure we only use unique tags from SA  
d.dplyr <- d.qc %>% 
  filter(tagging_project_name %in% "Gulf St Vincent monitoring") %>%
  dplyr::select(-transmitter_id) # Remove transmitter_id as we'll be using unique tag_ids as some (dual) sensor tags

# Add time catagories to location summary
location_summary <-  d.dplyr %>%
  mutate(transmitter_id = tag_id, # ensures the below code is compatible with using a unique tag_id (not the duplicate ping and sensor tag ids)
         species_scientific_name = as.factor(species_scientific_name),
         species_common_name = as.factor(species_common_name),
         Day = date(detection_datetime),
         Week = format(Day, "%Y-%W"),
         Month = format(Day, "%Y-%m")) %>%
  dplyr::select(transmitter_id,
                tagging_project_name,
                species_common_name,species_scientific_name,
                detection_datetime,
                installation_name,station_name,
                receiver_deployment_longitude, 
                receiver_deployment_latitude,
                Day,Week,Month)


### Days ------------
# Unite columns to get unique detections at receiver stations within a specified time interval
location_summary_day <- 
  location_summary %>%
  unite("z", species_common_name, transmitter_id, Day, remove = FALSE) %>%
  distinct(z,receiver_deployment_longitude,receiver_deployment_latitude) %>%
  arrange(z)

# Function to calculate great circle distances
fGCdist <- function(x) {
  tempsf <- SpatialPoints(location_summary_day[x, c("receiver_deployment_longitude","receiver_deployment_latitude")], 
                          proj4string = CRS("+proj=longlat +ellps=WGS84"))
  return(max(spDists(tempsf,longlat = TRUE)))
}
out2 <- tapply(1:nrow(location_summary_day), location_summary_day$z, fGCdist) # Run the function on out dayly dataset

dispersal_summary_day <- location_summary_day %>%
  group_by(z) %>%
  dplyr::summarize(n_stations=n()) %>%
  mutate(maxDistkm = round(as.numeric(as.vector(out2)),2)) %>%
  ungroup() %>%
  separate(z, c("species_common_name", "transmitter_id", "Day"), sep = "([._:])")

### Weeks ------------
# Unite columns to get unique detections at receiver stations within a specified time interval
location_summary_week <- 
  location_summary %>%
  unite("z", species_common_name, transmitter_id, Week, remove = FALSE) %>%
  distinct(z,receiver_deployment_longitude,receiver_deployment_latitude) %>%
  arrange(z)

# Function to calculate great circle distances
fGCdist <- function(x) {
  tempsf <- SpatialPoints(location_summary_week[x, c("receiver_deployment_longitude","receiver_deployment_latitude")], 
                          proj4string = CRS("+proj=longlat +ellps=WGS84"))
  return(max(spDists(tempsf,longlat = TRUE)))
}
out2 <- tapply(1:nrow(location_summary_week), location_summary_week$z, fGCdist) # Run the function on out weekly dataset

dispersal_summary_week <- location_summary_week %>%
  group_by(z) %>%
  dplyr::summarize(n_stations=n()) %>%
  mutate(maxDistkm = round(as.numeric(as.vector(out2)),2)) %>%
  ungroup() %>%
  separate(z, c("species_common_name", "transmitter_id", "Week"), sep = "([._:])")

### Months ------------
# Unite columns to get unique detections at receiver stations within a specified time interval
location_summary_month <- 
  location_summary %>%
  unite("z", species_common_name, transmitter_id, Month, remove = FALSE) %>%
  distinct(z,receiver_deployment_longitude,receiver_deployment_latitude) %>%
  arrange(z)

# Function to calculate great circle distances
fGCdist <- function(x) {
  tempsf <- SpatialPoints(location_summary_month[x, c("receiver_deployment_longitude","receiver_deployment_latitude")], 
                          proj4string = CRS("+proj=longlat +ellps=WGS84"))
  return(max(spDists(tempsf,longlat = TRUE)))
}
out2 <- tapply(1:nrow(location_summary_month), location_summary_month$z, fGCdist) # Run the function on out monthly dataset

dispersal_summary_month <- location_summary_month %>%
  group_by(z) %>%
  dplyr::summarize(n_stations=n()) %>%
  mutate(maxDistkm = round(as.numeric(as.vector(out2)),2)) %>%
  ungroup() %>%
  separate(z, c("species_common_name", "transmitter_id", "Month"), sep = "([._:])")

# Save the file as an RDS object
Dispersal_Timescales_Snapper <- list(Daily=dispersal_summary_day,
                                     Weekly=dispersal_summary_week,
                                     Monthly=dispersal_summary_month)
saveRDS(Dispersal_Timescales_Snapper, file = "Data/Dispersal_Timescales_Snapper.RDS") # Save to github

##############################################


# Compute a histogram of distance per month
disp.hist <- Dispersal_Timescales_Snapper$Monthly %>%
  ggplot(aes(maxDistkm)) + geom_histogram() +theme_bw()
# Compute a histogram of num receiver stations
stat.hist <- Dispersal_Timescales_Snapper$Monthly %>%
  ggplot(aes(n_stations)) + geom_histogram() +theme_bw()
ggarrange(disp.hist,stat.hist)
# 
# 
# lonlat <- "+proj=longlat +datum=WGS84 +no_defs"
# dispersal_summary.xy <- st_as_sf(dispersal_summary, 
#                                 coords = c("receiver_deployment_longitude", 
#                                            "receiver_deployment_latitude")) # Convert pts from a data frame into an sf object, using x and y as coordinates
# st_crs(dispersal_summary.xy) <- lonlat # Make sure that we assign the lon-lat CRS so that R knows how this object is projected
# 
# dispersal_summary.xy %>%
#   group_by(species_common_name, transmitter_id, Month) %>%
#   st_distance()
# 
# 
# ## Time of each detection
# dispersal_summary %>%
#   group_by(species_common_name, transmitter_id) %>%
#   summarise(
#     # num.years = year(max(Date.Time)) - year(min(Date.Time)),
#     num.days = difftime(date(max(detection_datetime)), 
#                         date(min(detection_datetime)), 
#                         units = "days")) %>%
#   arrange(desc(num.days))
# 
# ## Max distance between of daily detections
# dispersal_summary %>%
#   group_by(species_common_name, transmitter_id,date) %>%
#   summarise(
#     # num.years = year(max(Date.Time)) - year(min(Date.Time)),
#     num.days = difftime(date(max(detection_datetime)), 
#                         date(min(detection_datetime)), 
#                         units = "days")) %>%
#   arrange(desc(num.days))
# 
# 
# metadata <-
#   dispersal_summary %>%
#   filter(!is.na(species_common_name)) %>% 
#   mutate(transmitter_id = transmitter_id,
#          species_common_name=species_common_name,
#          species_scientific_name = species_scientific_name) %>%#,
#          #sex = Sex,
#          #length_mm = Bio) %>%
#   group_by(transmitter_id, species_common_name, 
#            species_scientific_name)%>% #, sex, length_mm) %>%
#   summarise(receiver_deployment_longitude = mean(receiver_deployment_longitude, na.rm=T),
#             receiver_deployment_latitude = mean(receiver_deployment_latitude, na.rm=T)) #%>% 
#   # mutate(region = 
#   #          case_when(Longitude < 0 ~ "Caribbean/Atlantic",
#   #                    TRUE ~ "Pacific"),
#          # installation =
#          #   case_when(region %in% "Caribbean/Atlantic" & Latitude > 20 ~ "Florida",
#          #             region %in% "Caribbean/Atlantic" & Latitude < 20 & Latitude > 0 ~ "Belize",
#          #             region %in% "Caribbean/Atlantic" & Latitude < 0 ~ "Brazil",
#          #             region %in% "Pacific" & Longitude < 115 ~ "Ningaloo",
#          #             region %in% "Pacific" & Longitude > 115 & Longitude < 120 ~ "Rowley",
#          #             region %in% "Pacific" & Longitude > 120 & Longitude < 130 ~ "Scott",
#          #             region %in% "Pacific" & Longitude > 130 & Longitude < 149 ~ "FNQ",
#          #             region %in% "Pacific" & Longitude > 149 ~ "Heron"),
#          # array_area_m2 = 
#          #   case_when(installation %in% "Belize" ~ 30790279,
#          #             installation %in% "Brazil" ~ 14074567,
#          #             installation %in% "Florida" ~ 110859213,
#          #             installation %in% "FNQ" ~ 104883087,
#          #             installation %in% "Heron" ~ 74363127,
#          #             installation %in% "Ningaloo" ~ 117306872,
#          #             installation %in% "Rowley" ~ 32621716,
#          #             installation %in% "Scott" ~ 24530935)
#          # ) %>% 
#   dplyr::select(-c(Longitude, Latitude)) %>% 
#   ungroup()
# 
# write_csv(metadata, "~/Desktop/2019-09-18_TagMetadata_installations.csv")
# 
# # array_area <-
# #   dispersal_summary %>%
# #   distinct(Longitude, Latitude) %>% 
# #   mutate(region = 
# #            case_when(Longitude < 0 ~ "Caribbean/Atlantic",
# #                      TRUE ~ "Pacific"),
# #          installation =
# #            case_when(region %in% "Caribbean/Atlantic" & Latitude > 20 ~ "Florida",
# #                      region %in% "Caribbean/Atlantic" & Latitude < 20 & Latitude > 0 ~ "Belize",
# #                      region %in% "Caribbean/Atlantic" & Latitude < 0 ~ "Brazil",
# #                      region %in% "Pacific" & Longitude < 115 ~ "Ningaloo",
# #                      region %in% "Pacific" & Longitude > 115 & Longitude < 120 ~ "Rowley",
# #                      region %in% "Pacific" & Longitude > 120 & Longitude < 130 ~ "Scott",
# #                      region %in% "Pacific" & Longitude > 130 & Longitude < 149 ~ "FNQ",
# #                      region %in% "Pacific" & Longitude > 149 ~ "Heron")) %>% 
# #   st_as_sf(coords=c("Longitude", "Latitude"), crs=4326) %>% 
# #   st_buffer(dist = 0.007) %>% 
# #   group_by(installation) %>% 
# #   summarise(geometry = st_union(geometry)) %>% 
# #   ungroup()
# # array_area %>% 
# #   mutate(area = st_area(.))
# 
# mapview(array_area, zcol = "installation", burst = T)@map %>% 
#   leaflet::addMeasure(primaryLengthUnit = "meters")
# 
# 
# 
# metadata
# 
# # ## Mapview of detection data
# # library(raster)
# # library(mapview)
# # dplot<- car_disp
# # coordinates(dplot)<-c("Longitude","Latitude")
# # projection(dplot)<-CRS("+init=epsg:4326")
# # m<-mapview(dplot, zcol="species", burst=T)
# # mapshot(m, "Raw_Detections.html")
#   
# daily<-
#   dispersal_summary %>%
#   mutate(date = date(detection_datetime),
#          transmitter_id = transmitter_id,
#          species_common_name=species_common_name,
#          species_scientific_name = species_scientific_name#,
#          # sex = Sex,
#          # length_mm = Bio
#          ) %>%
#   group_by(date, transmitter_id, species_common_name, species_scientific_name) %>%#, sex, length_mm) %>%
#   #summarise(max.dispersal_km = max(Consecutive.Dispersal, na.rm=T))
#   summarise(max.dispersal_km = max(Distance_QC, na.rm=T),
#             consecutive.dispersal_km = sum(Distance_QC, na.rm=T))
#   
# daily.un <- head(data.frame(daily),2)
# 
# weekly<-
#   dispersal_summary %>%
#   mutate(week = format(Date.Time, "%Y-%W"),
#          tag_id = Tag.ID,
#          common_name=Common.Name,
#          species = Sci.Name,
#          sex = Sex,
#          length_mm = Bio
#   ) %>%
#   group_by(week, tag_id, common_name, species, sex, length_mm) %>%
#   summarise(max.dispersal_km = max(Consecutive.Dispersal, na.rm=T))
# 
# 
# monthly<-
#   dispersal_summary %>%
#   mutate(month = format(Date.Time, "%Y-%m"),
#          tag_id = Tag.ID,
#          common_name=Common.Name,
#          species = Sci.Name,
#          sex = Sex,
#          length_mm = Bio
#   ) %>%
#   group_by(month, tag_id, common_name, species, sex, length_mm) %>%
#   summarise(max.dispersal_km = max(Consecutive.Dispersal, na.rm=T))
# 
# yearly<-
#   dispersal_summary %>%
#   mutate(year = year(Date.Time),
#          tag_id = Tag.ID,
#          common_name=Common.Name,
#          species = Sci.Name,
#          sex = Sex,
#          length_mm = Bio
#   ) %>%
#   group_by(year, tag_id, common_name, species, sex, length_mm) %>%
#   summarise(max.dispersal_km = max(Consecutive.Dispersal, na.rm=T))
# 
# Dispersal_ouput<-list(Daily=daily, Weekly=weekly, Monthly=monthly, Yearly=yearly)
# 
# saveRDS(Dispersal_ouput, "~/Desktop/Dispersal_Timescales.RDS")

