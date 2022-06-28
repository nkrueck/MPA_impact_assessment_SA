## Dispersal Kernel comparisons

library(tidyverse)
library(lubridate)

source("R/2018-10-15_displot.R")

pac_disp<-readRDS("data/Pacific/Dispersal Summaries/Pacific_dispersalSummary.RDS")
car_disp<-readRDS("data/Pacific/Dispersal Summaries/Caribbean_disperalSummary.RDS")

dispersal_summary<-
  bind_rows(
  pac_disp %>%
    mutate(Release.Date = date(Release.Date)), 
  car_disp) %>%
  filter(Consecutive.Dispersal > 0) %>%
  mutate(Sci.Name = as.factor(Sci.Name),
         Common.Name = as.factor(Common.Name))

## Time of each detection
dispersal_summary %>%
  group_by(Common.Name, Tag.ID) %>%
  summarise(
    # num.years = year(max(Date.Time)) - year(min(Date.Time)),
    num.days = difftime(date(max(Date.Time)), date(min(Date.Time)), units = "days"))

metadata <-
  dispersal_summary %>%
  filter(!is.na(Common.Name)) %>% 
  mutate(tag_id = Tag.ID,
         common_name=Common.Name,
         species = Sci.Name,
         sex = Sex,
         length_mm = Bio) %>%
  group_by(tag_id, common_name, species, sex, length_mm) %>%
  summarise(Longitude = mean(Longitude, na.rm=T),
            Latitude = mean(Latitude, na.rm=T)) %>% 
  mutate(region = 
           case_when(Longitude < 0 ~ "Caribbean/Atlantic",
                     TRUE ~ "Pacific"),
         installation =
           case_when(region %in% "Caribbean/Atlantic" & Latitude > 20 ~ "Florida",
                     region %in% "Caribbean/Atlantic" & Latitude < 20 & Latitude > 0 ~ "Belize",
                     region %in% "Caribbean/Atlantic" & Latitude < 0 ~ "Brazil",
                     region %in% "Pacific" & Longitude < 115 ~ "Ningaloo",
                     region %in% "Pacific" & Longitude > 115 & Longitude < 120 ~ "Rowley",
                     region %in% "Pacific" & Longitude > 120 & Longitude < 130 ~ "Scott",
                     region %in% "Pacific" & Longitude > 130 & Longitude < 149 ~ "FNQ",
                     region %in% "Pacific" & Longitude > 149 ~ "Heron"),
         array_area_m2 = 
           case_when(installation %in% "Belize" ~ 30790279,
                     installation %in% "Brazil" ~ 14074567,
                     installation %in% "Florida" ~ 110859213,
                     installation %in% "FNQ" ~ 104883087,
                     installation %in% "Heron" ~ 74363127,
                     installation %in% "Ningaloo" ~ 117306872,
                     installation %in% "Rowley" ~ 32621716,
                     installation %in% "Scott" ~ 24530935)
         ) %>% 
  dplyr::select(-c(Longitude, Latitude)) %>% 
  ungroup()

write_csv(metadata, "~/Desktop/2019-09-18_TagMetadata_installations.csv")

# array_area <-
#   dispersal_summary %>%
#   distinct(Longitude, Latitude) %>% 
#   mutate(region = 
#            case_when(Longitude < 0 ~ "Caribbean/Atlantic",
#                      TRUE ~ "Pacific"),
#          installation =
#            case_when(region %in% "Caribbean/Atlantic" & Latitude > 20 ~ "Florida",
#                      region %in% "Caribbean/Atlantic" & Latitude < 20 & Latitude > 0 ~ "Belize",
#                      region %in% "Caribbean/Atlantic" & Latitude < 0 ~ "Brazil",
#                      region %in% "Pacific" & Longitude < 115 ~ "Ningaloo",
#                      region %in% "Pacific" & Longitude > 115 & Longitude < 120 ~ "Rowley",
#                      region %in% "Pacific" & Longitude > 120 & Longitude < 130 ~ "Scott",
#                      region %in% "Pacific" & Longitude > 130 & Longitude < 149 ~ "FNQ",
#                      region %in% "Pacific" & Longitude > 149 ~ "Heron")) %>% 
#   st_as_sf(coords=c("Longitude", "Latitude"), crs=4326) %>% 
#   st_buffer(dist = 0.007) %>% 
#   group_by(installation) %>% 
#   summarise(geometry = st_union(geometry)) %>% 
#   ungroup()
# array_area %>% 
#   mutate(area = st_area(.))

mapview(array_area, zcol = "installation", burst = T)@map %>% 
  leaflet::addMeasure(primaryLengthUnit = "meters")



metadata

# ## Mapview of detection data
# library(raster)
# library(mapview)
# dplot<- car_disp
# coordinates(dplot)<-c("Longitude","Latitude")
# projection(dplot)<-CRS("+init=epsg:4326")
# m<-mapview(dplot, zcol="species", burst=T)
# mapshot(m, "Raw_Detections.html")
  
daily<-
  dispersal_summary %>%
  mutate(date = date(Date.Time),
         tag_id = Tag.ID,
         common_name=Common.Name,
         species = Sci.Name,
         sex = Sex,
         length_mm = Bio
         ) %>%
  group_by(date, tag_id, common_name, species, sex, length_mm) %>%
  summarise(max.dispersal_km = max(Consecutive.Dispersal, na.rm=T))

weekly<-
  dispersal_summary %>%
  mutate(week = format(Date.Time, "%Y-%W"),
         tag_id = Tag.ID,
         common_name=Common.Name,
         species = Sci.Name,
         sex = Sex,
         length_mm = Bio
  ) %>%
  group_by(week, tag_id, common_name, species, sex, length_mm) %>%
  summarise(max.dispersal_km = max(Consecutive.Dispersal, na.rm=T))


monthly<-
  dispersal_summary %>%
  mutate(month = format(Date.Time, "%Y-%m"),
         tag_id = Tag.ID,
         common_name=Common.Name,
         species = Sci.Name,
         sex = Sex,
         length_mm = Bio
  ) %>%
  group_by(month, tag_id, common_name, species, sex, length_mm) %>%
  summarise(max.dispersal_km = max(Consecutive.Dispersal, na.rm=T))

yearly<-
  dispersal_summary %>%
  mutate(year = year(Date.Time),
         tag_id = Tag.ID,
         common_name=Common.Name,
         species = Sci.Name,
         sex = Sex,
         length_mm = Bio
  ) %>%
  group_by(year, tag_id, common_name, species, sex, length_mm) %>%
  summarise(max.dispersal_km = max(Consecutive.Dispersal, na.rm=T))

Dispersal_ouput<-list(Daily=daily, Weekly=weekly, Monthly=monthly, Yearly=yearly)

saveRDS(Dispersal_ouput, "~/Desktop/Dispersal_Timescales.RDS")

