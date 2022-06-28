### Preprocessing Caribbean data

library(tidyverse)
library(lubridate)
library(ATT)
library(mapview)

######################
## Demians data
statinfo<-read_csv("data/Caribbean/Caribbean Reef Sharks/Demian/statinfo_car_demian.csv")
tag<-read_csv("data/Caribbean/Caribbean Reef Sharks/Demian/tags_car_demian.csv")
taginfo<-read_csv("data/Caribbean/Caribbean Reef Sharks/Demian/taginfo_car_demian.csv")

tagdata<-
  tag %>%
  mutate(Date.and.Time..UTC. = lubridate::ymd_hms(paste(date, time)),
         Transmitter = paste0("A69-1000-", Transmitter.Name)) %>%
  dplyr::select(-c(date, time, Latitude, Longitude, Station.Name)) %>%
  left_join(filter(statinfo[c("receiver_name","station_name","station_longitude","station_latitude")], receiver_name < 24),
            by=c("Receiver" = "receiver_name")) %>%
  rename(Longitude = station_longitude,
         Latitude = station_latitude,
         Station.Name = station_name) %>%
  filter(!is.na(Latitude))


ATTdata<-setupData(Tag.Detections = tagdata,
                   Tag.Metadata = taginfo,
                   Station.Information = statinfo,
                   source = "VEMCO")

disp_demian<-dispersalSummary(ATTdata)

max.daily.disp_demian <-
  disp_demian %>%
  filter(Consecutive.Dispersal > 0) %>%
  mutate(date = date(Date.Time)) %>%
  group_by(date, Tag.ID) %>%
  summarize(common_name = first(Common.Name),
            species = first(Sci.Name),
            sex = first(Sex),
            length_mm = NA,
            max.dispersal_km = max(Consecutive.Dispersal, na.rm = TRUE)) %>%
  rename(tag_id = Tag.ID)

## output file
write.csv(max.daily.disp_demian, "data/Caribbean/Max.Daily.Dispersal_Cperezi_Demian.csv", row.names = F)

######################
## Katies data
statinfo<-read_csv("data/Caribbean/Caribbean Reef Sharks/Katie/statinfo_car_katie.csv")
tag<-read_csv("data/Caribbean/Caribbean Reef Sharks/Katie/Shark acoustic data merged.csv")
taginfo<-read_csv("data/Caribbean/Caribbean Reef Sharks/Katie/taginfo_car_katie.csv")

tagdata<-
  tag %>%
  mutate(Date.and.Time.UTC = lubridate::ymd_hm(Date.and.Time.UTC)) %>%
  left_join(taginfo, by=c("Transmitter" = "transmitter_id")) %>%
  left_join(statinfo, by=c("Receiver" = "receiver_name")) %>%
  filter(Date.and.Time.UTC > deploymentdatetime_timestamp & Date.and.Time.UTC < recoverydatetime_timestamp) %>%
  filter(common_name %in% "Caribbean Reef Shark") %>%
  transmute(Date.and.Time..UTC. = Date.and.Time.UTC,
            Receiver = Receiver,
            Transmitter = Transmitter,
            Transmitter.Name = Transmitter,
            Transmitter.Serial = Transmitter,
            Sensor.Value = NA,
            Sensor.Unit = NA,
            Station.Name = station_name,
            Longitude = station_longitude,
            Latitude = station_latitude) %>%
  filter(!is.na(Latitude))

ATTdata<-setupData(Tag.Detections = tagdata,
                   Tag.Metadata = taginfo,
                   Station.Information = statinfo,
                   source = "VEMCO")

disp_katie<-dispersalSummary(ATTdata)

max.daily.disp_katie <-
  disp_katie %>%
  filter(Consecutive.Dispersal > 0) %>%
  mutate(date = date(Date.Time)) %>%
  group_by(date, Tag.ID) %>%
  summarize(common_name = first(Common.Name),
            species = first(Sci.Name),
            sex = first(Sex),
            length_mm = NA,
            max.dispersal_km = max(Consecutive.Dispersal, na.rm = TRUE)) %>%
  rename(tag_id = Tag.ID)

## output file
write.csv(max.daily.disp_katie , "data/Caribbean/Max.Daily.Dispersal_Cperezi_Katie.csv", row.names = F)

######################
## Ricardos data
## 2001-03
statinfo<-read_csv("data/Caribbean/Caribbean Reef Sharks/Ricardo/2001-03/statinfo01-03_car_ricardo.csv")
tag<-read_csv("data/Caribbean/Caribbean Reef Sharks/Ricardo/2001-03/tag2001-03.csv")
taginfo<-read_csv("data/Caribbean/Caribbean Reef Sharks/Ricardo/2001-03/taginfo01-03_car_ricardo.csv")

tagdata<-
  tag %>%
  left_join(taginfo, by=c("Transmitter" = "transmitter_id")) %>%
  left_join(statinfo, by=c("Receiver" = "station_name")) %>%
  transmute(Date.and.Time..UTC. = Date.Time,
            Receiver = Receiver,
            Transmitter = Transmitter,
            Transmitter.Name = Transmitter,
            Transmitter.Serial = Transmitter,
            Sensor.Value = NA,
            Sensor.Unit = NA,
            Station.Name = Receiver,
            Longitude = station_longitude,
            Latitude = station_latitude) %>%
  filter(!is.na(Latitude))


ATTdata<-setupData(Tag.Detections = tagdata,
                   Tag.Metadata = taginfo,
                   Station.Information = statinfo,
                   source = "VEMCO")

disp_ric0103<-dispersalSummary(ATTdata)

max.daily.disp_0103 <-
  disp_ric0103 %>%
  filter(Consecutive.Dispersal > 0) %>%
  mutate(date = date(Date.Time)) %>%
  group_by(date, Tag.ID) %>%
  summarize(common_name = first(Common.Name),
            species = first(Sci.Name),
            sex = first(Sex),
            length_mm = NA,
            max.dispersal_km = max(Consecutive.Dispersal, na.rm = TRUE)) %>%
  rename(tag_id = Tag.ID)

## output file
write.csv(max.daily.disp_0103, "data/Caribbean/Max.Daily.Dispersal_Cperezi_Ric0103.csv", row.names = F)


## 2011-12
statinfo<-read_csv("data/Caribbean/Caribbean Reef Sharks/Ricardo/2011-12/statinfo11-12_car_ricardo.csv")
taginfo<-read_csv("data/Caribbean/Caribbean Reef Sharks/Ricardo/2011-12/taginfo11-12_car_ricardo.csv")

library(readxl)
tag<-read_excel("data/Caribbean/Caribbean Reef Sharks/Ricardo/2011-12/FEN_2011_12.xlsx", sheet=1)

id<-c("S10","S13","S14","S17","S2!","S22","S23","S25","S26","S28","S3!","S30","S31","S32","S33","S34","S35","S36","S5","S6","S9")

for(i in 1:length(id)){
 colnum<-grep(id[i], names(tag))
 if(i == 1){
   out<-tag[colnum]
   names(out)<-c("Date.and.Time..UTC.","Receiver")
   out$Transmitter<-id[i]
 }else{
   dd<-tag[colnum]
   names(dd)<-c("Date.and.Time..UTC.","Receiver")
   dd$Transmitter<-id[i]
   out<-rbind(out,dd)
}}


tagdata<-
  out %>%
  left_join(taginfo, by=c("Transmitter" = "tag_id")) %>%
  left_join(statinfo, by=c("Receiver" = "station_name")) %>%
  transmute(Date.and.Time..UTC. = Date.and.Time..UTC.,
            Receiver = Receiver,
            Transmitter = transmitter_id,
            Transmitter.Name = Transmitter,
            Transmitter.Serial = transmitter_id,
            Sensor.Value = NA,
            Sensor.Unit = NA,
            Station.Name = Receiver,
            Longitude = station_longitude,
            Latitude = station_latitude) %>%
  filter(!is.na(Latitude))


ATTdata<-setupData(Tag.Detections = tagdata,
                   Tag.Metadata = taginfo,
                   Station.Information = statinfo,
                   source = "VEMCO")

disp_ric1112<-dispersalSummary(ATTdata)

max.daily.disp_1112 <-
  disp_ric1112 %>%
  filter(Consecutive.Dispersal > 0) %>%
  mutate(date = date(Date.Time)) %>%
  group_by(date, Tag.ID) %>%
  summarize(common_name = first(Common.Name),
            species = first(Sci.Name),
            sex = first(Sex),
            length_mm = NA,
            max.dispersal_km = max(Consecutive.Dispersal, na.rm = TRUE)) %>%
  rename(tag_id = Tag.ID)

## output file
write.csv(max.daily.disp_1112, "data/Caribbean/Max.Daily.Dispersal_Cperezi_Ric1112.csv", row.names = F)


################
## Nurse Shark
library(readxl)
path<-"data/Caribbean/Nurse Sharks/Wes Pratt/FWC_NurseSharks.xlsx"
sheetnames <- excel_sheets(path)
mylist <- lapply(excel_sheets(path), read_excel, path = path)
tag <- do.call("rbind", mylist)
  
tagdata <-
  tag %>%
  filter(!is.na(`Station Latitude`)) %>%
  transmute(
    Date.and.Time..UTC. = `Date/Time`,
    Receiver = `Receiver Name`,
    Transmitter = paste(`Code Space`,ID, sep="-"),
    Transmitter.Name = ID,
    Transmitter.Serial = ID,
    Sensor.Value = NA,
    Sensor.Unit = NA,
    Station.Name = `Station Name`,
    Longitude = `Station Longitude`,
    Latitude = `Station Latitude`
  )

tagdata[tagdata$Longitude >0,"Longitude"] <- tagdata[tagdata$Longitude >0,"Longitude"]*-1
tagdata[tagdata$Latitude < 0,"Latitude"] <- tagdata[tagdata$Latitude <0,"Latitude"] * -1

taginfo <-
  tag %>%
  filter(!is.na(`Station Latitude`)) %>%
  group_by(ID) %>%
  summarise(transmitter_id = paste(first(`Code Space`),first(ID), sep="-"),
            release_id = first(ID),
            tag_project_name = "FWC_NurseSharks",
            scientific_name = "Ginglymostoma cirratum",
            common_name = "Nurse Shark",
            release_latitude = NA,
            release_longitude = NA,
            ReleaseDate = NA,
            tag_expected_life_time_days = NA,
            tag_status = "deployed",
            sex = NA,
            measurement = NA) %>%
  rename(tag_id = ID)

statinfo<-
  tag %>%
  filter(!is.na(`Station Latitude`)) %>%
  group_by(`Station Name`) %>%
  summarise(project_name = "FWC_NurseShark",
            installation_name = "Caribbean",
            receiver_name = first(`Receiver Name`),
            deploymentdatetime_timestamp = NA,
            recoverydatetime_timestamp = NA,
            status = "deployed",
            station_latitude = first(`Station Latitude`),
            station_longitude = first(`Station Longitude`),
            imos_device = F) %>%
  rename(station_name = `Station Name`)
  
statinfo[statinfo$station_longitude >0,"station_longitude"] <- statinfo[statinfo$station_longitude >0,"station_longitude"]*-1
statinfo[statinfo$station_latitude < 0,"station_latitude"] <- statinfo[statinfo$station_latitude <0,"station_latitude"] * -1

ATTdata<-setupData(Tag.Detections = tagdata,
                   Tag.Metadata = taginfo,
                   Station.Information = statinfo,
                   source = "VEMCO")

disp_wes<-dispersalSummary(ATTdata)

max.daily.disp_wes <-
  disp_wes %>%
  filter(Consecutive.Dispersal > 0) %>%
  mutate(date = date(Date.Time)) %>%
  group_by(date, Tag.ID) %>%
  summarize(common_name = first(Common.Name),
            species = first(Sci.Name),
            sex = first(Sex),
            length_mm = NA,
            max.dispersal_km = max(Consecutive.Dispersal, na.rm = TRUE)) %>%
  rename(tag_id = Tag.ID)

## output file
write.csv(max.daily.disp_wes, "data/Caribbean/Max.Daily.Dispersal_NurseShark.csv", row.names = F)

######### 
# Merge

## combine raw dispersal output
Caribbean_dispersal<-
  bind_rows(
    disp_demian %>%
      mutate(
        Transmitter = as.character(Transmitter),
        Release.Latitude = as.numeric(Release.Latitude),
        Release.Longitude = as.numeric(Release.Longitude),
        Receiver = as.character(Receiver),
        Bio = as.numeric(Bio),
        source = "demian"
      ),
    disp_katie %>%
      mutate(
        Transmitter = as.character(Transmitter),
        Tag.ID = as.character(Tag.ID),
        Receiver = as.character(Receiver),
        source = "katie"
      ),
    disp_ric0103 %>%
      mutate(
        Transmitter = as.character(Transmitter),
        Release.Latitude = as.numeric(Release.Latitude),
        Release.Longitude = as.numeric(Release.Longitude),
        source = "ric0103"
      ),
    disp_ric1112 %>%
      mutate(
        Transmitter = as.character(Transmitter),
        Station.Name = as.character(Station.Name),
        Receiver = as.character(Receiver),
        Release.Latitude = as.numeric(Release.Latitude),
        Release.Longitude = as.numeric(Release.Longitude),
        source = "ric1112"
      ),
    disp_wes %>%
      mutate(Tag.ID = as.character(Tag.ID))
  ) %>%
  mutate(Common.Name = 
           case_when(
             Sci.Name %in% "Carcharhinus perezi" ~ "Caribbean Reef Shark",
             Sci.Name %in% "Ginglymostoma cirratum" ~ "Nurse Shark"
           ))

saveRDS(Caribbean_dispersal, "data/Dispersal Summaries/Caribbean_disperalSummary.RDS")

car1<-read_csv("data/Caribbean/Max.Daily.Dispersal_Cperezi_Demian.csv")
car2<-read_csv("data/Caribbean/Max.Daily.Dispersal_Cperezi_Katie.csv")
car3<-read_csv("data/Caribbean/Max.Daily.Dispersal_Cperezi_Ric0103.csv")
car4<-read_csv("data/Caribbean/Max.Daily.Dispersal_Cperezi_Ric1112.csv")

nur<-read_csv("data/Caribbean/Max.Daily.Dispersal_NurseShark.csv")

dispCar<-do.call("rbind", list(car1,car2,car3,car4,nur))
write_csv(dispCar, "data/Max.Daily.Dispersal_Caribbean.csv")


## Displots
source("R/displot.R")

disp_Pacific<-read_csv("data/Max.Daily.Dispersal_Pacific.csv")

keepPac<-
  disp_Pacific %>%
  group_by(common_name, tag_id) %>%
  summarise(n = n()) %>%
  filter(n > 1)

disp_Caribbean<-read_csv("data/Max.Daily.Dispersal_Caribbean.csv")

keepCar<-
  disp_Caribbean %>%
  group_by(common_name, tag_id) %>%
  summarise(n = n()) %>%
  filter(n > 1)

dispPac<-
  disp_Pacific %>%
  filter(tag_id %in% unique(keepPac$tag_id))

dispCar<-
  disp_Caribbean %>%
  filter(tag_id %in% unique(keepCar$tag_id))

quartz(width=6, height=4)
displot(dispPac, cn="Blacktip Reef Shark", var="max.dispersal_km", plotit=T, lcol="steelblue4", lty=1, bars=F, ylim=c(0,1), lab=F, xlab="Maximum Daily Dispersal distance (km)", stand=T)
displot(dispPac, cn="Whitetip Reef Shark", var="max.dispersal_km", plotit=T, lcol="steelblue4", lty=2, bars=F, ylim=c(0,1), add=T, lab=F, stand=T)
displot(dispPac, cn="Grey Reef Shark", var="max.dispersal_km", plotit=T, lcol="steelblue4", lty=3, bars=F, ylim=c(0,1), add=T, lab=F, stand=T)

displot(dispCar, cn="Caribbean Reef Shark", var="max.dispersal_km", plotit=T, bars=F, ylim=c(0,1), add=T, lab=F, lcol="coral", stand=T)
displot(dispCar, cn="Nurse Shark", var="max.dispersal_km", plotit=T, bars=F, ylim=c(0,1), add=T, lab=F, lcol="coral", lty=2, stand=T)

legend('topright', lty=c(1,2,3,NA,1,2), col=c(rep("steelblue4",3),NA,rep("coral",2)), 
       legend=c("Blacktip Reef Shark","Whitetip Reef Shark","Grey Reef Shark",
                NA, "Caribbean Reef Shark","Nurse Shark"), bty="n", lwd=2)






