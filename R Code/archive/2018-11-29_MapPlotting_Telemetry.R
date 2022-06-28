## Telemetry maps for STAR methods.

library(mapview)
library(leaflet)
library(tidyverse)
library(sf)
library(sp)
library(ggpubr)


### Input data
maxndata <- read_csv("data/BRUVS/maxndata.csv")
clusters <- read_csv("data/BRUVS/kmeans_clusters_lowhigh.csv")

plotdat <- 
  maxndata %>%
  filter(common_name %in% c("Whitetip reef shark",
                            "Blacktip reef shark",
                            "Grey reef shark",
                            "Caribbean reef shark",
                            "Nurse shark")) %>%
  group_by(site_name) %>%
  summarise(mean_lat = mean(set_lat, na.rm=T),
            mean_lon = mean(set_long, na.rm=T),
            min_maxn = min(maxn, na.rm=T),
            max_maxn = max(maxn, na.rm=T),
            mean_maxn = mean(maxn, na.rm=T),
            se_maxn = sd(maxn, na.rm=T)/sqrt(length(maxn))) %>%
  full_join(clusters) %>%
  mutate(mean_lon =
           case_when(
             mean_lon < -110 ~ mean_lon + 360,
             mean_lon > -110 ~ mean_lon
           ),
         region = 
           case_when(
             mean_lon < 70 ~ "Atlantic",
             mean_lon > 70 ~ "Pacific"
           ))

Atl<-
  plotdat %>%
  filter(region %in% "Atlantic") %>%
  st_as_sf(coords=c("mean_lon","mean_lat"), crs=4326) %>%
  mapview(cex = "mean_maxn", zcol="mean_maxn", legend = T, layer.name = "Mean MaxN")

Atl  
# mapshot(Atl, file="AtlanticBRUVS.png", remove_controls = c("zoomControl", "layersControl", "homeButton"))


Pac<-
  plotdat %>%
  filter(region %in% "Pacific") %>%
  st_as_sf(coords=c("mean_lon","mean_lat"), crs=4326) %>%
  mapview(cex = "mean_maxn", zcol="mean_maxn", legend = T, layer.name = "Mean MaxN")

Pac
# mapshot(Pac, file="PacificBRUVS.png", remove_controls = c("zoomControl", "layersControl", "homeButton"), vwidth=800, vheight=550)

### Plotting the clusters for each species
plotdat <- 
  maxndata %>%
  filter(common_name %in% c("Whitetip reef shark",
                            "Blacktip reef shark",
                            "Grey reef shark",
                            "Caribbean reef shark",
                            "Nurse shark")) %>%
  group_by(site_name, common_name) %>%
  summarise(mean_lat = mean(set_lat, na.rm=T),
            mean_lon = mean(set_long, na.rm=T),
            mean_maxn = mean(maxn, na.rm=T)) %>%
  full_join(
    clusters %>%
      gather(common_name, cluster, Whitetip.reef.shark:Nurse.shark) %>%
      mutate(common_name = gsub("[.]", " ", common_name))
    ) %>%
  mutate(mean_lon =
           case_when(
             mean_lon < -110 ~ mean_lon + 360,
             mean_lon > -110 ~ mean_lon
           ),
         region =
           case_when(
             mean_lon < 70 ~ "Atlantic",
             mean_lon > 70 ~ "Pacific"
           ),
         Abundance = factor(
           case_when (
             cluster %in% 1 ~ "Low",
             cluster %in% 2 ~ "High"
           ), levels=c("Low", "High"))
         ) %>% 
  ungroup()

### Correct cluster classification for Whitetip and Caribbean

plotdat[plotdat$common_name %in% "Whitetip reef shark",]$Abundance <- fct_recode(plotdat[plotdat$common_name %in% "Whitetip reef shark",]$Abundance, High = "Low", Low = "High")
plotdat[plotdat$common_name %in% "Caribbean reef shark",]$Abundance <- fct_recode(plotdat[plotdat$common_name %in% "Caribbean reef shark",]$Abundance, High = "Low", Low = "High")


## ggplot but with grey land and global projection
library(rworldmap)
library(rworldxtra)
library(cowplot)

#Pacific
world2 <- map_data("world2")
base_pac <-
  ggplot() + 
  geom_map(data=world2, map=world2, aes(x=long, y=lat, map_id=region), fill="#CCCCCC") +
  coord_map("ortho",
                 ylim=c(-40,40),
                 orientation = c(-20,  175.0, 2)) + # -9.40, 183.0
  scale_x_continuous(breaks = seq(-360, 360, by=10)) +
  scale_y_continuous(breaks = seq(90,-90, by=-10)) +
  geom_hline(yintercept=0, size = rel(0.5)) +
  theme_map() %+replace%
  theme(panel.grid.major = element_line(size = rel(0.05)),
        legend.position = "bottom",
        legend.box.spacing = unit(3, "lines"),
        # strip.text = element_blank(),
        panel.spacing = unit(2, "lines"))
  
pacific <-
  base_pac +
  geom_point(
    aes(
      x = mean_lon,
      y = mean_lat,
      color = I("white"),
      fill = Abundance,
      size = mean_maxn
    ),
    shape = 21,
    data =
      plotdat %>%
      filter(!common_name %in% c("Caribbean reef shark", "Nurse shark")) %>%
      filter(!is.na(Abundance)) %>%
      droplevels()
  ) +
  scale_size_continuous(limits=c(1,3),breaks=seq(1,3, by=0.5)) +
  facet_wrap(~ common_name)

pacific

#Caribbean
world <- map_data("world")
base_atl <-
  ggplot() + 
  geom_map(data=world, map=world, aes(x=long, y=lat, map_id=region), fill="#CCCCCC") +
  coord_map("ortho",
            ylim=c(-40,40),
            orientation = c(10, -65, 0)) + # -9.40, 183.0
  scale_x_continuous(breaks = seq(-360, 360, by=10)) +
  scale_y_continuous(breaks = seq(90,-90, by=-10)) +
  geom_hline(yintercept=0, size = rel(0.5)) +
  theme_map() %+replace%
  theme(panel.grid.major = element_line(size = rel(0.05)),
        legend.position = "bottom",
        legend.box.spacing = unit(3, "lines"),
        # strip.text = element_blank(),
        panel.spacing = unit(2, "lines"))

atlantic <-
  base_atl +
  geom_point(
    aes(
      x = mean_lon,
      y = mean_lat,
      color = I("white"),
      fill = Abundance,
      size = mean_maxn
    ),
    shape = 21,
    data =
      plotdat %>%
      filter(common_name %in% c("Caribbean reef shark", "Nurse shark")) %>%
      filter(!is.na(Abundance)) %>%
      droplevels()
  ) +
  scale_size_continuous(limits=c(1,3),breaks=seq(1,3, by=0.5)) +
  facet_wrap( ~ common_name)

atlantic

## Save separately
ggsave("PacificMaps.pdf", plot = pacific,  width = 10, height = 3, units ="in")
ggsave("AtlanticMaps.pdf", plot = atlantic,  width = 10, height = 3, units ="in")


## Put together in one plot
both <- ggarrange(pacific, atlantic, nrow=2, common.legend = T)
both  
## for some reason ggarrange removes minor graticules!

ggsave("AbundanceMaps.pdf", plot=both, width = 10, height = 7, units="in")

ggsave("AbundanceMaps.png", plot=both, width = 20, height = 14, units="in")


## mapview version
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## Customised leaflet legend
addLegendCustom <- function(map, colors, labels, sizes, shapes, borders, opacity = 0.5){

  make_shapes <- function(colors, sizes, borders, shapes) {
    shapes <- gsub("circle", "50%", shapes)
    shapes <- gsub("square", "0%", shapes)
    paste0(colors, "; width:", sizes, "px; height:", sizes, "px; border:1px solid ", borders, "; border-radius:", shapes)
  }
  make_labels <- function(sizes, labels) {
    paste0("<div style='display: inline-block;height: ",
           sizes, "px;margin-top: 4px;line-height: ",
           sizes, "px;'>", labels, "</div>")
  }

  legend_colors <- make_shapes(colors, sizes*2, borders, shapes)
  legend_labels <- make_labels(sizes, labels)

  return(addLegend(map, colors = legend_colors, labels = legend_labels, opacity = opacity))
}

whitetip<- plotdat %>% filter(region %in% "Pacific", common_name %in% "Whitetip reef shark", !is.na(Abundance))
blacktip<- plotdat %>% filter(region %in% "Pacific", common_name %in% "Blacktip reef shark", !is.na(Abundance))
greyreef<- plotdat %>% filter(region %in% "Pacific", common_name %in% "Grey reef shark", !is.na(Abundance))
caribbean<- plotdat %>% filter(region %in% "Atlantic", common_name %in% "Caribbean reef shark", !is.na(Abundance))
nurse<- plotdat %>% filter(region %in% "Atlantic", common_name %in% "Nurse shark", !is.na(Abundance))

mul<-5 ## multiplier to increase size of circle markers and size legend equally (keep constant across all maps)

m<-blacktip ## select which species to plot 

map<-
leaflet() %>%
  addProviderTiles(provider=providers$CartoDB.Positron) %>%
  addCircleMarkers(data = m,
                   lng= ~mean_lon, lat= ~mean_lat,
                   radius = ~mean_maxn * mul, fillColor = gg_color_hue(2)[m$Abundance],
                   color = "white",
                   weight = 2, fillOpacity=1,
                   label = as.character(round(m$mean_maxn,2))) %>%
    addLegend(colors = gg_color_hue(2), labels = levels(m$Abundance), opacity = 1) %>%
    addLegendCustom(colors = "grey",
                    # labels = pretty(seq(min(m$mean_maxn, na.rm=T), max(m$mean_maxn, na.rm=T), l=4)),
                    # sizes = pretty(seq(min(m$mean_maxn, na.rm=T), max(m$mean_maxn, na.rm=T), l=4)) * mul,
                    labels = seq(1,3, by=0.5),
                    sizes = seq(1,3, by=0.5) * mul,
                    shapes = "circle",
                    borders = "white") %>%
  addScaleBar("bottomleft")

map

mapshot(map, file="BlacktipBRUVS.png", remove_controls = c("zoomControl", "layersControl", "homeButton"))

# # Simple mapview
# m <- blacktip
# m %>%
#   st_as_sf(coords=c("mean_lon", "mean_lat"), crs=4326) %>%
#   mapview(zcol="Abundance", col.region=c(gg_color_hue(2)), #burst=T,
#           color=NA, lwd=1, legend=T, layer.name = "Abundance", alpha.region=1,
#           cex="mean_maxn") %>%
#   .@map %>%
#   addLegendCustom(colors = "grey",
#                   labels = pretty(seq(min(m$mean_maxn, na.rm=T), max(m$mean_maxn, na.rm=T), l=4)),
#                   sizes = pretty(seq(min(m$mean_maxn, na.rm=T), max(m$mean_maxn, na.rm=T), l=4)) * mul,
#                   shapes = "circle",
#                   borders = "grey")


