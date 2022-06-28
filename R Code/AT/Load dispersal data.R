## Load in the dispersal data files for $Daily, $Weekly and $Monthly subsets
Dispersal_Timescales_EagleRay <- readRDS(file = "Data/Dispersal_Timescales_EagleRay.RDS")
Dispersal_Timescales_Harlequin <- readRDS(file = "Data/Dispersal_Timescales_Harlequin.RDS")
Dispersal_Timescales_Kingfish <- readRDS(file = "Data/Dispersal_Timescales_Kingfish.RDS")
Dispersal_Timescales_Snapper <- readRDS(file = "Data/Dispersal_Timescales_Snapper.RDS")
Dispersal_Timescales_Trevally <- readRDS(file = "Data/Dispersal_Timescales_Trevally.RDS")
Dispersal_Timescales_BTWrasse<- readRDS(file = "Data/Dispersal_Timescales_BTWrasse.RDS")
Dispersal_Timescales_BlueGroper <- readRDS(file = "Data/Dispersal_Timescales_BlueGroper.RDS")
Dispersal_Timescales_BronzeWhaler <- readRDS(file = "Data/Dispersal_Timescales_BronzeWhaler.RDS")
Dispersal_Timescales_DuskyWhaler <- readRDS(file = "Data/Dispersal_Timescales_DuskyWhaler.RDS")



# Compute a histogram of distance per month
disp.hist <- Dispersal_Timescales_Trevally$Weekly %>%
  ggplot(aes(maxDistkm)) + geom_histogram() +theme_bw()
# Compute a histogram of num receiver stations
stat.hist <- Dispersal_Timescales_Trevally$Weekly %>%
  ggplot(aes(n_stations)) + geom_histogram() +theme_bw()
ggarrange(disp.hist,stat.hist)