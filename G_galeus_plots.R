
# De Wysiecki et al. 
# "Population-scale habitat use of school sharks (Triakidae: Galeorhinus galeus) in the Southwest Atlantic: 
# insights from temporally explicit niche modelling and habitat associations"

# Plots

library(viridis)
library(ggplot2)
library(SDMtune)
library(ellipsenm)
library(kuenm)
library(raster)
library(rgdal)
library(rgeos)
library(rgl)
library(gridExtra)
library(ggforce)
library(dplyr)

setwd('SET YOUR WORKING DIRECTORY')

# Coastline (spatial polygons) - taken from GSHHG coastline database of 15 June 2017
coastline <- 'DOWNLOAD AND READ'
coast <- readOGR(dsn = coastline, layer = 'GSHHS_f_L1_SouthAmerica')

# Projection
crs <- CRS('+init=epsg:4326')

# Species
Species <- 'Galeorhinus galeus'

# Season
season <- c('summer', 'autumn', 'winter', 'spring')

# Seed
set.seed(111)

# Function for 5% and 10% Minimum Training Presence calculation (taken from https://babichmorrowc.github.io/post/2019-04-12-sdm-threshold/)
sdm_threshold.5 <- function(sdm, occs, type = 'mtp', binary = FALSE){
  occPredVals <- raster::extract(sdm, occs)
  if(type == 'mtp'){
    thresh <- min(na.omit(occPredVals))
  } else if(type == 'p05'){
    if(length(occPredVals) < 10){
      p05 <- floor(length(occPredVals) * 0.95)
    } else {
      p05 <- ceiling(length(occPredVals) * 0.95)
    }
    thresh <- rev(sort(occPredVals))[p05]
  }
  sdm_thresh <- sdm
  sdm_thresh[sdm_thresh < thresh] <- NA
  if(binary){
    sdm_thresh[sdm_thresh >= thresh] <- 1
  }
  return(sdm_thresh)
}#MTP5
sdm_threshold.10 <- function(sdm, occs, type = 'mtp', binary = FALSE){
  occPredVals <- raster::extract(sdm, occs)
  if(type == 'mtp'){
    thresh <- min(na.omit(occPredVals))
  } else if(type == 'p10'){
    if(length(occPredVals) < 10){
      p10 <- floor(length(occPredVals) * 0.9)
    } else {
      p10 <- ceiling(length(occPredVals) * 0.9)
    }
    thresh <- rev(sort(occPredVals))[p10]
  }
  sdm_thresh <- sdm
  sdm_thresh[sdm_thresh < thresh] <- NA
  if(binary){
    sdm_thresh[sdm_thresh >= thresh] <- 1
  }
  return(sdm_thresh)
}#MTP10

# Functions for response plots
.get_presence <- function(swd) {
  return(swd@data[swd@pa == 1, , drop = FALSE])
}
.get_absence <- function(swd) {
  return(swd@data[swd@pa == 0, , drop = FALSE])
}

# Spatial polygons to erase
df_mar <- data.frame(xmin = -40, xmax = -30, ymin = -70, ymax = -50)
df_tierra <- data.frame(xmin = -75, xmax = -65, ymin = -25, ymax = -20)
df_rdp <- data.frame(y = c(-34.2257, -32.5165, -33.9072, -36.0317, -34.2257),
                     x = c(-60.2998, -58.4014, -55.8800, -58.1127, -60.2998))
df_lagoa <- data.frame(y = c(-29.7969, -30.2501, -31.1866, -31.7214, -31.9466, -32.1433, -32.2448, -29.7969),
                       x = c(-51.3938, -50.2811, -50.8794, -51.4996, -51.9246, -52.0806, -52.8042, -51.3938))

#------------------------------------------- Figure 1 ------------------------------------------

# Predictors relevant to the niche
var_names <- c('Depth', 'Surface_temperature', 'SST_fronts', 'Turbidity')

# Get calibration areas for each season

# Summer
env <- stack('summer.tif')[[2:5]]
names(env) <- var_names
dat_summer <- read.csv('2. Analyses/summer/G_galeus_joint.csv')
occ.tot <- SpatialPoints(cbind(dat_summer$longitude, dat_summer$latitude), proj4string = crs) 
occ.buff <- buffer(occ.tot, width = 1000000) 
env.M <- crop(env, occ.buff)
env.M <- mask(env.M, occ.buff)
env.M_summer <- stack(env.M)

# Autumn
env <- stack('autumn.tif')[[2:5]]
names(env) <- var_names
dat_autumn <- read.csv('2. Analyses/autumn/G_galeus_joint.csv')
occ.tot <- SpatialPoints(cbind(dat_autumn$longitude, dat_autumn$latitude), proj4string = crs) 
occ.buff <- buffer(occ.tot, width = 1000000) 
env.M <- crop(env, occ.buff)
env.M <- mask(env.M, occ.buff)
env.M_autumn <- stack(env.M)

# Winter
env <- stack('winter.tif')[[2:5]]
names(env) <- var_names
dat_winter <- read.csv('2. Analyses/winter/G_galeus_joint.csv')
occ.tot <- SpatialPoints(cbind(dat_winter$longitude, dat_winter$latitude), proj4string = crs) 
occ.buff <- buffer(occ.tot, width = 1000000) 
env.M <- crop(env, occ.buff)
env.M <- mask(env.M, occ.buff)
env.M_winter <- stack(env.M)

# Spring
env <- stack('spring.tif')[[2:5]]
names(env) <- var_names
dat_spring <- read.csv('2. Analyses/spring/G_galeus_joint.csv')
occ.tot <- SpatialPoints(cbind(dat_spring$longitude, dat_spring$latitude), proj4string = crs) 
occ.buff <- buffer(occ.tot, width = 1000000) 
env.M <- crop(env, occ.buff)
env.M <- mask(env.M, occ.buff)
env.M_spring <- stack(env.M)

# PCA to reduce environmental dimensions
pcas_summer <- kuenm_rpca(variables = env.M_summer, var.scale = T, write.result = F, n.pcs = 3)
pcas_autumn <- kuenm_rpca(variables = env.M_autumn, var.scale = T, write.result = F, n.pcs = 3)
pcas_winter <- kuenm_rpca(variables = env.M_winter, var.scale = T, write.result = F, n.pcs = 3)
pcas_spring <- kuenm_rpca(variables = env.M_spring, var.scale = T, write.result = F, n.pcs = 3)

# Preparing overlap objects to perform analyses
niche_summer <- overlap_object(dat_summer, species = 'species', longitude = 'longitude', latitude = 'latitude', method = 'mve1',
                               level = 95, variables = pcas_summer$PCRasters_initial)
niche_summer0 <- data.frame(season = 'summer', getValues(niche_summer@variables))
niche_summer0 <- sample_n(niche_summer0[which(complete.cases(niche_summer0)), ], 10000)
niche_autumn <- overlap_object(dat_autumn, species = 'species', longitude = 'longitude', latitude = 'latitude', method = 'mve1',
                               level = 95, variables = pcas_autumn$PCRasters_initial)
niche_autumn0 <- data.frame(season = 'autumn', getValues(niche_autumn@variables))
niche_autumn0 <- sample_n(niche_autumn0[which(complete.cases(niche_autumn0)), ], 10000)
niche_winter <- overlap_object(dat_winter, species = 'species', longitude = 'longitude', latitude = 'latitude', method = 'mve1',
                               level = 95, variables = pcas_winter$PCRasters_initial)
niche_winter0 <- data.frame(season = 'winter', getValues(niche_winter@variables))
niche_winter0 <- sample_n(niche_winter0[which(complete.cases(niche_winter0)), ], 10000)
niche_spring <- overlap_object(dat_spring, species = 'species', longitude = 'longitude', latitude = 'latitude', method = 'mve1',
                               level = 95, variables = pcas_spring$PCRasters_initial)
niche_spring0 <- data.frame(season = 'spring', getValues(niche_spring@variables))
niche_spring0 <- sample_n(niche_spring0[which(complete.cases(niche_spring0)), ], 10000)

final_back_pcs0 <- rbind(niche_summer0, niche_autumn0, niche_winter0, niche_spring0)
final_back_pcs0 <- transform(final_back_pcs0, season = factor(season, levels = c('spring', 'autumn', 'winter', 'summer'))) # ordering of factors

# Overlap
N.test1 <- ellipsoid_overlap(niche_summer, niche_autumn, overlap_type = 'back_union', significance_test = F)
summer_occ_pcs <- data.frame(season = 'summer', N.test1@data[[1]])
N.test2 <- ellipsoid_overlap(niche_autumn, niche_winter, overlap_type = 'back_union', significance_test = F)
autumn_occ_pcs <- data.frame(season = 'autumn', N.test2@data[[1]])
N.test3 <- ellipsoid_overlap(niche_winter, niche_spring, overlap_type = 'back_union', significance_test = F)
winter_occ_pcs <- data.frame(season = 'winter', N.test3@data[[1]])
N.test4 <- ellipsoid_overlap(niche_spring, niche_summer, overlap_type = 'back_union', significance_test = F)
spring_occ_pcs <- data.frame(season = 'spring', N.test4@data[[1]])

final_occ_pcs0 <- rbind(summer_occ_pcs, autumn_occ_pcs, winter_occ_pcs, spring_occ_pcs)
final_occ_pcs0 <- transform(final_occ_pcs0, season = factor(season, levels = c('spring', 'autumn', 'winter', 'summer'))) # ordering of factors

Labs <- as_labeller(c(summer_autumn = 'summer~vs~autumn', autumn_winter = 'autumn~vs~winter',
                    winter_spring = 'winter~vs~spring', spring_summer = 'spring~vs~summer'), default = label_parsed)
  
# Plot
ggplot() + 
  geom_vline(xintercept = 0, linetype = 'dashed', size = 0.2) + geom_hline(yintercept = 0, linetype = 'dashed', size = 0.2) +
  geom_point(data = final_back_pcs0, aes(x = pc_1, y = pc_2), color = 'grey10', size = 2.5, shape = 1) +
  geom_point(data = final_occ_pcs0, aes(x = pc_1, y = pc_2, color = season), size = 1) +
  geom_mark_hull(data = final_occ_pcs0, aes(x = pc_1, y = pc_2, color = season, fill = season), expand = 0.02, size = 0.9) +
  scale_fill_viridis_d(option = 'E', alpha = 0.3) + scale_color_viridis_d(option = 'E') +
  xlim(-8, 2) + ylim(-10, 10) + xlab('PC 1') + ylab('PC 2') + 
  theme(panel.grid = element_blank()) #, legend.position = 'none'
ggsave('Figure 1a.pdf', width = 10, height = 6, units = 'cm')

ggplot() + 
  geom_vline(xintercept = 0, linetype = 'dashed', size = 0.2) + geom_hline(yintercept = 0, linetype = 'dashed', size = 0.2) +
  geom_point(data = final_back_pcs0, aes(x = pc_2, y = pc_3), color = 'grey10', size = 2.5, shape = 1) +
  geom_point(data = final_occ_pcs0, aes(x = pc_2, y = pc_3, color = season), size = 1) +
  geom_mark_hull(data = final_occ_pcs0, aes(x = pc_2, y = pc_3, color = season, fill = season), expand = 0.02, size = 0.9) +
  scale_fill_viridis_d(option = 'E', alpha = 0.3) + scale_color_viridis_d(option = 'E') +
  xlim(-10, 5) + ylim(-10, 5) + xlab('PC 2') + ylab('PC 3') + 
  theme(panel.grid = element_blank()) #, legend.position = 'none'
ggsave('Figure 3b.pdf', width = 10, height = 6, units = 'cm')


#------------------------------------------- Figure 2 ------------------------------------------

# Occurrences by sex

# Annual
final_mod0 <- raster('2. Analyses/annual/Final_Model_Stats/Statistics_E/calibration_median.tif') # model
final_mod0 <- crop(final_mod0, extent(final_mod0) + c(0.5, 0, 0, 0))
final_mod <- extend(final_mod0, extent(final_mod0) + 1.75)
final_mod[final_mod >= 0] <- 1

# records
dat <- read.csv('Cazon_corrected.csv') # all occurrences
data_f <- subset(dat, Sex == 'female')
data_m <- subset(dat, Sex == 'male')
data_x <- subset(dat, Sex == '')

df <- data.frame(coordinates(final_mod), as.data.frame(final_mod))

ggplot() +
  geom_polygon(data = coast, aes(x = long, y = lat, group = group), color = 'black', fill = 'grey80', size = 0.15) +
  geom_point(data = data_x, aes(x = Longitude, y = Latitude), shape = 21, size = 1.5, stroke = 0.01, color = 'black', fill = 'grey95') +
  geom_point(data = data_f, aes(x = Longitude, y = Latitude), shape = 21, size = 1.5, stroke = 0.1, color = 'black', fill = '#FFEA46FF') +
  geom_point(data = data_m, aes(x = Longitude, y = Latitude), shape = 4, size = 0.3, stroke = 0.4, color = '#00204DFF') +  
  scale_y_continuous(name = NULL, breaks = c(-50, -40, -30), labels = c('50º', '40º', '30º')) + 
  scale_x_continuous(name = NULL, breaks = c(-65, -60, -55, -50), labels = c('65º', '60º', '55º', '50º')) +
  ggsn::scalebar(x.min = min(df$x), x.max = max(df$x), y.min = min(df$y), y.max = max(df$y), transform = T,
                 dist = 200, st.size = 2, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.25, anchor = c(x = -47, y = -49.35)) +
  coord_equal(xlim = c(-69.4, -45.5), ylim = c(-50.85, -26), expand = 0) + 
  theme(panel.background = element_rect(fill = NA),
        panel.grid = element_blank(), legend.position = 'none',
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.5)) 
ggsave('Figure 2a_annual.pdf', width = 6.7, height = 6.5, units = 'cm')

# Seasonal
final_mod00 <- raster('2. Analyses/annual/Final_Model_Stats/Statistics_E/calibration_median.tif') # model
dat0 <- read.csv('Cazon_corrected.csv') # all occurrences

for(i in 1:4){
  final_mod0 <- crop(final_mod00, extent(final_mod00) + c(0.5, 0, 0, 0))
  final_mod <- extend(final_mod0, extent(final_mod0) + 1.75)
  final_mod[final_mod >= 0] <- 1
  
  # records
  dat <- subset(dat0, Season == season[i])
  data_f <- subset(dat, Sex == 'female')
  data_m <- subset(dat, Sex == 'male')
  data_x <- subset(dat, Sex == '')
  
  # Calibration areas
  occ.tot <- SpatialPoints(cbind(dat$Longitude, dat$Latitude), proj4string = crs)
  M.buff <- buffer(occ.tot, width = 100000)
  M.buff <- crop(M.buff, extent(final_mod0))
  
  df <- data.frame(coordinates(final_mod), as.data.frame(final_mod))
  
  ggplot() +
    geom_polygon(data = coast, aes(x = long, y = lat, group = group), color = 'black', fill = 'grey80', size = 0.15) +
    geom_point(data = data_x, aes(x = Longitude, y = Latitude), shape = 21, size = 1.5, stroke = 0.01, color = 'black', fill = 'grey95') +
    geom_point(data = data_f, aes(x = Longitude, y = Latitude), shape = 21, size = 1.5, stroke = 0.1, color = 'black', fill = '#FFEA46FF') +
    geom_point(data = data_m, aes(x = Longitude, y = Latitude), shape = 4, size = 0.3, stroke = 0.4, color = '#00204DFF') +  
    scale_y_continuous(name = NULL, breaks = c(-50, -40, -30), labels = c('50º', '40º', '30º')) + 
    scale_x_continuous(name = NULL, breaks = c(-65, -60, -55, -50), labels = c('65º', '60º', '55º', '50º')) +
    coord_equal(xlim = c(-69.4, -45.5), ylim = c(-50.85, -26), expand = 0) + 
    theme(panel.background = element_rect(fill = NA),
          panel.grid = element_blank(), legend.position = 'none',
          panel.border = element_rect(colour = 'black', fill = NA, size = 0.5)) 
  ggsave(paste('Figure 2a_', season[i], '.pdf', sep = ''), width = 6.7, height = 6.5, units = 'cm')
}


# Occurrences by stage

# Annual
final_mod0 <- raster('2. Analyses/annual/Final_Model_Stats/Statistics_E/calibration_median.tif') # model
final_mod0 <- crop(final_mod0, extent(final_mod0) + c(0.5, 0, 0, 0))
final_mod <- extend(final_mod0, extent(final_mod0) + 1.75)
final_mod[final_mod >= 0] <- 1

# records
dat <- read.csv('Cazon_corrected.csv') # all occurrences
data_a <- subset(dat, Stage == 'adult')
data_j <- subset(dat, Stage == 'juvenile')
data_x <- subset(dat, Stage == '')

# Calibration areas
occ.tot <- SpatialPoints(cbind(dat$Longitude, dat$Latitude), proj4string = crs)
M.buff <- buffer(occ.tot, width = 100000)
M.buff <- crop(M.buff, extent(final_mod0))

df <- data.frame(coordinates(final_mod), as.data.frame(final_mod))

ggplot() +
  geom_polygon(data = coast, aes(x = long, y = lat, group = group), color = 'black', fill = 'grey80', size = 0.15) +
  geom_point(data = data_x, aes(x = Longitude, y = Latitude), shape = 21, size = 1.5, stroke = 0.01, color = 'black', fill = 'grey95') +
  geom_point(data = data_j, aes(x = Longitude, y = Latitude), shape = 21, size = 1.5, stroke = 0.1, color = 'black', fill = '#FFEA46FF') +
  geom_point(data = data_a, aes(x = Longitude, y = Latitude), shape = 4, size = 0.3, stroke = 0.4, color = '#00204DFF') +  
  scale_y_continuous(name = NULL, breaks = c(-50, -40, -30), labels = c('50º', '40º', '30º')) + 
  scale_x_continuous(name = NULL, breaks = c(-65, -60, -55, -50), labels = c('65º', '60º', '55º', '50º')) +
  ggsn::scalebar(x.min = min(df$x), x.max = max(df$x), y.min = min(df$y), y.max = max(df$y), transform = T,
                 dist = 200, st.size = 2, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.25, anchor = c(x = -47, y = -49.35)) +
  coord_equal(xlim = c(-69.4, -45.5), ylim = c(-50.85, -26), expand = 0) + 
  theme(panel.background = element_rect(fill = NA),
        panel.grid = element_blank(), legend.position = 'none',
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.5)) 
ggsave('Figure 4b_annual.pdf', width = 6.7, height = 6.5, units = 'cm')

# Seasonal
final_mod00 <- raster('2. Analyses/annual/Final_Model_Stats/Statistics_E/calibration_median.tif') # model
dat0 <- read.csv('Cazon_corrected.csv') # all occurrences

for(i in 1:4){
  final_mod0 <- crop(final_mod00, extent(final_mod00) + c(0.5, 0, 0, 0))
  final_mod <- extend(final_mod0, extent(final_mod0) + 1.75)
  final_mod[final_mod >= 0] <- 1
  
  # records
  dat <- subset(dat0, Season == season[i])
  data_a <- subset(dat, Stage == 'adult')
  data_j <- subset(dat, Stage == 'juvenile')
  data_x <- subset(dat, Stage == '')
  
  # Calibration areas
  occ.tot <- SpatialPoints(cbind(dat$Longitude, dat$Latitude), proj4string = crs)
  M.buff <- buffer(occ.tot, width = 100000)
  M.buff <- crop(M.buff, extent(final_mod0))
  
  df <- data.frame(coordinates(final_mod), as.data.frame(final_mod))
  
  ggplot() +
    geom_polygon(data = coast, aes(x = long, y = lat, group = group), color = 'black', fill = 'grey80', size = 0.15) +
    geom_point(data = data_x, aes(x = Longitude, y = Latitude), shape = 21, size = 1.5, stroke = 0.01, color = 'black', fill = 'grey95') +
    geom_point(data = data_j, aes(x = Longitude, y = Latitude), shape = 21, size = 1.5, stroke = 0.1, color = 'black', fill = '#FFEA46FF') +
    geom_point(data = data_a, aes(x = Longitude, y = Latitude), shape = 4, size = 0.3, stroke = 0.4, color = '#00204DFF') +  
    scale_y_continuous(name = NULL, breaks = c(-50, -40, -30), labels = c('50º', '40º', '30º')) + 
    scale_x_continuous(name = NULL, breaks = c(-65, -60, -55, -50), labels = c('65º', '60º', '55º', '50º')) +
    coord_equal(xlim = c(-69.4, -45.5), ylim = c(-50.85, -26), expand = 0) + 
    theme(panel.background = element_rect(fill = NA),
          panel.grid = element_blank(), legend.position = 'none',
          panel.border = element_rect(colour = 'black', fill = NA, size = 0.5)) 
  ggsave(paste('Figure 4b_', season[i], '.pdf', sep = ''), width = 6.7, height = 6.5, units = 'cm')
}


#------------------------------------------- Figure 3 ------------------------------------------

# Go to script 'G_galeus_habitat_associations.R'


#------------------------------------------- Figure 4 -------------------------------------------

# Core area of distribution

# Annual calibration areas
final_mod000 <- raster('2. Analyses/annual/Final_Model_Stats/Statistics_E/calibration_median.tif') # model
final_mod000 <- crop(final_mod000, extent(final_mod000) + c(0.5, 0, 0, 0))
final_mod00 <- extend(final_mod000, extent(final_mod000) + 1.75)
occ_cal <- read.csv('2. Analyses/annual/G_galeus_joint.csv') # calibration occurrences
occ.tot <- SpatialPoints(cbind(occ_cal$longitude, occ_cal$latitude), proj4string = crs)
M.buff <- buffer(occ.tot, width = 1000000)
M.buff <- crop(M.buff, extent(final_mod000))
bg_binary <- final_mod00
bg_binary[bg_binary < 0.1] <- 1
bg_binary[bg_binary] <- 2

# Summer
final_mod0 <- raster('2. Analyses/summer/Final_Model_Stats/Statistics_E/calibration_median.tif') # model
final_mod0 <- crop(final_mod0, extent(final_mod0) + c(0.5, 0, 0, 0))
final_mod <- extend(final_mod0, extent(final_mod00) + 1.75)
occ_cal <- read.csv('2. Analyses/summer/G_galeus_joint.csv') # calibration occurrences
# Thresholds (%5 minimumm training presence)
mtp5 <- sdm_threshold.5(final_mod, occ_cal[, c('longitude', 'latitude')], 'p05', binary = F)
bin_5_summer <- final_mod
bin_5_summer[bin_5_summer < minValue(mtp5)] <- NA
bin_5_summer[bin_5_summer >= minValue(mtp5)] <- 1
# Thresholds (%10 minimumm training presence)
mtp10 <- sdm_threshold.10(final_mod, occ_cal[, c('longitude', 'latitude')], 'p10', binary = F)
bin_10_summer <- final_mod
bin_10_summer[bin_10_summer < minValue(mtp10)] <- NA
bin_10_summer[bin_10_summer >= minValue(mtp10)] <- 1

# Autumn
final_mod0 <- raster('2. Analyses/autumn/Final_Model_Stats/Statistics_E/calibration_median.tif') # model
final_mod0 <- crop(final_mod0, extent(final_mod0) + c(0.5, 0, 0, 0))
final_mod <- extend(final_mod0, extent(final_mod00) + 1.75)
occ_cal <- read.csv('2. Analyses/autumn/G_galeus_joint.csv') # calibration occurrences
# Thresholds (%5 minimumm training presence)
mtp5 <- sdm_threshold.5(final_mod, occ_cal[, c('longitude', 'latitude')], 'p05', binary = F)
bin_5_autumn <- final_mod
bin_5_autumn[bin_5_autumn < minValue(mtp5)] <- NA
bin_5_autumn[bin_5_autumn >= minValue(mtp5)] <- 1
# Thresholds (%10 minimumm training presence)
mtp10 <- sdm_threshold.10(final_mod, occ_cal[, c('longitude', 'latitude')], 'p10', binary = F)
bin_10_autumn <- final_mod
bin_10_autumn[bin_10_autumn < minValue(mtp10)] <- NA
bin_10_autumn[bin_10_autumn >= minValue(mtp10)] <- 1

# Winter
final_mod0 <- raster('2. Analyses/winter/Final_Model_Stats/Statistics_E/calibration_median.tif') # model
final_mod0 <- crop(final_mod0, extent(final_mod0) + c(0.5, 0, 0, 0))
final_mod <- extend(final_mod0, extent(final_mod00) + 1.75)
occ_cal <- read.csv('2. Analyses/winter/G_galeus_joint.csv') # calibration occurrences
# Thresholds (%5 minimumm training presence)
mtp5 <- sdm_threshold.5(final_mod, occ_cal[, c('longitude', 'latitude')], 'p05', binary = F)
bin_5_winter <- final_mod
bin_5_winter[bin_5_winter < minValue(mtp5)] <- NA
bin_5_winter[bin_5_winter >= minValue(mtp5)] <- 1
# Thresholds (%10 minimumm training presence)
mtp10 <- sdm_threshold.10(final_mod, occ_cal[, c('longitude', 'latitude')], 'p10', binary = F)
bin_10_winter <- final_mod
bin_10_winter[bin_10_winter < minValue(mtp10)] <- NA
bin_10_winter[bin_10_winter >= minValue(mtp10)] <- 1

# Spring
final_mod0 <- raster('2. Analyses/spring/Final_Model_Stats/Statistics_E/calibration_median.tif') # model
final_mod0 <- crop(final_mod0, extent(final_mod0) + c(0.5, 0, 0, 0))
final_mod <- extend(final_mod0, extent(final_mod00) + 1.75)
occ_cal <- read.csv('2. Analyses/spring/G_galeus_joint.csv') # calibration occurrences
# Thresholds (%5 minimumm training presence)
mtp5 <- sdm_threshold.5(final_mod, occ_cal[, c('longitude', 'latitude')], 'p05', binary = F)
bin_5_spring <- final_mod
bin_5_spring[bin_5_spring < minValue(mtp5)] <- NA
bin_5_spring[bin_5_spring >= minValue(mtp5)] <- 1
# Thresholds (%10 minimumm training presence)
mtp10 <- sdm_threshold.10(final_mod, occ_cal[, c('longitude', 'latitude')], 'p10', binary = F)
bin_10_spring <- final_mod
bin_10_spring[bin_10_spring < minValue(mtp10)] <- NA
bin_10_spring[bin_10_spring >= minValue(mtp10)] <- 1

Col <- c('#00204DFF', '#414D6BFF', '#7C7B78FF', '#BCAF6FFF', '#FFEA46FF')

# Mosaic 10% MTP
mod_mosaic_10 <- mosaic(bg_binary, bin_10_summer, bin_10_autumn, bin_10_winter, bin_10_spring, fun = sum)
mod_mosaic_10[mod_mosaic_10 < 2] <- NA
df_10 <- data.frame(coordinates(mod_mosaic_10), as.data.frame(mod_mosaic_10))

ggplot() +
  geom_tile(data = df_10, aes(x = x, y = y, fill = as.factor(layer))) + 
  scale_fill_manual(values = Col, na.value = 'white') + 
  geom_polygon(data = M.buff, aes(x = long, y = lat, group = group), size = 0.5, color = 'black', fill = NA) +
  geom_polygon(data = df_rdp, aes(x = x, y = y), color = 'black', fill = 'white', size = 0.5) +
  geom_polygon(data = df_lagoa, aes(x = x, y = y), color = 'black', fill = 'white',size = 0.1) +
  geom_polygon(data = coast, aes(x = long, y = lat, group = group), color = 'black', fill = 'grey80', size = 0.15) +
  geom_rect(data = df_mar, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'white') +
  geom_rect(data = df_tierra, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'grey80') +
  scale_y_continuous(name = NULL, breaks = c(-55, -45, -35, -25), labels = c('55º', '45º', '35º', '25º')) + 
  scale_x_continuous(name = NULL, breaks = c(-70, -60, -50, -40), labels = c('70º', '60º', '50º', '40º')) +
  ggsn::scalebar(x.min = min(df_10$x), x.max = max(df_10$x), y.min = min(df_10$y), y.max = max(df_10$y), transform = T,
                 dist = 200, st.size = 2, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.25, anchor = c(x = -39.5, y = -57.35)) +
  coord_equal(xlim = c(-70.2, -37.1), ylim = c(-59.05, -21.5), expand = 0) + 
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(), legend.position = 'none',
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.5)) 
ggsave('Figure 4a.tiff', dpi = 900, width = 6.2, height = 6.6, units = 'cm')

# Mosaic 5% MTP
mod_mosaic_5 <- mosaic(bg_binary, bin_5_summer, bin_5_autumn, bin_5_winter, bin_5_spring, fun = sum)
mod_mosaic_5[mod_mosaic_5 < 2] <- NA
df_5 <- data.frame(coordinates(mod_mosaic_5), as.data.frame(mod_mosaic_5))

ggplot() +
  geom_tile(data = df_5, aes(x = x, y = y, fill = as.factor(layer))) + 
  scale_fill_manual(values = Col, na.value = 'white') + 
  geom_polygon(data = M.buff, aes(x = long, y = lat, group = group), size = 0.5, color = 'black', fill = NA) +
  geom_polygon(data = df_rdp, aes(x = x, y = y), color = 'black', fill = 'white', size = 0.5) +
  geom_polygon(data = df_lagoa, aes(x = x, y = y), color = 'black', fill = 'white',size = 0.1) +
  geom_polygon(data = coast, aes(x = long, y = lat, group = group), color = 'black', fill = 'grey80', size = 0.15) +
  geom_rect(data = df_mar, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'white') +
  geom_rect(data = df_tierra, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'grey80') +
  scale_y_continuous(name = NULL, breaks = c(-55, -45, -35, -25), labels = c('55º', '45º', '35º', '25º')) + 
  scale_x_continuous(name = NULL, breaks = c(-70, -60, -50, -40), labels = c('70º', '60º', '50º', '40º')) +
  coord_equal(xlim = c(-70.2, -37.1), ylim = c(-59.05, -21.5), expand = 0) + 
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(), legend.position = 'none',
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.5)) 
ggsave('Figure 4b.tiff', dpi = 900, width = 6.2, height = 6.6, units = 'cm')


#------------------------------------------- Figure S1 --------------------------------------------

# Annual continuous habitat suitability

final_mod0 <- raster('2. Analyses/annual/Final_Model_Stats/Statistics_E/calibration_median.tif') # model
final_mod0 <- crop(final_mod0, extent(final_mod0) + c(0.5, 0, 0, 0))
final_mod <- extend(final_mod0, extent(final_mod0) + 1.75)
occ_cal <- read.csv('2. Analyses/annual/G_galeus_joint.csv') # calibration occurrences
occ_ind <- read.csv('2. Analyses/annual/G_galeus_indep.csv') # independent occurrences

# Calibration areas
occ.tot <- SpatialPoints(cbind(occ_cal$longitude, occ_cal$latitude), proj4string = crs)
M.buff <- buffer(occ.tot, width = 1000000)
M.buff <- crop(M.buff, extent(final_mod0))

df <- data.frame(coordinates(final_mod), as.data.frame(final_mod))

# Plot
ggplot() +
  geom_tile(data = df, aes(x = x, y = y, fill = calibration_median)) + 
  scale_fill_viridis_c(option = 'E', na.value = 'white') + 
  geom_polygon(data = M.buff, aes(x = long, y = lat, group = group), size = 0.5, color = 'black', fill = NA) +
  geom_polygon(data = df_rdp, aes(x = x, y = y), color = 'black', fill = 'white', size = 0.5) +
  geom_polygon(data = df_lagoa, aes(x = x, y = y), color = 'black', fill = 'white',size = 0.1) +
  geom_polygon(data = coast, aes(x = long, y = lat, group = group), color = 'black', fill = 'grey80', size = 0.15) +
  geom_rect(data = df_mar, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'white') +
  geom_rect(data = df_tierra, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'grey80') +
  geom_point(data = occ_cal, aes(x = longitude, y = latitude), size  =  0.05, color = 'grey5', fill = 'grey5') +
  geom_point(data = occ_ind, aes(x = longitude, y = latitude), size  =  0.05, color = 'white', fill = 'white') +
  scale_y_continuous(name = NULL, breaks = c(-55, -45, -35, -25), labels = c('55º', '45º', '35º', '25º')) + 
  scale_x_continuous(name = NULL, breaks = c(-70, -60, -50, -40), labels = c('70º', '60º', '50º', '40º')) +
  ggsn::scalebar(x.min = min(df$x), x.max = max(df$x), y.min = min(df$y), y.max = max(df$y), transform = T, 
                 dist = 200, st.size = 2, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.25, anchor = c(x = -47, y = -57.4)) +
  coord_equal(xlim = c(-70.2, -37.1), ylim = c(-59.05, -21.5), expand = 0) + 
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(), legend.position = 'none',
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.5)) 
ggsave('Figure S1a.tiff', dpi = 900, width = 6.2, height = 6.6, units = 'cm')


# Annual binary habitat suitability

bg_binary <- final_mod
bg_binary[bg_binary < 0.1] <- 1
bg_binary[bg_binary] <- 2

# Thresholds (%5 minimumm training presence)
mtp5 <- sdm_threshold.5(final_mod, occ_cal[, c('longitude', 'latitude')], 'p05', binary = F)
bin_5 <- final_mod
bin_5[bin_5 < minValue(mtp5)] <- NA
bin_5[bin_5 >= minValue(mtp5)] <- 3

# Thresholds (%10 minimumm training presence)
mtp10 <- sdm_threshold.10(final_mod, occ_cal[, c('longitude', 'latitude')], 'p10', binary = F)
bin_10 <- final_mod
bin_10[bin_10 < minValue(mtp10)] <- NA
bin_10[bin_10 >= minValue(mtp10)] <- 4

# Mosaic
mod_mosaic <- mosaic(bg_binary, bin_5, bin_10, fun = sum)
mod_mosaic[mod_mosaic < 1] <- NA

Col <- c('#00204DFF', '#7C7B78FF', '#FFEA46FF')
df <- data.frame(coordinates(mod_mosaic), as.data.frame(mod_mosaic))

# Plot
ggplot() +
  geom_tile(data = df, aes(x = x, y = y, fill = as.factor(layer))) + 
  scale_fill_manual(values = Col, na.value = 'white') + 
  geom_polygon(data = M.buff, aes(x = long, y = lat, group = group), size = 0.5, color = 'black', fill = NA) +
  geom_polygon(data = df_rdp, aes(x = x, y = y), color = 'black', fill = 'white', size = 0.5) +
  geom_polygon(data = df_lagoa, aes(x = x, y = y), color = 'black', fill = 'white',size = 0.1) +
  geom_polygon(data = coast, aes(x = long, y = lat, group = group), color = 'black', fill = 'grey80', size = 0.15) +
  geom_rect(data = df_mar, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'white') +
  geom_rect(data = df_tierra, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'grey80') +
  scale_y_continuous(name = NULL, breaks = c(-55, -45, -35, -25), labels = c('55º', '45º', '35º', '25º')) + 
  scale_x_continuous(name = NULL, breaks = c(-70, -60, -50, -40), labels = c('70º', '60º', '50º', '40º')) + 
  coord_equal(xlim = c(-70.2, -37.1), ylim = c(-59.05, -21.5), expand = 0) + 
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(), legend.position = 'none',
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.5)) 
ggsave('Figure S1b.tiff', dpi = 900, width = 6.2, height = 6.6, units = 'cm')


#------------------------------------------- Figure S2 ------------------------------------------

# Response plots for the annual model

# Read predictors
env <- stack('env.tif')
var_names <- c('Distance_to_coast', 'Depth', 'Surface_temperature', 'SST_fronts', 'Turbidity', 
               'Primary_productivity')
names(env) <- var_names

# Subset predictors relevant to the final model
var_names <- c('Depth', 'Surface_temperature', 'SST_fronts', 'Turbidity')
env <- env[[var_names]]
env[['Depth']] <- env[['Depth']] * -1

p_df <- data.frame()
a_df <- data.frame()
plot_data_df <- data.frame()

# Occurrences
dat <- read.csv('2. Analyses/annual/G_galeus_joint.csv')
dat <- dat[, c('longitude', 'latitude')]

# Rio de la Plata and Lagoa dos Patos spatial polygons
df_rdp <- matrix(c(-60.2998, -58.4014, -55.8800, -58.1127, -60.2998, 
                   -34.2257, -32.5165, -33.9072, -36.0317, -34.2257), 5, 2)
df_rdp <- SpatialPolygons(list(Polygons(list(Polygon(df_rdp)), ID = 1)), proj4string = crs)
df_lagoa <- matrix(c(-51.3938, -50.2811, -50.8794, -51.4996, -51.9246, -52.0806, -52.8042, -51.3938,
                     -29.7969, -30.2501, -31.1866, -31.7214, -31.9466, -32.1433, -32.2448, -29.7969), 8, 2)
df_lagoa <- SpatialPolygons(list(Polygons(list(Polygon(df_lagoa)), ID = 1)), proj4string = crs)

for(j in var_names){
  
  # prepare data
  Env <- env[[j]]

  # calibration areas
  occ.tot <- SpatialPoints(dat, proj4string = crs) 
  occ.buff <- buffer(occ.tot, width = 1000000) 
  env.M <- crop(Env, extent(occ.buff))
  env.M <- mask(env.M, occ.buff) 
  env.M <- mask(env.M, df_rdp, inverse = T)
  env.M <- mask(env.M, df_lagoa, inverse = T)
  env.M <- stack(env.M)
  
  # background points
  set.seed(111)
  notna <- which(complete.cases(values(env.M)))
  samp <- sample(notna, 10000, replace = F)
  samplocs <- as.data.frame(xyFromCell(env.M, samp))
  
  # SWD object
  data <- prepareSWD(species = Species, p = dat, a = samplocs, env = env.M)
  
  # run maxent replicates
  folds = randomFolds(data, k = 10, only_presence = T, seed = 111)
  default_model <- train(method = 'Maxent', data = data, fc = 'lq', reg = 0.1, iter = 1000, folds = folds)
  
  # presences and absences with variable data
  p <- .get_presence(default_model@data)
  p$var <- j
  names(p)[names(p) == j] <- 'values'
  a <- .get_absence(default_model@data)
  a$var <- j
  names(a)[names(a) == j] <- 'values'

  pred <- as.data.frame(matrix(data = NA, nrow = dim(data@data)[1], ncol = 10))
  for(i in 1:10){
    pred[, i] <- predict(default_model@models[[i]], data = data@data, type = 'logistic')
  }
  
  # plot data
  plot_data <- as.data.frame(matrix(data = NA, nrow = dim(data@data)[1], ncol = 4))
  names(plot_data) <- c('mean', 'sd', 'max', 'min')
  plot_data$mean <- rowMeans(pred)
  plot_data$sd <- apply(pred, 1, sd)
  plot_data$max <- plot_data$mean + plot_data$sd
  plot_data$min <- plot_data$mean - plot_data$sd
  plot_data$var <- j
  plot_data$values <- data@data[, j]
  
  # data frames
  p_df <- rbind(p_df, p)
  a_df <- rbind(a_df, a)
  plot_data_df <- rbind(plot_data_df, plot_data)
}

# Plot preparation 
plot_data_df <- plot_data_df[!(plot_data_df$var == 'Depth' & plot_data_df$values > 2000), ]
plot_data_df <- plot_data_df[!(plot_data_df$var == 'Surface_temperature' & plot_data_df$values > 30), ]
plot_data_df <- plot_data_df[!(plot_data_df$var == 'Surface_temperature' & plot_data_df$values < 5), ]
plot_data_df <- plot_data_df[!(plot_data_df$var == 'SST_fronts' & plot_data_df$values > 1), ]
plot_data_df <- plot_data_df[!(plot_data_df$var == 'Turbidity' & plot_data_df$values > 1.2), ]
plot_data_df <- transform(plot_data_df, var = factor(var, levels = c('Depth', 'Surface_temperature', 'SST_fronts', 'Turbidity')))

a_df <- a_df[!(a_df$var == 'Depth' & a_df$values > 2000), ]
a_df <- a_df[!(a_df$var == 'Surface_temperature' & a_df$values > 30), ]
a_df <- a_df[!(a_df$var == 'Surface_temperature' & a_df$values < 5), ]
a_df <- a_df[!(a_df$var == 'SST_fronts' & a_df$values > 1), ]
a_df <- a_df[!(a_df$var == 'Turbidity' & a_df$values > 1.2), ]
a_df <- transform(a_df, var = factor(var, levels = c('Depth', 'Surface_temperature', 'SST_fronts', 'Turbidity')))

p_df <- p_df[!(p_df$var == 'Depth' & p_df$values > 2000), ]
p_df <- p_df[!(p_df$var == 'Surface_temperature' & p_df$values > 30), ]
p_df <- p_df[!(p_df$var == 'Surface_temperature' & p_df$values < 5), ]
p_df <- p_df[!(p_df$var == 'SST_fronts' & p_df$values > 1), ]
p_df <- p_df[!(p_df$var == 'Turbidity' & p_df$values > 1.2), ]
p_df <- transform(p_df, var = factor(var, levels = c('Depth', 'Surface_temperature', 'SST_fronts', 'Turbidity')))

# Labeller
var_names <- as_labeller(c(Depth = 'Bathymetry~(m)', Turbidity = 'Coefficient~Kd490~(m^-1)',
                           Surface_temperature = 'Surface~Temperature~(ºC)', SST_fronts = 'Surface~Temperature~fronts~(ºC)'), default = label_parsed)

# Plot
ggplot(data = plot_data_df, aes(x = values, y = mean, ymin = min, ymax = max)) + 
  geom_line(colour = '#00204DFF') + 
  geom_ribbon(fill = '#00204DFF', alpha = 0.2) +
  geom_rug(data = p_df, inherit.aes = F, aes(values), sides = 't', color = '#FFEA46FF', size = 0.3) + 
  geom_rug(data = a_df, inherit.aes = F, aes(values), sides = 'b', color = '#7C7B78FF', size = 0.3) + 
  labs(x = NULL, y = 'Logistic output') + ylim(0, 1) +
  theme(panel.background = element_blank(), 
        panel.grid.minor = element_blank(),                                       
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.35),
        panel.grid.major = element_line(size = 0.2, colour = 'grey90'),
        strip.background = element_rect(fill = 'transparent'),
        strip.text = element_text(vjust = -0.5, size = 8),
        panel.spacing.y = unit(-0.3, 'lines'),
        panel.spacing.x = unit(0.7, 'lines'),
        plot.margin = unit(c(-0.2, 0.4, 0.2, 0.2), 'cm'),
        axis.title = element_text(size = 10)) +
  facet_wrap(~ var, scales = 'free_x', labeller = var_names, ncol = 4) 
ggsave('Figure S2.pdf', width = 19, height = 4.6, units = 'cm')


#------------------------------------------- Figure S3 ------------------------------------------

# Seasonal binary habitat suitability

final_mod00 <- raster('2. Analyses/annual/Final_Model_Stats/Statistics_E/calibration_median.tif') # model

for(i in 1:4){
  
  final_mod0 <- raster(paste('2. Analyses/', season[i], '/Final_Model_Stats/Statistics_E/calibration_median.tif', sep = '')) # model
  final_mod0 <- crop(final_mod0, extent(final_mod0) + c(0.5, 0, 0, 0))
  final_mod <- extend(final_mod0, extent(final_mod00) + 1.75)
  occ_cal <- read.csv(paste('2. Analyses/', season[i], '/G_galeus_joint.csv', sep = '')) # calibration occurrences
  
  # Calibration areas
  occ.tot <- SpatialPoints(cbind(occ_cal$longitude, occ_cal$latitude), proj4string = crs)
  M.buff <- buffer(occ.tot, width = 1000000)
  M.buff <- crop(M.buff, extent(final_mod0))
  
  # Binary
  bg_binary <- final_mod
  bg_binary[bg_binary < 0.1] <- 1
  bg_binary[bg_binary] <- 2
  
  # Thresholds (%5 minimumm training presence)
  mtp5 <- sdm_threshold.5(final_mod, occ_cal[, c('longitude', 'latitude')], 'p05', binary = F)
  bin_5 <- final_mod
  bin_5[bin_5 < minValue(mtp5)] <- NA
  bin_5[bin_5 >= minValue(mtp5)] <- 3
  
  # Thresholds (%10 minimumm training presence)
  mtp10 <- sdm_threshold.10(final_mod, occ_cal[, c('longitude', 'latitude')], 'p10', binary = F)
  bin_10 <- final_mod
  bin_10[bin_10 < minValue(mtp10)] <- NA
  bin_10[bin_10 >= minValue(mtp10)] <- 4
  
  # Mosaic
  mod_mosaic <- mosaic(bg_binary, bin_5, bin_10, fun = sum)
  mod_mosaic[mod_mosaic < 1] <- NA
  
  Col <- c('#00204DFF', '#7C7B78FF', '#FFEA46FF')
  df <- data.frame(coordinates(mod_mosaic), as.data.frame(mod_mosaic))
  
  # Plot
  ggplot() +
    geom_tile(data = df, aes(x = x, y = y, fill = as.factor(layer))) + 
    scale_fill_manual(values = Col, na.value = 'white') + 
    geom_polygon(data = M.buff, aes(x = long, y = lat, group = group), size = 0.5, color = 'black', fill = NA) +
    geom_polygon(data = df_rdp, aes(x = x, y = y), color = 'black', fill = 'white', size = 0.5) +
    geom_polygon(data = df_lagoa, aes(x = x, y = y), color = 'black', fill = 'white',size = 0.1) +
    geom_polygon(data = coast, aes(x = long, y = lat, group = group), color = 'black', fill = 'grey80', size = 0.15) +
    geom_rect(data = df_mar, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'white') +
    geom_rect(data = df_tierra, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'grey80') +
    scale_y_continuous(name = NULL, breaks = c(-55, -45, -35, -25), labels = c('55º', '45º', '35º', '25º')) + 
    scale_x_continuous(name = NULL, breaks = c(-70, -60, -50, -40), labels = c('70º', '60º', '50º', '40º')) + 
    # ggsn::scalebar(x.min = min(df$x), x.max = max(df$x), y.min = min(df$y), y.max = max(df$y), transform = T,
    #                dist = 200, st.size = 2, height = 0.013, model = 'WGS84', dist_unit = 'km',
    #                border.size = 0.25, anchor = c(x = -50.5, y = -57.35)) +
    coord_equal(xlim = c(-70.2, -36.3), ylim = c(-59.05, -20.9), expand = 0) + 
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid = element_blank(), legend.position = 'none',
          panel.border = element_rect(colour = 'black', fill = NA, size = 0.5)) 
  ggsave(paste('Figure S3_', season[i], '.tiff', sep = ''), dpi = 900, width = 6.2, height = 6.6, units = 'cm')
}  


#------------------------------------------- Figure S4 ------------------------------------------

# Seasonal continuos habitat suitability

final_mod00 <- raster('2. Analyses/annual/Final_Model_Stats/Statistics_E/calibration_median.tif') # model

for(i in 1:4){
  
  final_mod0 <- raster(paste('2. Analyses/', season[i], '/Final_Model_Stats/Statistics_E/calibration_median.tif', sep = '')) # model
  final_mod0 <- crop(final_mod0, extent(final_mod0) + c(0.5, 0, 0, 0))
  final_mod <- extend(final_mod0, extent(final_mod00) + 1.75)
  occ_cal <- read.csv(paste('2. Analyses/', season[i], '/G_galeus_joint.csv', sep = '')) # calibration occurrences
  occ_ind <- read.csv(paste('2. Analyses/', season[i], '/G_galeus_indep.csv', sep = '')) # independent occurrences
  
  # Calibration areas
  occ.tot <- SpatialPoints(cbind(occ_cal$longitude, occ_cal$latitude), proj4string = crs)
  M.buff <- buffer(occ.tot, width = 1000000)
  M.buff <- crop(M.buff, extent(final_mod0))
  
  df <- data.frame(coordinates(final_mod), as.data.frame(final_mod))
  
  # Plot
  ggplot() +
    geom_tile(data = df, aes(x = x, y = y, fill = calibration_median)) + 
    scale_fill_viridis_c(option = 'E', na.value = 'white') + 
    geom_polygon(data = M.buff, aes(x = long, y = lat, group = group), size = 0.5, color = 'black', fill = NA) +
    geom_polygon(data = df_rdp, aes(x = x, y = y), color = 'black', fill = 'white', size = 0.5) +
    geom_polygon(data = df_lagoa, aes(x = x, y = y), color = 'black', fill = 'white',size = 0.1) +
    geom_polygon(data = coast, aes(x = long, y = lat, group = group), color = 'black', fill = 'grey80', size = 0.15) +
    geom_rect(data = df_mar, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'white') +
    geom_rect(data = df_tierra, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'grey80') +
    geom_point(data = occ_cal, aes(x = longitude, y = latitude), size  =  0.05, color = 'grey5', fill = 'grey5') +
    geom_point(data = occ_ind, aes(x = longitude, y = latitude), size  =  0.05, color = 'white', fill = 'white') +
    scale_y_continuous(name = NULL, breaks = c(-55, -45, -35, -25), labels = c('55º', '45º', '35º', '25º')) + 
    scale_x_continuous(name = NULL, breaks = c(-70, -60, -50, -40), labels = c('70º', '60º', '50º', '40º')) +
    ggsn::scalebar(x.min = min(df$x), x.max = max(df$x), y.min = min(df$y), y.max = max(df$y), transform = T,
                   dist = 200, st.size = 2, height = 0.013, model = 'WGS84', dist_unit = 'km',
                   border.size = 0.25, anchor = c(x = -50.5, y = -57.35)) +
    coord_equal(xlim = c(-70.2, -36.3), ylim = c(-59.05, -20.9), expand = 0) + 
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid = element_blank(), legend.position = 'none',
          panel.border = element_rect(colour = 'black', fill = NA, size = 0.5))
  ggsave(paste('Figure S4_', season[i], '.tiff', sep = ''), dpi = 900, width = 6.2, height = 6.6, units = 'cm')
}  


#------------------------------------------- Figure S5 ------------------------------------------

# Response plots for the seasonal models

# Get feature classes and regularization multipliers from results
fc <- c('lqp', 'lqp', 'lqp', 'lqp')
rm <- c(0.4, 0.1, 0.1, 0.1)

p_df <- data.frame()
a_df <- data.frame()
plot_data_df <- data.frame()

# Rio de la Plata and Lagoa dos Patos spatial polygons
df_rdp <- matrix(c(-60.2998, -58.4014, -55.8800, -58.1127, -60.2998, 
                   -34.2257, -32.5165, -33.9072, -36.0317, -34.2257), 5, 2)
df_rdp <- SpatialPolygons(list(Polygons(list(Polygon(df_rdp)), ID = 1)), proj4string = crs)
df_lagoa <- matrix(c(-51.3938, -50.2811, -50.8794, -51.4996, -51.9246, -52.0806, -52.8042, -51.3938,
                     -29.7969, -30.2501, -31.1866, -31.7214, -31.9466, -32.1433, -32.2448, -29.7969), 8, 2)
df_lagoa <- SpatialPolygons(list(Polygons(list(Polygon(df_lagoa)), ID = 1)), proj4string = crs)

for(k in 1:4){

  # Read predictors
  env <- stack(paste(season[k], '.tif', sep = ''))
  var_names <- c('Distance_to_coast', 'Depth', 'Surface_temperature', 'SST_fronts', 'Turbidity', 
                 'Primary_productivity')
  names(env) <- var_names
  
  # Subset predictors relevant to the final model
  var_names <- c('Depth', 'Surface_temperature', 'SST_fronts', 'Turbidity')
  env <- env[[var_names]]
  env[['Depth']] <- env[['Depth']] * -1
  
  # occurrences
  dat <- read.csv(paste('2. Analyses/', season[k], '/G_galeus_joint.csv', sep = ''))
  dat <- dat[, c('longitude', 'latitude')]
  
  for(j in var_names){
    
    # prepare data
    Env <- env[[j]]
    
    # calibration areas
    occ.tot <- SpatialPoints(dat, proj4string = crs) 
    occ.buff <- buffer(occ.tot, width = 1000000) 
    env.M <- crop(Env, extent(occ.buff))
    env.M <- mask(env.M, occ.buff) 
    env.M <- mask(env.M, df_rdp, inverse = T)
    env.M <- mask(env.M, df_lagoa, inverse = T)
    env.M <- stack(env.M)
    
    # background points
    set.seed(111)
    notna <- which(complete.cases(values(env.M)))
    samp <- sample(notna, 10000, replace = F)
    samplocs <- as.data.frame(xyFromCell(env.M, samp))
    
    # SWD object
    data <- prepareSWD(species = Species, p = dat, a = samplocs, env = env.M)
    
    # run maxent replicates
    folds = randomFolds(data, k = 10, only_presence = T, seed = 111)
    default_model <- train(method = 'Maxent', data = data, fc = fc[k], reg = rm[k], iter = 1000, folds = folds)
    
    # presences and absences with variable data
    p <- .get_presence(default_model@data)
    p$var <- j
    names(p)[names(p) == j] <- 'values'
    p$Season <- season[k]
    a <- .get_absence(default_model@data)
    a$var <- j
    names(a)[names(a) == j] <- 'values'
    a$Season <- season[k]
    
    pred <- as.data.frame(matrix(data = NA, nrow = dim(data@data)[1], ncol = 10))
    for(i in 1:10){
      pred[, i] <- predict(default_model@models[[i]], data = data@data, type = 'logistic')
    }
    
    # plot data
    plot_data <- as.data.frame(matrix(data = NA, nrow = dim(data@data)[1], ncol = 4))
    names(plot_data) <- c('mean', 'sd', 'max', 'min')
    plot_data$mean <- rowMeans(pred)
    plot_data$sd <- apply(pred, 1, sd)
    plot_data$max <- plot_data$mean + plot_data$sd
    plot_data$min <- plot_data$mean - plot_data$sd
    plot_data$var <- j
    plot_data$values <- data@data[, j]
    plot_data$Season <- season[k]

    # data frames
    p_df <- rbind(p_df, p)
    a_df <- rbind(a_df, a)
    plot_data_df <- rbind(plot_data_df, plot_data)
    
  }
}

# plot preparation 
order_i_want <- c('Depth', 'Surface_temperature', 'SST_fronts', 'Turbidity')

plot_data_df <- plot_data_df[!(plot_data_df$var == 'Depth' & plot_data_df$values > 1500), ]
plot_data_df <- plot_data_df[!(plot_data_df$var == 'Surface_temperature' & plot_data_df$values > 30), ]
plot_data_df <- plot_data_df[!(plot_data_df$var == 'Surface_temperature' & plot_data_df$values < 5), ]
plot_data_df <- plot_data_df[!(plot_data_df$var == 'SST_fronts' & plot_data_df$values > 1.2), ]
plot_data_df <- plot_data_df[!(plot_data_df$var == 'Turbidity' & plot_data_df$values > 1.2), ]
plot_data_df <- transform(plot_data_df, var = factor(var, levels = order_i_want), Season = factor(Season, levels = season))

a_df <- a_df[!(a_df$var == 'Depth' & a_df$values > 1500), ]
a_df <- a_df[!(a_df$var == 'Surface_temperature' & a_df$values > 30), ]
a_df <- a_df[!(a_df$var == 'Surface_temperature' & a_df$values < 5), ]
a_df <- a_df[!(a_df$var == 'SST_fronts' & a_df$values > 1.2), ]
a_df <- a_df[!(a_df$var == 'Turbidity' & a_df$values > 1.2), ]
a_df <- transform(a_df, var = factor(var, levels = order_i_want), Season = factor(Season, levels = season))

p_df <- p_df[!(p_df$var == 'Depth' & p_df$values > 1500), ]
p_df <- p_df[!(p_df$var == 'Surface_temperature' & p_df$values > 30), ]
p_df <- p_df[!(p_df$var == 'Surface_temperature' & p_df$values < 5), ]
p_df <- p_df[!(p_df$var == 'SST_fronts' & p_df$values > 1.2), ]
p_df <- p_df[!(p_df$var == 'Turbidity' & p_df$values > 1.2), ]
p_df <- transform(p_df, var = factor(var, levels = order_i_want), Season = factor(Season, levels = season))

# labeller
var_names <- as_labeller(c(Depth = 'Bathymetry~(m)', Turbidity = 'Coefficient~Kd490~(m^-1)', Surface_temperature = 'Surface~Temperature~(ºC)', 
                           SST_fronts = 'Surface~Temperature~fronts~(ºC)', `summer` = 'Summer', `autumn` = 'Autumn', `winter` = 'Winter', 
                           `spring` = 'Spring'), default = label_parsed)
#plot_data_df=plot_data_df[!(plot_data_df$var=='sstf' & plot_data_df$bimonth%in%c('Dec-Jan','Feb-Mar','Apr-May') & plot_data_df$values>0.5),]

# plot
ggplot(data = plot_data_df, aes(x = values, y = mean, ymin = min, ymax = max)) + 
  geom_line(colour = '#00204DFF') + 
  geom_ribbon(fill = '#00204DFF', alpha = 0.2) +
  geom_rug(data = p_df, inherit.aes = F, aes(values), sides = 't', color = '#FFEA46FF', size = 0.3) + 
  geom_rug(data = a_df, inherit.aes = F, aes(values), sides = 'b', color = '#7C7B78FF', size = 0.3) + 
  labs(x = NULL, y = 'Logistic output') + ylim(0, 1) +
  theme(panel.background = element_blank(), 
        panel.grid.minor = element_blank(),                                       
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.35),
        panel.grid.major = element_line(size = 0.2, colour = 'grey90'),
        strip.background = element_rect(fill = 'transparent'),
        strip.text = element_text(vjust = -0.3, size = 8),
        plot.margin = unit(c(0, 0, 0.2, 0.2), 'cm'),
        axis.title = element_text(size = 10)) +
  facet_grid(c('Season', 'var'), scales = 'free', labeller = var_names) 
ggsave('Figure S5.pdf', width = 19, height = 15, units = 'cm')


#------------------------------------------- Figure S6 ------------------------------------------

# Backgorund overlap tatistical significance test

# Predictors relevant to the niche
var_names <- c('Depth', 'Surface_temperature', 'SST_fronts', 'Turbidity')

# Get calibration areas for each season

# Summer
env <- stack('summer.tif')[[2:5]]
names(env) <- var_names
dat_summer <- read.csv('2. Analyses/summer/G_galeus_joint.csv')
occ.tot <- SpatialPoints(cbind(dat_summer$longitude, dat_summer$latitude), proj4string = crs) 
occ.buff <- buffer(occ.tot, width = 1000000) 
env.M <- crop(env, occ.buff)
env.M <- mask(env.M, occ.buff)
env.M_summer <- stack(env.M)

# Autumn
env <- stack('autumn.tif')[[2:5]]
names(env) <- var_names
dat_autumn <- read.csv('2. Analyses/autumn/G_galeus_joint.csv')
occ.tot <- SpatialPoints(cbind(dat_autumn$longitude, dat_autumn$latitude), proj4string = crs) 
occ.buff <- buffer(occ.tot, width = 1000000) 
env.M <- crop(env, occ.buff)
env.M <- mask(env.M, occ.buff)
env.M_autumn <- stack(env.M)

# Winter
env <- stack('winter.tif')[[2:5]]
names(env) <- var_names
dat_winter <- read.csv('2. Analyses/winter/G_galeus_joint.csv')
occ.tot <- SpatialPoints(cbind(dat_winter$longitude, dat_winter$latitude), proj4string = crs) 
occ.buff <- buffer(occ.tot, width = 1000000) 
env.M <- crop(env, occ.buff)
env.M <- mask(env.M, occ.buff)
env.M_winter <- stack(env.M)

# Spring
env <- stack('spring.tif')[[2:5]]
names(env) <- var_names
dat_spring <- read.csv('2. Analyses/spring/G_galeus_joint.csv')
occ.tot <- SpatialPoints(cbind(dat_spring$longitude, dat_spring$latitude), proj4string = crs) 
occ.buff <- buffer(occ.tot, width = 1000000) 
env.M <- crop(env, occ.buff)
env.M <- mask(env.M, occ.buff)
env.M_spring <- stack(env.M)

# PCA to reduce environmental dimensions
pcas_summer <- kuenm_rpca(variables = env.M_summer, var.scale = T, write.result = F, n.pcs = 3)
pcas_autumn <- kuenm_rpca(variables = env.M_autumn, var.scale = T, write.result = F, n.pcs = 3)
pcas_winter <- kuenm_rpca(variables = env.M_winter, var.scale = T, write.result = F, n.pcs = 3)
pcas_spring <- kuenm_rpca(variables = env.M_spring, var.scale = T, write.result = F, n.pcs = 3)

# Preparing overlap objects to perform analyses
niche_summer <- overlap_object(dat_summer, species = 'species', longitude = 'longitude', latitude = 'latitude', method = 'mve1',
                               level = 95, variables = pcas_summer$PCRasters_initial)
niche_autumn <- overlap_object(dat_autumn, species = 'species', longitude = 'longitude', latitude = 'latitude', method = 'mve1',
                               level = 95, variables = pcas_autumn$PCRasters_initial)
niche_winter <- overlap_object(dat_winter, species = 'species', longitude = 'longitude', latitude = 'latitude', method = 'mve1',
                               level = 95, variables = pcas_winter$PCRasters_initial)
niche_spring <- overlap_object(dat_spring, species = 'species', longitude = 'longitude', latitude = 'latitude', method = 'mve1',
                               level = 95, variables = pcas_spring$PCRasters_initial)

# Overlap
N.test1 <- ellipsoid_overlap(niche_summer, niche_autumn, overlap_type = 'back_union', significance_test = T, replicates = 1000)
N.test2 <- ellipsoid_overlap(niche_autumn, niche_winter, overlap_type = 'back_union', significance_test = T, replicates = 1000)
N.test3 <- ellipsoid_overlap(niche_winter, niche_spring, overlap_type = 'back_union', significance_test = T, replicates = 1000)
N.test4 <- ellipsoid_overlap(niche_spring, niche_summer, overlap_type = 'back_union', significance_test = T, replicates = 1000)

# Check p-values
summary(N.test1)
summary(N.test2)
summary(N.test3)
summary(N.test4)

# Plots
Grob1 <- ggplot() +
  geom_histogram(aes(N.test1@significance_results$union_random$Niche_1_vs_2$overlap), 
                 fill = '#00204DFF',bins = 60) + xlab('Overlap') + ylab('Frequency') +
  geom_vline(xintercept  =  quantile(N.test1@significance_results$union_random$Niche_1_vs_2$overlap, 0.05),
             colour = '#FFEA46FF', size = 0.75, linetype = 'dashed') +  
  geom_vline(xintercept  =  N.test1@union_overlap$overlap,
             colour = '#7C7B78FF', size = 0.75) + theme_bw()
Grob2 <- ggplot() +
  geom_histogram(aes(N.test2@significance_results$union_random$Niche_1_vs_2$overlap), 
                 fill = '#00204DFF',bins = 60) + xlab('Overlap') + ylab('Frequency') +
  geom_vline(xintercept  =  quantile(N.test2@significance_results$union_random$Niche_1_vs_2$overlap, 0.05),
             colour = '#FFEA46FF', size = 0.75, linetype = 'dashed') +  
  geom_vline(xintercept  =  N.test2@union_overlap$overlap,
             colour = '#7C7B78FF', size = 0.75) + theme_bw()
Grob3 <- ggplot() +
  geom_histogram(aes(N.test3@significance_results$union_random$Niche_1_vs_2$overlap), 
                 fill = '#00204DFF',bins = 60) + xlab('Overlap') + ylab('Frequency') +
  geom_vline(xintercept  =  quantile(N.test3@significance_results$union_random$Niche_1_vs_2$overlap, 0.05),
             colour = '#FFEA46FF', size = 0.75, linetype = 'dashed') +  
  geom_vline(xintercept  =  N.test3@union_overlap$overlap,
             colour = '#7C7B78FF', size = 0.75) + theme_bw()
Grob4 <- ggplot() +
  geom_histogram(aes(N.test4@significance_results$union_random$Niche_1_vs_2$overlap), 
                 fill = '#00204DFF',bins = 60) + xlab('Overlap') + ylab('Frequency') +
  geom_vline(xintercept  =  quantile(N.test4@significance_results$union_random$Niche_1_vs_2$overlap, 0.05),
             colour = '#FFEA46FF', size = 0.75, linetype = 'dashed') +  
  geom_vline(xintercept  =  N.test4@union_overlap$overlap,
             colour = '#7C7B78FF', size = 0.75) + theme_bw()

ggsave('Figure S6.pdf', plot = grid.arrange(Grob1, Grob2, Grob3, Grob4, ncol = 4), 
       width = 30, height = 8, units = 'cm')


#------------------------------------------- Figure S7 ------------------------------------------

# Go to script 'G_galeus_habitat_associations.R'


#------------------------------------------- Figure S8 ------------------------------------------

# Go to script 'G_galeus_habitat_associations.R'


#------------------------------------------- Appendix 1 ------------------------------------------

# Recaptures map

# Calibration areas
occ_cal <- read.csv('2. Analyses/annual/G_galeus_joint.csv') # calibration occurrences
occ.tot <- SpatialPoints(cbind(occ_cal$longitude, occ_cal$latitude), proj4string = crs)
M.buff <- buffer(occ.tot, width = 1000000)
df_mar2 <- data.frame(xmin = -80, xmax = -30, ymin = -70, ymax = -10)

recaptures <- data.frame(lon = c(-64.88, -62.11, -57.3, -52.5), 
                         lat = c(-42.72, -40.33, -38.4, -35.3))

# Plot
ggplot() +
  geom_rect(data = df_mar2, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'white') +
  geom_polygon(data = M.buff, aes(x = long, y = lat, group = group), size = 0.75, color = 'black', fill = '#00204DFF') +
  geom_polygon(data = df_rdp, aes(x = x, y = y), color = 'black', fill = 'white', size = 0.75) +
  geom_polygon(data = df_lagoa, aes(x = x, y = y), color = 'black', fill = 'white',size = 0.1) +
  geom_polygon(data = coast, aes(x = long, y = lat, group = group), color = 'black', fill = 'grey80', size = 0.15) +
  geom_rect(data = df_mar, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'white') +
  geom_rect(data = df_tierra, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'grey80') +
  geom_point(data = recaptures, aes(x = lon, y = lat), size  =  0.0001, color = '#FFEA46FF', fill = '#FFEA46FF') +
  scale_y_continuous(name = NULL, breaks = c(-45, -40, -35, -30), labels = c('45º', '40º', '35º', '30º')) + 
  scale_x_continuous(name = NULL, breaks = c(-65, -60, -55, -50), labels = c('65º', '60º', '55º', '50º')) +
  ggsn::scalebar(x.min = -66, x.max = -49, y.min = -46, y.max = -29, transform = T,
                 dist = 100, st.size = 2.5, height = 0.013, model = 'WGS84', dist_unit = 'km',
                 border.size = 0.5, anchor = c(x = -55, y = -45)) +
  coord_equal(xlim = c(-66, -49), ylim = c(-46, -29), expand = 0) + 
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(), legend.position = 'none',
        panel.border = element_rect(colour = 'black', fill = NA, size = 0.5)) 
ggsave('Appendix 1.pdf', width = 9.1, height = 10, units = 'cm')


#--------------------------------------------- END -------------------------------------





