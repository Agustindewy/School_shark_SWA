
# De Wysiecki et al. 
# "Population-scale habitat use of school sharks (Triakidae: Galeorhinus galeus) in the Southwest Atlantic: 
# insights from temporally explicit niche modelling and habitat associations"

# Predictors download and preparation

library(raster)
library(rgdal)
library(sqldf)
library(maps)
library(ggplot2)
library(grec)
library(ncdf4)

setwd('SET YOUR WORKING DIRECTORY')

# Projection
crs <- CRS('+init=epsg:4326')

# Species
Species <- 'Galeorhinus galeus'

# SWA coastline (spatial polygons) - taken from GSHHG coastline database
coast <- readOGR(dsn = 'DOWNLOAD AND READ THE RASTER', layer = 'GSHHS_f_L1_SouthAmerica')

# Season
season <- c('summer', 'autumn', 'winter', 'spring')

#-------------------------------------------- Extent ------------------------------------------------

# Data download extent based on calibration areas (lon.min, lon.max, lat.min, lat.max)
ext <- c(-71, -35, -60, -12)

# Working grid based on extreme records
bb <- matrix(c(-70, -60, -35, -12), 2, 2) # bounding box
cs <- c(0.041, 0.041) # cell size
cc <- bb[, 1] + (cs / 2) # extra space to include all points
cd <- ceiling(diff(t(bb)) / cs) # x-y cell dimension
gri <- GridTopology(cellcentre.offset = cc, cellsize = cs, cells.dim = cd) # grid
sp.gri <- SpatialGrid(gri, proj4string = crs) # spatial grid
sp.ras <- raster(sp.gri) # base raster


#-------------------------------------- Static predictors ------------------------------------------------

# Distance to coast (kilometers) - taken from MARSPEC (30 arcseconds, ~1 km)
dc <- raster('DOWNLOAD AND READ THE RASTER/Dist_30s.tif')
dc <- crop(dc, ext) # crop to broader study area
dc <- projectRaster(dc, sp.ras)
m <- matrix(1, ncol = 5, nrow = 5) # moving window
dc <- focal(dc, m, fun = mean, NAonly = T, pad = T, na.rm = T) # interpolate for missing values
dc <- mask(dc, coast, inverse = T) # tidy up with coastline

# Bathymetry (meters) - taken from MARSPEC (30 arcseconds, ~1 km)
dep <- raster('DOWNLOAD AND READ THE RASTER/Bathy_30s.tif')
dep <- crop(dep, dc) # crop to study area
dep <- projectRaster(dep, dc)
dep <- focal(dep, m, fun = mean, NAonly = T, pad = T, na.rm = T) # interpolate for missing values
dep <- mask(dep, dc) #mask to study area

static_pred <- stack(dc, dep)
#writeRaster(static_pred, filename = 'static_pred.tif') #write as tiff


#-------------------------------------- Dynamic predictors ------------------------------------------------

#------------------------------------ Sea surface temperature ----------------------------------------

# ºC, 01-2003 to 12-2020 monthly, ~4 km)

# Gotta download and read in 2 periods because of big data size (1 = 2003-2011, 2 = 2012-2020)
sst <- read.csv('DOWNLOAD AND READ THE CSV/sst_1.csv', skip = 2, header = F)
sst <- sst[which(!is.na(sst$V4)), ] # erase NA values
sst$Year <- substr(sst$V1, 1, 4)
sst$Month <- substr(sst$V1, 6, 7)
sst$Month <- as.numeric(sst$Month)
# sst$Season1 <- ''
sst$Season2 <- ''
sst <- sst[, -1]

sst$Season1[which(sst$Month < 4)] <- 'Summer'
sst$Season1[which(sst$Month > 3 & sst$Month < 7)] <- 'Autumn'
sst$Season1[which(sst$Month > 6 & sst$Month < 10)] <- 'Winter'
sst$Season1[which(sst$Month > 9)] <- 'Spring'

sst$Season2[which(sst$Month > 11 | sst$Month < 3)] <- 'summer'
sst$Season2[which(sst$Month > 2 & sst$Month < 6)] <- 'autumn'
sst$Season2[which(sst$Month > 5 & sst$Month < 9)] <- 'winter'
sst$Season2[which(sst$Month > 8 & sst$Month < 12)] <- 'spring'

# Annual
sst0 <- data.frame(lat = sst$V2, lon = sst$V3, sst = sst$V4)
sst.ras <- SpatialPoints(cbind(sst0$lon, sst0$lat), proj4string = crs) # transform to spatial points
sst.ras <- rasterize(sst.ras, sp.ras, sst0$sst, fun = mean) # rasterize
sst.ras <- projectRaster(sst.ras, sp.ras)
m <- matrix(1, ncol = 5, nrow = 5)
sst.ras <- focal(sst.ras, m, fun = mean, NAonly = T, pad = T, na.rm = T) # interpolate for missing values
sst.ras <- mask(sst.ras, static_pred[[1]]) # mask to study area
writeRaster(sst.ras, filename = 'sst_annual_1.tif') # write as tif

# Season 1
sst0 <- data.frame(Season = sst$Season1, lat = sst$V2, lon = sst$V3, sst = sst$V4)
sst0 <- aggregate(sst ~ Season + lon + lat, data = sst0, mean)
sst.list <- vector('list', 4)
for(i in 1:4){
  sub <- subset(sst0, Season == season[i])
  sst.ras <- SpatialPoints(cbind(sub$lon, sub$lat), proj4string = crs)
  sst.ras <- rasterize(sst.ras, sp.ras, sub$sst, fun = mean)
  sst.ras <- projectRaster(sst.ras, sp.ras)
  sst.ras <- focal(sst.ras, m, fun = mean, NAonly = T, pad = T, na.rm = T)
  sst.ras <- mask(sst.ras, static_pred[[1]])
  sst.list[[i]] <- sst.ras
}
sst.season <- stack(sst.list)
writeRaster(sst.season, filename = 'sst_season1_1.tif') # write as tiff

# Season 2
sst0 <- data.frame(Season = sst$Season2, lat = sst$V2, lon = sst$V3, sst = sst$V4)
sst0 <- aggregate(sst ~ Season + lon + lat, data = sst0, mean)
sst.list <- vector('list', 4)
for(i in 1:4){
  sub <- subset(sst0, Season == season[i])
  sst.ras <- SpatialPoints(cbind(sub$lon, sub$lat), proj4string = crs)
  sst.ras <- rasterize(sst.ras, sp.ras, sub$sst, fun = mean) 
  sst.ras <- projectRaster(sst.ras, sp.ras)
  sst.ras <- focal(sst.ras, m, fun = mean, NAonly = T, pad = T, na.rm = T)
  sst.ras <- mask(sst.ras, static_pred[[1]]) 
  sst.list[[i]] <- sst.ras
}
sst.season <- stack(sst.list)
writeRaster(sst.season, filename = 'sst_season2_1.tif') # write as tiff


#------------------------------------- SST temporal merging ----------------------------------------

# Annual
sst_1 <- raster('sst_annual_1.tif')
sst_2 <- raster('sst_annual_2.tif')
sst <- mean(sst_1, sst_2)
writeRaster(sst, filename = 'sst.tif')

# Season 1
sst_1 <- stack('sst_season1_1.tif')
sst_2 <- stack('sst_season1_2.tif')
sst_season <- stack(mean(sst_1[[1]], sst_2[[1]]), mean(sst_1[[2]], sst_2[[2]]),
                    mean(sst_1[[3]], sst_2[[3]]), mean(sst_1[[4]], sst_2[[4]]))
writeRaster(sst_season, filename = 'sst_season1.tif')

# Season 2
sst_1 <- stack('sst_season2_1.tif')
sst_2 <- stack('sst_season2_2.tif')
sst_season <- stack(mean(sst_1[[1]], sst_2[[1]]), mean(sst_1[[2]], sst_2[[2]]), 
                    mean(sst_1[[3]], sst_2[[3]]), mean(sst_1[[4]], sst_2[[4]]))
writeRaster(sst_season, filename = 'sst_season2.tif') 


#------------------------------------- Season 1 or Season 2 ? ------------------------------------

# Check which seasonal aggrupation reflects greater differences between seasonal occurrences
# Idea is to show the greatest possible niche differenciation between seasons

# Read sea surface temperature and distance to coast
sst_1 <- stack('sst_season1.tif')
sst_2 <- stack('sst_season2.tif')
dc <- stack('static_pred.tif')[[1]]
#env <- stack(sst_1, sst_2, dc)

# Read occurrences
occs <- read.csv('Presences.csv')
dups <- duplicated(occs[c('Longitude', 'Latitude')]) # Get rid of duplicates
occs <- occs[!dups, ]

# Processing occurrences
# Save any points to nearest pixel with data that may be falling outside a layer of the M space raster 
library(rSDM) # package in GitHub (Pakillo/rSDM), install and call 'rSDM'
spp <- SpatialPoints(occs[, c('Longitude', 'Latitude')], crs)
spp_corrected <- points2nearestcell(locs = spp, ras = dc, layer = 1, move = T) # plots will appear if any correction applies
occs$Longitude <- round(spp_corrected@coords[, 1], 2) # replace coordinates including those new if any
occs$Latitude <- round(spp_corrected@coords[, 2], 2) 

# Environmental filtering and extraccion
source('envSample.R')
sst1.df <- data.frame()
sst2.df <- data.frame()
env.M1.df <- data.frame()
env.M2.df <- data.frame()
for(i in 1:4){
  # Occurrences
  sub1 <- subset(occs, Season1 == season[i])
  sub2 <- subset(occs, Season2 == season[i])
  coords1 <- sub1[, c('Longitude', 'Latitude')]
  coords2 <- sub2[, c('Longitude', 'Latitude')]
  sst1.data <- data.frame(sst = extract(sst_1[[i]], coords1))
  sst2.data <- data.frame(sst = extract(sst_2[[i]], coords2))
  dc.data1 <- data.frame(dist = extract(dc, coords1))
  dc.data2 <- data.frame(dist = extract(dc, coords2))
  coords1 <- envSample(coords1, filters = list(sst1.data$sst, dc.data1$dist),
                       res = list(0.5, 5), do.plot = T)
  coords2 <- envSample(coords2, filters = list(sst2.data$sst, dc.data2$dist),
                       res = list(0.5, 5), do.plot = T)
  sst1.extract <- data.frame(sst = extract(sst_1[[i]], coords1), Longitude = coords1$lon, Latitude = coords1$lat)
  sst2.extract <- data.frame(sst = extract(sst_2[[i]], coords2), Longitude = coords2$lon, Latitude = coords2$lat)
  sst1.extract$season <- season[[i]]
  sst2.extract$season <- season[[i]]
  sst1.df <- rbind(sst1.df, sst1.extract)
  sst2.df <- rbind(sst2.df, sst2.extract)
  
  # M space 1 
  occ.tot <- SpatialPoints(cbind(sst1.df$Longitude, sst1.df$Latitude), proj4string = crs)
  occ.buff <- buffer(occ.tot, width = 1000000)
  env.M1 <- crop(sst_1[[i]], occ.buff)
  env.M1 <- mask(env.M1, occ.buff)
  variables_values <- data.frame(na.omit(values(env.M1)))
  variables_values <- data.frame(sst = variables_values[sample(1:dim(variables_values)[1], 10000), ])
  variables_values$season <- season[[i]]
  env.M1.df <- rbind(env.M1.df, variables_values)
  
  # M space 2
  occ.tot <- SpatialPoints(cbind(sst2.df$Longitude, sst2.df$Latitude), proj4string = crs)
  occ.buff <- buffer(occ.tot, width = 1000000)
  env.M2 <- crop(sst_2[[i]], occ.buff)
  env.M2 <- mask(env.M2, occ.buff)
  variables_values <- data.frame(na.omit(values(env.M2)))
  variables_values <- data.frame(sst = variables_values[sample(1:dim(variables_values)[1], 10000), ])
  variables_values$season <- season[[i]]
  env.M2.df <- rbind(env.M2.df, variables_values)
}

sst1.df$logical <- 'jan-feb-mar'
sst2.df$logical <- 'dec-jan-feb'
sst.df <- rbind(sst1.df, sst2.df)
sst.df <- transform(sst.df, season = factor(season, levels = c('summer', 'autumn', 'winter', 'spring')))

env.M1.df$logical <- 'jan-feb-mar'
env.M2.df$logical <- 'dec-jan-feb'
env.M.df <- rbind(env.M1.df, env.M2.df)
env.M.df <- transform(env.M.df, season = factor(season, levels = c('summer', 'autumn', 'winter', 'spring')))

p <- ggplot(sst.df, aes(factor(season), sst, fill = factor(logical)))
p + geom_boxplot()

p <- ggplot(env.M.df, aes(factor(season), sst, fill = factor(logical)))
p + geom_boxplot()

# Use season 2 (dec-jan-feb) because is the one with the most inter-seasonal variability


#--------------------------------------- Coefficient Kd490 ---------------------------------------

# 1/m, 01-2003 to 12-2020 monthly, ~4 km)

# Gotta download and read in 2 periods because of big data size (1 = 2003-2011, 2 = 2012-2020)
tur <- read.csv('DOWNLOAD AND READ THE CSV/tur_1.csv', skip = 2, header = F)
tur <- tur[which(!is.na(tur$V4)), ] # erase NA values
tur$Year <- substr(tur$V1, 1, 4)
tur$Month <- substr(tur$V1, 6, 7)
tur$Month <- as.numeric(tur$Month)
tur$Season2 <- ''
tur <- tur[, -1]

tur$Season2[which(tur$Month > 11 | tur$Month < 3)] <- 'summer'
tur$Season2[which(tur$Month > 2 & tur$Month < 6)] <- 'autumn'
tur$Season2[which(tur$Month > 5 & tur$Month < 9)] <- 'winter'
tur$Season2[which(tur$Month > 8 & tur$Month < 12)] <- 'spring'

# Annual
tur0 <- data.frame(lat = tur$V2, lon = tur$V3, tur = tur$V4)
tur.ras <- SpatialPoints(cbind(tur0$lon, tur0$lat), proj4string = crs) # transform to spatial points
tur.ras <- rasterize(tur.ras, sp.ras, tur0$tur, fun = mean) # rasterize
tur.ras <- projectRaster(tur.ras, sp.ras)
m <- matrix(1, ncol = 5, nrow = 5)
tur.ras <- focal(tur.ras, m, fun = mean, NAonly = T, pad = T, na.rm = T) # interpolate for missing values
tur.ras <- mask(tur.ras, static_pred[[1]]) # mask to study area
writeRaster(tur.ras, filename = 'tur_annual_1.tif') # write as tif

# Season 2
tur0 <- data.frame(Season = tur$Season2, lat = tur$V2, lon = tur$V3, tur = tur$V4)
tur0 <- aggregate(tur ~ Season + lon + lat, data = tur0, mean)
tur.list <- vector('list', 4)
for(i in 1:4){
  sub <- subset(tur0, Season == season[i])
  tur.ras <- SpatialPoints(cbind(sub$lon, sub$lat), proj4string = crs)
  tur.ras <- rasterize(tur.ras, sp.ras, sub$tur, fun = mean) 
  tur.ras <- projectRaster(tur.ras, sp.ras)
  tur.ras <- focal(tur.ras, m, fun = mean, NAonly = T, pad = T, na.rm = T)
  tur.ras <- mask(tur.ras, static_pred[[1]]) 
  tur.list[[i]] <- tur.ras
}
tur.season <- stack(tur.list)
writeRaster(tur.season, filename = 'tur_season2_1.tif') # write as tiff


#------------------------------------- Primary productivity ------------------------------------

# mgC/m2/day, 01-2003 to 12-2020 monthly, ~4 km)

# Gotta download and read in 2 periods because of big data size (1 = 2003-2011, 2 = 2012-2020)
pro <- read.csv('DOWNLOAD AND READ THE CSV/pp_1.csv', skip = 2, header = F)
pro <- pro[which(!is.na(pro$V4)), ] # erase NA values
pro$Year <- substr(pro$V1, 1, 4)
pro$Month <- substr(pro$V1, 6, 7)
pro$Month <- as.numeric(pro$Month)
pro$Season2 <- ''
pro <- pro[, -c(1:2)]

pro$Season2[which(pro$Month > 11 | pro$Month < 3)] <- 'summer'
pro$Season2[which(pro$Month > 2 & pro$Month < 6)] <- 'autumn'
pro$Season2[which(pro$Month > 5 & pro$Month < 9)] <- 'winter'
pro$Season2[which(pro$Month > 8 & pro$Month < 12)] <- 'spring'

# Annual
pro0 <- data.frame(lat = pro$V3, lon = pro$V4, pro = pro$V5)
pro.ras <- SpatialPoints(cbind(pro0$lon, pro0$lat), proj4string = crs) # transform to spatial points
pro.ras <- rasterize(pro.ras, sp.ras, pro0$pro, fun = mean) # rasterize
pro.ras <- projectRaster(pro.ras, sp.ras)
m <- matrix(1, ncol = 5, nrow = 5)
pro.ras <- focal(pro.ras, m, fun = mean, NAonly = T, pad = T, na.rm = T) # interpolate for missing values
pro.ras <- mask(pro.ras, static_pred[[1]]) # mask to study area
writeRaster(pro.ras/1000, filename = 'pro_annual_1.tif') # write as tif

# Season 2
pro0 <- data.frame(Season = pro$Season2, lat = pro$V3, lon = pro$V4, pro = pro$V5)
pro0 <- aggregate(pro ~ Season + lon + lat, data = pro0, mean)
pro.list <- vector('list', 4)
for(i in 1:4){
  sub <- subset(pro0, Season == season[i])
  pro.ras <- SpatialPoints(cbind(sub$lon, sub$lat), proj4string = crs)
  pro.ras <- rasterize(pro.ras, sp.ras, sub$pro, fun = mean) 
  pro.ras <- projectRaster(pro.ras, sp.ras)
  pro.ras <- focal(pro.ras, m, fun = mean, NAonly = T, pad = T, na.rm = T)
  pro.ras <- mask(pro.ras, static_pred[[1]]) 
  pro.list[[i]] <- pro.ras
}
pro.season <- stack(pro.list)
writeRaster(pro.season/1000, filename = 'pro_season2_1.tif') # write as tiff


#----------------------------------- TUR & PRO Temporal merging -------------------------------------

# Annual
tur_1 <- raster('tur_annual_1.tif')
tur_2 <- raster('tur_annual_2.tif')
tur <- mean(tur_1, tur_2)
writeRaster(tur, filename = 'tur.tif')
pro_1 <- raster('pro_annual_1.tif')
pro_2 <- raster('pro_annual_2.tif')
pro <- mean(pro_1, pro_2)
writeRaster(pro, filename = 'pro.tif')

# Season 2
tur_1 <- stack('tur_season2_1.tif')
tur_2 <- stack('tur_season2_2.tif')
tur_season <- stack(mean(tur_1[[1]], tur_2[[1]]), mean(tur_1[[2]], tur_2[[2]]), 
                    mean(tur_1[[3]], tur_2[[3]]), mean(tur_1[[4]], tur_2[[4]]))
writeRaster(tur_season, filename = 'tur_season2.tif') 
pro_1 <- stack('pro_season2_1.tif')
pro_2 <- stack('pro_season2_2.tif')
pro_season <- stack(mean(pro_1[[1]], pro_2[[1]]), mean(pro_1[[2]], pro_2[[2]]), 
                    mean(pro_1[[3]], pro_2[[3]]), mean(pro_1[[4]], pro_2[[4]]))
writeRaster(pro_season, filename = 'pro_season2.tif') 


#------------------------------------------ SST fronts ------------------------------------------------

# Derived from sea surface temperature, in ºC

# Annual
sst <- raster('sst.tif') # read annual sst raster
Front <- detectFronts(sst, method = 'median_filter')
Front <- focal(Front, m, fun = mean, NAonly = T, pad = T, na.rm = T) # interpolate for missing values
Front <- mask(Front, sst) # mask to study area
writeRaster(Front, filename = 'sstf.tif') # write as tiff

# Season 2
sst_season <- stack('sst_season2.tif') # read seasonal sst rasters
Front0 <- stack()
for(i in 1:4){
  sst <- sst_season[[i]]
  Front <- detectFronts(sst)
  Front <- focal(Front, m, fun = mean, NAonly = T, pad = T, na.rm = T) # interpolate for missing values
  Front <- mask(Front, sst) # mask to study area
  Front0 <- stack(Front0, Front)
}
writeRaster(Front0, filename = 'sstf_season2.tif') # write as tiff


#----------------------------------------- Final stacks ------------------------------------------

Names <- c('Distance_to_coast', 'Depth', 'Surface_temperature', 'SST_fronts', 'Turbidity', 'Primary_productivity')

# Annual
static_pred <- stack('static_pred.tif')
sst <- raster('sst.tif')
sstf <- raster('sstf.tif')
tur <- raster('tur.tif')
pro <- raster('pro.tif')
annual.stack <- stack(static_pred, sst, sstf, tur, pro)
names(annual.stack) <- Names
writeRaster(annual.stack, filename = 'env.tif') # write as tiff

# Season 2
sst <- stack('sst_season2.tif')
sstf <- stack('sstf_season2.tif')
tur <- stack('tur_season2.tif')
pro <- stack('pro_season2.tif')
for(i in 1:4){
  sst0 <- sst[[i]]
  sstf0 <- sstf[[i]]
  tur0 <- tur[[i]]
  pro0 <- pro[[i]]
  seasonal.stack <- stack(static_pred, sst0, sstf0, tur0, pro0)
  names(seasonal.stack) <- Names
  writeRaster(seasonal.stack, filename = paste(season[i], '.tif', sep = '')) # write as tiff
}

#---------------------------------------------- END ------------------------------------------
