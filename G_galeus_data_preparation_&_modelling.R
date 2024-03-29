
# De Wysiecki et al. 
# "Population-scale habitat use of school sharks (Triakidae: Galeorhinus galeus) in the Southwest Atlantic: 
# insights from temporally explicit niche modelling and habitat associations"

# Data preparation and modelling

library(raster)
library(rgdal)
library(sqldf)
library(maps)
library(ggplot2)
library(rgeos)
library(kuenm)

setwd('SET YOUR WORKING DIRECTORY')

# Projection
crs <- CRS('+init=epsg:4326')

# Species
Species <- 'Galeorhinus galeus'

# Season
season <- c('summer', 'autumn', 'winter', 'spring')

# Seed
set.seed(111)

# Read occurrence data
dat <- read.csv('Presences.csv')

# Rio de la Plata and Lagoa dos Patos spatial polygons
df_rdp <- matrix(c(-60.2998, -58.4014, -55.8800, -58.1127, -60.2998, 
                   -34.2257, -32.5165, -33.9072, -36.0317, -34.2257), 5, 2)
df_rdp <- SpatialPolygons(list(Polygons(list(Polygon(df_rdp)), ID = 1)), proj4string = crs)
df_lagoa <- matrix(c(-51.3938, -50.2811, -50.8794, -51.4996, -51.9246, -52.0806, -52.8042, -51.3938,
                     -29.7969, -30.2501, -31.1866, -31.7214, -31.9466, -32.1433, -32.2448, -29.7969), 8, 2)
df_lagoa <- SpatialPolygons(list(Polygons(list(Polygon(df_lagoa)), ID = 1)), proj4string = crs)

#-------------------------------------- Annual analysis --------------------------------------------

temporal_block <- 'annual'
dir.create(temporal_block)

#------------------------------------------ G space --------------------------------------------

# Read predictors
env <- stack('env.tif')
var_names <- c('Distance_to_coast', 'Depth', 'Surface_temperature', 'SST_fronts', 'Turbidity', 
               'Primary_productivity')
names(env) <- var_names

# Subset predictors relevant to G. galeus
var_names <- c('Distance_to_coast', 'Depth', 'Surface_temperature', 'SST_fronts', 'Turbidity', 
               'Primary_productivity')

# Rio de la Plata and Lagoa dos Patos not suitable
env <- mask(env, df_rdp, inverse = T)
env <- mask(env, df_lagoa, inverse = T)

#-------------------------------------- Data preparation ----------------------------------

# Initial M space for processing occurrences
coord <- cbind(dat$Longitude, dat$Latitude)
occ.tot <- SpatialPoints(coord, proj4string = crs) # transform to spatial points
occ.buff <- buffer(occ.tot, width = 1000000) # calibration areas: 1000 km buffer around each point
env.M <- crop(env[[var_names]], occ.buff) # crop raster stack to buffer polygon
env.M <- mask(env.M, occ.buff) # mask raster stack to buffer polygon
env.M <- stack(env.M)

# Check correlation in M space
layerStats(env.M, 'pearson', na.rm = T) 
# depth and dist to coast > 0.8, keep depth
# turbidity and primary productivity > 0.9, keep turbidity
var_set <- var_names[-c(1, 6)]
env.M <- env.M[[var_set]]

# Processing occurrences
# Save any points to nearest pixel with data that may be falling outside a layer of the M space raster 
library(rSDM) #package in GitHub (Pakillo/rSDM), install and call 'rSDM'
for(i in 1:length(env.M@layers)){
  spp <- SpatialPoints(dat[, c('Longitude', 'Latitude')], crs)
  spp_corrected <- points2nearestcell(locs = spp, ras = env.M, layer = i, move = T) 
  dat$Longitude <- round(spp_corrected@coords[, 1], 2)  # replace coordinates including those new if any
  dat$Latitude <- round(spp_corrected@coords[, 2], 2)
} # plots will appear if any correction applies

# Independent data
occ_ind <- subset(dat, Type == 'independent')
dups <- duplicated(occ_ind[c('Longitude', 'Latitude')])# Get rid of duplicates
occ_ind <- occ_ind[!dups, ]
occ_ind0 <- data.frame(species = Species, longitude = occ_ind$Longitude, latitude = occ_ind$Latitude, Type = 'independent')
write.csv(occ_ind0[, 1:3], paste(temporal_block, '/G_galeus_indep.csv', sep = ''), row.names = F)

# Calibration records
occ_cal <- subset(dat, Type == 'calibration')
dups <- duplicated(occ_cal[c('Longitude', 'Latitude')])# Get rid of duplicates
occ_cal <- occ_cal[!dups, ]

aaa <- subset(occ_cal, Latitude > -31)
aaa <- aaa[sample(nrow(aaa), 10), ]
occ_cal <- subset(occ_cal, Latitude < -31)
occ_cal <- rbind(occ_cal, aaa)

# Environmental filtering
source('envSample.R')
coords <- occ_cal[, c('Longitude', 'Latitude')]
env.data <- extract(env.M, coords)
env.data <- as.data.frame(env.data)
coords <- envSample(coords, filters = list(env.data$Surface_temperature, env.data$Depth),
                    res = list(0.5, 10), do.plot = T)
occ_cal <- merge(occ_cal, coords, by.x = c('Longitude', 'Latitude'), by.y = c('lon', 'lat'))

# Detecting environmental outliers in species occurrences (Cobos et al., 2018)
# Getting data from the variables
variables_values <- na.omit(values(env.M)) # for the region of interest
occ_variables <- na.omit(cbind(as.numeric(row.names(occ_cal)), extract(env.M, occ_cal[, c('Longitude', 'Latitude')])))
occ_variables <- as.data.frame(occ_variables)
if(dim(variables_values)[1] > 10000) { # sample of 10000 values if more pixels exist
  variables_values <- variables_values[sample(1:dim(variables_values)[1], 10000), ]
}

# Plot for searching for potential environmental outliers
i <- 2
j <- 4
par(mfcol = c(1, 2), mar = c(4.5, 4, 0.5, 0.5), cex = 0.9)
plot(variables_values[, i], variables_values[, j], col = 'grey65', pch = 1,
     xlab = colnames(variables_values)[i], ylab = colnames(variables_values)[j])
points(occ_variables[, i + 1], occ_variables[, j + 1], col = 'blue', pch = 19, cex = 1.5)
legend('bottomright', legend = c('Region of interest', 'Occurrences'),
       pch = c(1, 19), col = c('gray65', 'black'), bty = 'n')
plot(variables_values[, i], variables_values[, j], col = 'gray65',
     pch = 1, xlab = colnames(variables_values)[i], ylab = colnames(variables_values)[j])
text(occ_variables[, i + 1], occ_variables[, j + 1], occ_variables[, 1], cex = 1, col = 'blue')
legend('bottomright', legend = 'Occurrence ID', bty = 'n')

# Remove outliers from calibration records
occ_variables <- subset(occ_variables, !V1 %in% c(2, 3, 1))
#c(138, 136, 233, 173, 305)
occ_cal <- occ_cal[which(row.names(occ_cal) %in% as.character(occ_variables$V1)), ]

# Final corrected points (calibration + independent)
occ_corrected <- rbind(occ_cal, occ_ind)
write.csv(occ_corrected, paste(temporal_block, '/G_galeus_corrected.csv', sep = ''), row.names = F)

# Occurrence subsets for traning, testing and independent evaluation
# Create csv files for training and testing (75% and 25%)
occ_cal <- data.frame(species = occ_cal$Species, longitude = occ_cal$Longitude, latitude = occ_cal$Latitude)
occ_cal$check <- paste(occ_cal[, 'longitude'], occ_cal[, 'latitude'], sep = '_')
train <- occ_cal[sample(nrow(occ_cal), round((length(occ_cal[, 1]) / 10 * 7.5))), ]
test <- occ_cal[!occ_cal[, 4] %in% train[, 4], ]
occ_cal$check <- NULL; train$check <- NULL; test$check <- NULL
write.csv(occ_cal[, 1:3], paste(temporal_block, '/G_galeus_joint.csv', sep = ''), row.names = F)
write.csv(train[, 1:3], paste(temporal_block, '/G_galeus_train.csv', sep = ''), row.names = F)
write.csv(test[, 1:3], paste(temporal_block, '/G_galeus_test.csv', sep = ''), row.names = F)

# Final M space (calibration areas)
coord <- cbind(occ_cal$longitude, occ_cal$latitude)
occ.tot <- SpatialPoints(coord, proj4string = crs) 
occ.buff <- buffer(occ.tot, width = 1000000) 
env.M <- crop(env[[var_set]], extent(occ.buff) + 1)
env.M <- mask(env.M, occ.buff) 
env.M <- mask(env.M, df_rdp, inverse = T)
env.M <- mask(env.M, df_lagoa, inverse = T)
env.M <- stack(env.M)

# Create M space for modelling
dir.create(paste(temporal_block, '/M_variables', sep = ''))
dir.create(paste(temporal_block, '/M_variables/Set1', sep = ''))
writeRaster(env.M[[var_set]], filename = paste(temporal_block, '/M_variables/Set1/env.asc', sep = ''),
            format = 'ascii', bylayer = T, suffix = var_set)
dir.create(paste(temporal_block, '/G_variables', sep = ''))
dir.create(paste(temporal_block, '/G_variables/Set1', sep = ''))
dir.create(paste(temporal_block, '/G_variables/Set1/current', sep = ''))
writeRaster(env[[var_set]], filename = paste(temporal_block, '/G_variables/Set1/current/env.asc', sep = ''),
            format = 'ascii', bylayer = T, suffix = var_set)

#------------------------------------------ Modelling -------------------------------------------

# Seed
set.seed(111)

#--------------------------------------- Candidate models --------------------------------------------

occ_joint <- paste(temporal_block, '/G_galeus_joint.csv', sep = '')
occ_tra <- paste(temporal_block, '/G_galeus_train.csv', sep = '')
M_var_dir <- paste(temporal_block, '/M_variables', sep = '')
batch_cal <- paste(temporal_block, '/Candidate_Models', sep = '')
out_dir <- paste(temporal_block, '/Candidate_Models', sep = '')
reg_mult <- c(seq(0.1, 2, 0.1), 2.5, 3, 4, 5, 7, 10) # finer set of values
f_clas <- c('lq', 'lp', 'lqp')
args <- NULL
maxent_path <- 'MAXENT JAR DIRECTORY'
wait <- TRUE
run <- TRUE

kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal, 
          out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, args = args, maxent.path = maxent_path, 
          wait = wait, run = run)

#----------------------------- Evaluation and selection of best models ------------------------

occ_test <- paste(temporal_block, '/G_galeus_test.csv', sep = '')
out_eval <- paste(temporal_block, '/Calibration_Results', sep = '')
threshold <- 5
rand_percent <- 50
iterations <- 500
kept <- FALSE
selection <- 'OR_AICc' # OR_AICc
paral_proc <- FALSE 

kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test, 
            batch = batch_cal, out.eval = out_eval, threshold = threshold, rand.percent = rand_percent, 
            iterations = iterations, kept = kept, selection = selection, parallel.proc = paral_proc)

#-------------------------------------- Final model creation ----------------------------

batch_fin <- paste(temporal_block, '/Final_Models', sep = '')
mod_dir <- paste(temporal_block, '/Final_Models', sep = '')
rep_n <- 10
rep_type <- 'Bootstrap'
jackknife <- FALSE
out_format <- 'logistic'
project <- TRUE
G_var_dir <- paste(temporal_block, '/G_variables', sep = '')
ext_type <- 'ext'
write_mess <- FALSE
write_clamp <- FALSE
wait1 <- TRUE
run1 <- TRUE
args <- NULL

kuenm_mod(occ.joint = occ_joint, M.var.dir = M_var_dir, out.eval = out_eval, batch = batch_fin, 
          rep.n = rep_n, rep.type = rep_type, jackknife = jackknife, out.dir = mod_dir, 
          out.format = out_format, project = project, G.var.dir = G_var_dir, ext.type = ext_type, 
          write.mess = write_mess, write.clamp = write_clamp, maxent.path = maxent_path,args = args, 
          wait = wait1, run = run1)

#--------------------------------- Evaluation with independent data ------------------------------

occ_ind <- paste(temporal_block, '/G_galeus_indep.csv', sep = '')
replicates <- TRUE
out_feval <- paste(temporal_block, '/Final_Evaluation', sep = '')

fin_eval <- kuenm_feval(path = mod_dir, occ.joint = occ_joint, occ.ind = occ_ind, replicates = replicates,
                        out.eval = out_feval, threshold = threshold, rand.percent = rand_percent,
                        iterations = iterations, parallel.proc = paral_proc)

#------------------------ Uncertainty: median and range of best selected models ---------------------

sp_name <- 'Galeorhinus_galeus'
format <- 'asc'
project <- TRUE
stats <- c('median', 'range')
rep <- TRUE
scenarios <- 'current'
ext_type <- 'E'
out_dir <- paste(temporal_block, '/Final_Model_Stats', sep = '')

kuenm_modstats(sp.name = sp_name, fmod.dir = mod_dir, format = format, project = project,
               statistics = stats, replicated = rep, proj.scenarios = scenarios,
               ext.type = ext_type, out.dir = out_dir)


#--------------------------------------- Seasonal analysis --------------------------------------------

# Run each season at a time
i <- 1 # summer = 1, autumn = 2, winter = 3, spring = 4
temporal_block <- season[i]
dir.create(paste(temporal_block, sep = ''))

#------------------------------------------ G space --------------------------------------------

# Read predictors
env <- stack(paste(season[i], '.tif', sep = ''))
var_names <- c('Distance_to_coast', 'Depth', 'Surface_temperature', 'SST_fronts', 'Turbidity', 
               'Primary_productivity')
names(env) <- var_names

# Subset predictors relevant to G. galeus
var_names <- c('Distance_to_coast', 'Depth', 'Surface_temperature', 'SST_fronts', 'Turbidity', 
               'Primary_productivity')

# Rio de la Plata and Lagoa dos Patos not suitable
env <- mask(env, df_rdp, inverse = T)
env <- mask(env, df_lagoa, inverse = T)

#-------------------------------------- Data preparation ----------------------------------

# Initial M space for processing occurrences
dat_sub <- subset(dat, Season == season[i])
coord <- cbind(dat_sub$Longitude, dat_sub$Latitude)
occ.tot <- SpatialPoints(coord, proj4string = crs) # transform to spatial points
occ.buff <- buffer(occ.tot, width = 1000000) # calibration areas: 1000 km buffer around each point
env.M <- crop(env[[var_names]], occ.buff) # crop raster stack to buffer polygon
env.M <- mask(env.M, occ.buff) # mask raster stack to buffer polygon
env.M <- stack(env.M)

# Check correlation in M space
layerStats(env.M, 'pearson', na.rm = T) 
# Summer: depth and dist to coast > 0.8, keep depth; turbidity and primary productivity > 0.9, keep turbidity
# Autumn: depth and dist to coast > 0.8, keep depth; turbidity and primary productivity > 0.9, keep turbidity
# Winter: depth and dist to coast > 0.8, keep depth; turbidity and primary productivity > 0.9, keep turbidity
# Spring: depth and dist to coast > 0.8, keep depth; turbidity and primary productivity > 0.9, keep turbidity
var_set <- var_names[-c(1, 6)]
env.M <- env.M[[var_set]]

# Processing occurrences
# Save any points to nearest pixel with data that may be falling outside a layer of the M space raster 
library(rSDM) #package in GitHub (Pakillo/rSDM), install and call 'rSDM'
for(i in 1:length(env.M@layers)){
  spp <- SpatialPoints(dat_sub[, c('Longitude', 'Latitude')], crs)
  spp_corrected <- points2nearestcell(locs = spp, ras = env.M, layer = i, move = T) 
  dat_sub$Longitude <- round(spp_corrected@coords[, 1], 2)  # replace coordinates including those new if any
  dat_sub$Latitude <- round(spp_corrected@coords[, 2], 2)
} # plots will appear if any correction applies

# Independent data
occ_ind <- subset(dat_sub, Type == 'independent')
dups <- duplicated(occ_ind[c('Longitude', 'Latitude')])# Get rid of duplicates
occ_ind <- occ_ind[!dups, ]
occ_ind0 <- data.frame(species = Species, longitude = occ_ind$Longitude, latitude = occ_ind$Latitude, Type = 'independent')
write.csv(occ_ind0[, 1:3], paste(temporal_block, '/G_galeus_indep.csv', sep = ''), row.names = F)

# Calibration records
occ_cal <- subset(dat_sub, Type == 'calibration')
dups <- duplicated(occ_cal[c('Longitude', 'Latitude')])# Get rid of duplicates
occ_cal <- occ_cal[!dups, ]

if(temporal_block %in% c('autumn')){
  aaa <- subset(occ_cal, Latitude > -33)
  aaa <- aaa[sample(nrow(aaa), 10), ]
  occ_cal <- subset(occ_cal, Latitude < -33)
  occ_cal <- rbind(occ_cal, aaa)
}

# Environmental filtering
source('envSample.R')
coords <- occ_cal[, c('Longitude', 'Latitude')]
env.data <- extract(env.M, coords)
env.data <- as.data.frame(env.data)
coords <- envSample(coords, filters = list(env.data$Surface_temperature, env.data$Depth),
                    res = list(0.5, 10), do.plot = T)
occ_cal <- merge(occ_cal, coords, by.x = c('Longitude', 'Latitude'), by.y = c('lon', 'lat'))

# Detecting environmental outliers in species occurrences (Cobos et al., 2018)
# Getting data from the variables
variables_values <- na.omit(values(env.M)) # for the region of interest
occ_variables <- na.omit(cbind(as.numeric(row.names(occ_cal)), extract(env.M, occ_cal[, c('Longitude', 'Latitude')])))
occ_variables <- as.data.frame(occ_variables)
if(dim(variables_values)[1] > 10000) { # sample of 10000 values if more pixels exist
  variables_values <- variables_values[sample(1:dim(variables_values)[1], 10000), ]
}

# Plot for searching for potential environmental outliers
i <- 1
j <- 2
par(mfcol = c(1, 2), mar = c(4.5, 4, 0.5, 0.5), cex = 0.9)
plot(variables_values[, i], variables_values[, j], col = 'grey65', pch = 1,
     xlab = colnames(variables_values)[i], ylab = colnames(variables_values)[j])
points(occ_variables[, i + 1], occ_variables[, j + 1], col = 'blue', pch = 19, cex = 1.5)
legend('bottomright', legend = c('Region of interest', 'Occurrences'),
       pch = c(1, 19), col = c('gray65', 'black'), bty = 'n')
plot(variables_values[, i], variables_values[, j], col = 'gray65',
     pch = 1, xlab = colnames(variables_values)[i], ylab = colnames(variables_values)[j])
text(occ_variables[, i + 1], occ_variables[, j + 1], occ_variables[, 1], cex = 1, col = 'blue')
legend('bottomright', legend = 'Occurrence ID', bty = 'n')

# Remove outliers from calibration records
occ_variables <- subset(occ_variables, !V1 %in% c())
# Summer: c(118, 53, 28, 3, 2, 1, 4)
# Autumn: c(78, 85, 47, 62)
# Winter: c(63, 139, 159, 160)
# Spring: c(172)
occ_cal <- occ_cal[which(row.names(occ_cal) %in% as.character(occ_variables$V1)), ]

# Final corrected points (calibration + independent)
occ_corrected <- rbind(occ_cal, occ_ind)
write.csv(occ_corrected, paste(temporal_block, '/G_galeus_corrected.csv', sep = ''), row.names = F)

# Occurrence subsets for traning, testing and independent evaluation
# Create csv files for training and testing (75% and 25%)
occ_cal <- data.frame(species = occ_cal$Species, longitude = occ_cal$Longitude, latitude = occ_cal$Latitude)
occ_cal$check <- paste(occ_cal[, 'longitude'], occ_cal[, 'latitude'], sep = '_')
train <- occ_cal[sample(nrow(occ_cal), round((length(occ_cal[, 1]) / 10 * 7.5))), ]
test <- occ_cal[!occ_cal[, 4] %in% train[, 4], ]
occ_cal$check <- NULL; train$check <- NULL; test$check <- NULL
write.csv(occ_cal[, 1:3], paste(temporal_block, '/G_galeus_joint.csv', sep = ''), row.names = F)
write.csv(train[, 1:3], paste(temporal_block, '/G_galeus_train.csv', sep = ''), row.names = F)
write.csv(test[, 1:3], paste(temporal_block, '/G_galeus_test.csv', sep = ''), row.names = F)

# Final M space (calibration areas)
coord <- cbind(occ_cal$longitude, occ_cal$latitude)
occ.tot <- SpatialPoints(coord, proj4string = crs) 
occ.buff <- buffer(occ.tot, width = 1000000) 
env.M <- crop(env[[var_set]], extent(occ.buff) + 1)
env.M <- mask(env.M, occ.buff) 
env.M <- mask(env.M, df_rdp, inverse = T)
env.M <- mask(env.M, df_lagoa, inverse = T)
env.M <- stack(env.M)

# Create M space for modelling
dir.create(paste(temporal_block, '/M_variables', sep = ''))
dir.create(paste(temporal_block, '/M_variables/Set1', sep = ''))
writeRaster(env.M[[var_set]], filename = paste(temporal_block, '/M_variables/Set1/env.asc', sep = ''),
            format = 'ascii', bylayer = T, suffix = var_set)
dir.create(paste(temporal_block, '/G_variables', sep = ''))
dir.create(paste(temporal_block, '/G_variables/Set1', sep = ''))
dir.create(paste(temporal_block, '/G_variables/Set1/current', sep = ''))
writeRaster(env[[var_set]], filename = paste(temporal_block, '/G_variables/Set1/current/env.asc', sep = ''),
            format = 'ascii', bylayer = T, suffix = var_set)

#------------------------------------------ Modelling -------------------------------------------

# Seed
set.seed(111)

#--------------------------------------- Candidate models --------------------------------------------

occ_joint <- paste(temporal_block, '/G_galeus_joint.csv', sep = '')
occ_tra <- paste(temporal_block, '/G_galeus_train.csv', sep = '')
M_var_dir <- paste(temporal_block, '/M_variables', sep = '')
batch_cal <- paste(temporal_block, '/Candidate_Models', sep = '')
out_dir <- paste(temporal_block, '/Candidate_Models', sep = '')
reg_mult <- c(seq(0.1, 2, 0.1), 2.5, 3, 4, 5, 7, 10) # finer set of values
f_clas <- c('lq', 'lp', 'lqp')
args <- NULL
maxent_path <- 'MAXENT JAR DIRECTORY'
wait <- TRUE
run <- TRUE

kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal, 
          out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, args = args, maxent.path = maxent_path, 
          wait = wait, run = run)

#----------------------------- Evaluation and selection of best models ------------------------

occ_test <- paste(temporal_block, '/G_galeus_test.csv', sep = '')
out_eval <- paste(temporal_block, '/Calibration_Results', sep = '')
threshold <- 5
rand_percent <- 50
iterations <- 500
kept <- FALSE
selection <- 'OR_AICc' # OR_AICc
paral_proc <- FALSE 

kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test, 
            batch = batch_cal, out.eval = out_eval, threshold = threshold, rand.percent = rand_percent, 
            iterations = iterations, kept = kept, selection = selection, parallel.proc = paral_proc)

#-------------------------------------- Final model creation ----------------------------

batch_fin <- paste(temporal_block, '/Final_Models', sep = '')
mod_dir <- paste(temporal_block, '/Final_Models', sep = '')
rep_n <- 10
rep_type <- 'Bootstrap'
jackknife <- FALSE
out_format <- 'logistic'
project <- TRUE
G_var_dir <- paste(temporal_block, '/G_variables', sep = '')
ext_type <- 'ext'
write_mess <- FALSE
write_clamp <- FALSE
wait1 <- TRUE
run1 <- TRUE
args <- NULL

kuenm_mod(occ.joint = occ_joint, M.var.dir = M_var_dir, out.eval = out_eval, batch = batch_fin, 
          rep.n = rep_n, rep.type = rep_type, jackknife = jackknife, out.dir = mod_dir, 
          out.format = out_format, project = project, G.var.dir = G_var_dir, ext.type = ext_type, 
          write.mess = write_mess, write.clamp = write_clamp, maxent.path = maxent_path,args = args, 
          wait = wait1, run = run1)

#--------------------------------- Evaluation with independent data ------------------------------

occ_ind <- paste(temporal_block, '/G_galeus_indep.csv', sep = '')
replicates <- TRUE
out_feval <- paste(temporal_block, '/Final_Evaluation', sep = '')

fin_eval <- kuenm_feval(path = mod_dir, occ.joint = occ_joint, occ.ind = occ_ind, replicates = replicates,
                        out.eval = out_feval, threshold = threshold, rand.percent = rand_percent,
                        iterations = iterations, parallel.proc = paral_proc)

#------------------------ Uncertainty: median and range of best selected models ---------------------

sp_name <- 'Galeorhinus_galeus'
format <- 'asc'
project <- TRUE
stats <- c('median', 'range')
rep <- TRUE
scenarios <- 'current'
ext_type <- 'E'
out_dir <- paste(temporal_block, '/Final_Model_Stats', sep = '')

kuenm_modstats(sp.name = sp_name, fmod.dir = mod_dir, format = format, project = project,
               statistics = stats, replicated = rep, proj.scenarios = scenarios,
               ext.type = ext_type, out.dir = out_dir)

#--------------------------------------------- END -------------------------------------
