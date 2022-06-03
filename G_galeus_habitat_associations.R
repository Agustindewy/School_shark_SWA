
# De Wysiecki et al. 
# "Population-scale habitat use of school sharks (Triakidae: Galeorhinus galeus) in the Southwest Atlantic: 
# insights from temporally explicit niche modelling and habitat associations"

# Habitat associations

library(raster)

setwd('SET YOUR WORKING DIRECTORY')

# Season
season <- c('summer', 'autumn', 'winter', 'spring')

#useful vectors
Vars <- c('Latitude', 'Bathymetry', 'Surface_temperature', 'SST_fronts', 'Turbidity')
var_names <- c('Bathymetry', 'Surface_temperature', 'SST_fronts', 'Turbidity')

# Functions
source('Perry_&_Smith_(1994).R') # Perry & Smith custom function. Author Federico Cortés (cortes.federico@gmail.com)
my.quantile <- function(X) { quantile(X, probs = c(0.025, 0.5, 0.975)) }
closest <- function(xv, sv) { xv[which(abs(xv - sv) == min(abs(xv - sv)))] } # which xv value is the closest to a sv value
rangeAMB.ft <- function(X) { with(perry_sex.stage, Xambiente[which(Ft == closest(Ft, X)[1])[1]]) } # function to compute variable ranges of f(t)
rangeAMB.gt <- function(X) { with(perry_sex.stage, Xambiente[which(Gt == closest(Gt, X)[1])[1]]) } # function to compute variable ranges of g(t)
rangeAMB.ft_females <- function(X) {with(perry_females, Xambiente[which(Ft == closest(Ft, X)[1])[1]]) }
rangeAMB.gt_females <- function(X) {with(perry_females, Xambiente[which(Gt == closest(Gt, X)[1])[1]]) }
rangeAMB.ft_males <- function(X) {with(perry_males, Xambiente[which(Ft == closest(Ft, X)[1])[1]]) }
rangeAMB.gt_males <- function(X) {with(perry_males, Xambiente[which(Gt == closest(Gt, X)[1])[1]]) }
rangeAMB.ft_adults <- function(X) {with(perry_adults, Xambiente[which(Ft == closest(Ft, X)[1])[1]]) }
rangeAMB.gt_adults <- function(X) {with(perry_adults, Xambiente[which(Gt == closest(Gt, X)[1])[1]]) }
rangeAMB.ft_juveniles <- function(X) {with(perry_juveniles, Xambiente[which(Ft == closest(Ft, X)[1])[1]]) }
rangeAMB.gt_juveniles <- function(X) {with(perry_juveniles, Xambiente[which(Gt == closest(Gt, X)[1])[1]]) }
detachAllPackages <- function() { # detaching function to avoid conflicts between raster and tydiverse packages

  basic.packages.blank <- c(    
    'stats',    
    'graphics',    
    'grDevices',    
    'utils',   
    'datasets',  
    'methods',    
    'base'    
  )    
  basic.packages <- paste('package:', basic.packages.blank, sep = '')   
  package.list <- search()[ifelse(unlist(gregexpr('package:', search())) == 1, TRUE, FALSE)]   
  package.list <- setdiff(package.list, basic.packages)   
  if (length(package.list) > 0) {   
    for (package in package.list) {   
      detach(package, character.only = TRUE)   
    }   
  }    
}

#------------------------------------- Prepare data -------------------------------------

# Read occurrences
dat0 <- read.csv('Presences.csv')
dat0$Stage <- as.factor(dat0$Stage)

# Associate seasonal environmental values
dat <- data.frame()
for(i in season) {
  
  subdat0 <- subset(dat0, Season == i)
  env <- stack(paste(i, '.tif', sep = ''))[[2:5]]
  names(env) <- var_names
  
  # extract and clean
  subdat <- extract(env, subdat0[, c('Longitude', 'Latitude')])
  subdat <- cbind(subdat0, subdat)
  subdat <- subdat[which(complete.cases(subdat[, var_names])), ]
  
  dat <- rbind(dat, subdat)
}

# Latitude and bathymetry as positive
dat$Latitude <- -dat$Latitude
dat$Bathymetry <- -dat$Bathymetry

# create dummy variables
# binary column: 1 if female, 0 if male, -1 if no sex
dat$female <- ifelse(is.na(dat$Sex), -1, ifelse(dat$Sex == 'female', 1, 0))
# binary column: 1 if male, 0 if female, -1 if no Sex
dat$male <- ifelse(is.na(dat$Sex), -1, ifelse(dat$Sex == 'male', 1, 0))
# binary column: 1 if adult, 0 if juvenile or neonate, -1 if no Stage
dat$adult <- ifelse(is.na(dat$Stage), -1, ifelse(dat$Stage == 'adult', 1, 0))
# binary column: 1 if juvenile, 0 if adult or neonate, -1 if no Stage
dat$juvenile <- ifelse(is.na(dat$Stage), -1, ifelse(dat$Stage == 'juvenile', 1, 0))
# binary column: 1 if neonate, 0 if adult or juvenile, -1 if no Stage
dat$neonate <- ifelse(is.na(dat$Stage), -1, ifelse(dat$Stage == 'neonate', 1, 0))

# Collapse to unique observation points and seasons to clean true duplicates
dat_final <- aggregate(cbind(female, male, adult, juvenile, neonate) 
                    ~ round(Latitude, 4) + round(Longitude, 4) + Season + Bathymetry + Surface_temperature + SST_fronts + Turbidity,
                    data = dat, FUN = max, na.action = na.pass)

# Change Latitude and Longitude names
names(dat_final)[names(dat_final) == 'round(Latitude, 4)'] <- 'Latitude'
names(dat_final)[names(dat_final) == 'round(Longitude, 4)'] <- 'Longitude'

# Create dummy variables to get g(t)
# binary column: 1 if sex, 0 if no sex
dat_final$Sex <- ifelse(dat_final$female > 0 | dat_final$male > 0, 1, 0)
# binary column: 1 if stage, 0 if no stage
dat_final$Stage <- ifelse(dat_final$adult > 0 | dat_final$juvenile > 0 | dat_final$neonate > 0, 1, 0)
# binary column: 1 if sex and stage, 0 if otherwise 
dat_final$sex.stage <- ifelse(dat_final$Sex > 0 & dat_final$Stage > 0, 1, 0) # g(t)


#----------------------------------------- Aim 1 -------------------------------------

# Explore bias in data with biological information

list_final1 <- list() # empty lists
N.iter = 1000 # number of iterations

# Perry & Smith analysis
for(j in Vars) {

  dat_final0 <- dat_final
  
  # clean outliers
  if(j %in% Vars) {
    lower_bound <- quantile(dat_final0[, j], 0.01)
    upper_bound <- quantile(dat_final0[, j], 0.99)
    outlier_ind <- which(dat_final0[, j] > lower_bound&dat_final0[, j] < upper_bound)
    dat_final0 <- dat_final0[outlier_ind, ]
  }

  # define var range interval, min and max var values, and N iterations for bootstrapping
  min_var <- min(dat_final0[, j])
  max_var <- max(dat_final0[, j])
  rango_var <- (max_var - min_var) / 500
  
  # Perry & Smith analysis, f(t), g(t) and Dmax calculation for presences with biological data
  perry_sex.stage <- perry_D(AMB = dat_final0[, j], CAP = dat_final0$sex.stage, rango_t = rango_var, min_t = min_var, max_t = max_var)

  # Perry & Smith analysis, bootstrapping of f(t) and g(t) for presences with biological data
  boot_sex.stage <- perry_BOOT(AMB = dat_final0[, j], CAP = dat_final0$sex.stage, rango_t = rango_var, min_t = min_var, max_t = max_var, N_iter = N.iter)

  # 95% environmental range of f(t) and g(t)
  rangeAMB <- matrix(c(0.025, 0.5, 0.975), ncol = 1)
  rangeAMB <- data.frame(percentil = rangeAMB, ranFt = apply(rangeAMB, 1, rangeAMB.ft),
                         ranGt = apply(rangeAMB, 1, rangeAMB.gt))

  # proportion of environment covered by 95% of presences with biological information
  perAMBcap <- with(perry_sex.stage, Ft[which(Gt == closest(Gt, 0.975)[1])[1]] -
                   ifelse(Ft[1] > 0.025 & Gt[1] > 0.025, 0, Ft[which(Gt == closest(Gt, 0.025)[1])[1]]))
  
  # results
  list_final1[[j]] <- c(perry_sex.stage, boot_sex.stage[c('p', 'distD', 'randomFt', 'randomGt')], rangeAMB, prop = perAMBcap)
}

# Plot
detachAllPackages()
library(tidyverse)
library(viridis)

# Prepare data
final_df <- data.frame()

for(j in Vars) {
  
  dtaBOOT_sex.stage <- rbind(data.frame(VAR = j, category = 'sex.stage', curva = 'Ft', cdfOBS = list_final1[[j]]$Ft, Xamb = list_final1[[j]]$Xambiente,
                                     as.data.frame(t(apply(list_final1[[j]]$randomFt, 1, my.quantile)))),
                          data.frame(VAR = j, category = 'sex.stage', curva = 'Gt', cdfOBS = list_final1[[j]]$Gt, Xamb = list_final1[[j]]$Xambiente,
                                     as.data.frame(t(apply(list_final1[[j]]$randomGt, 1, my.quantile)))))
  names(dtaBOOT_sex.stage)[6:8] <- c('ICinf', 'median', 'ICsup')
  
  final_df <- rbind(final_df, dtaBOOT_sex.stage)
}

# Get p-values
test.df <-data.frame()
for(j in Vars) {

  test.BOOT <- data.frame(VAR = j, category = 'sex.stage', D = round(list_final1[[j]]$Dmax, 3), p_value = round(list_final1[[j]]$p, 4),
                       AMBmin = round(list_final1[[j]]$ranFt[[1]], 2), AMBmax = round(list_final1[[j]]$ranFt[[3]], 2),
                       CAPmin = round(list_final1[[j]]$ranGt[[1]], 2), CAPmax = round(list_final1[[j]]$ranGt[[3]], 2),
                       perAMBcap = round(list_final1[[j]]$prop, 2))
  test.df <- rbind(test.df, test.BOOT)
}

# Final plot, g(t) line, f(t) line plus bootstrapping 95% confidence interval of g(t)
final_df <- transform(final_df, VAR = factor(VAR, levels = Vars)) # ordering of factors

# Labeller
var_names <- as_labeller(c(Latitude = 'Latitude~(ºS)', Bathymetry = 'Bathymetry~(m)', Turbidity = 'Coefficient~Kd490~(m^-1)',
                           Surface_temperature = 'Surface~Temperature~(ºC)', SST_fronts = 'Surface~Temperature~fronts~(ºC)'), default = label_parsed)

cols <- viridis_pal(option = 'E')(2)

ggplot(final_df, aes(x = Xamb)) + ylab('Cumulative frequency') + xlab('Environmental variable') +
  geom_ribbon(aes(ymin = ICinf, ymax = ICsup), data = subset(final_df, curva == 'Gt'), alpha = 0.35, fill = cols[1]) +
  geom_line(aes(y = cdfOBS), data = subset(final_df, curva == 'Gt'), size = 0.75, col = cols[2]) +
  geom_line(aes(y = cdfOBS), data = subset(final_df, curva == 'Ft'), size = 0.5, col = cols[1], linetype = 'dashed') +
  theme(panel.grid.minor = element_blank(),
      panel.grid.major = element_line(size = 0.3)) +
  facet_wrap(~ VAR, scales = 'free_x', ncol = 3, labeller = var_names) 
ggsave('Figure 3.pdf', width = 19, height = 12, units = 'cm')

write.csv(test.df, 'Table_aim_1.csv', row.names = F)


#----------------------------------------- Aim 2 -------------------------------------

# Intraspecific habitat associations

# Use only presences with sex and/or stage information
dat_final_sex <- subset(dat_final, Sex > 0)
dat_final_stage <- subset(dat_final, Stage > 0)

df_final <- data.frame() # empty data frame
N.iter <- 1000 # number of iterations

# Perry & Smith analysis
for(k in season) {

  pre_df_final <- data.frame()
  
  # Perry & Smith analysis
  for(j in Vars){

    # subset by season
    dat_final_sex0 <- subset(dat_final_sex, Season == k)
    dat_final_stage0 <- subset(dat_final_stage, Season == k)
    
    # clean outliers
    if(j %in% Vars[-6]) {
      
      # sexes
      lower_bound <- quantile(dat_final_sex0[, j], 0.01)
      upper_bound <- quantile(dat_final_sex0[, j], 0.99)
      outlier_ind <- which(dat_final_sex0[, j] > lower_bound&dat_final_sex0[, j] < upper_bound)
      dat_final_sex0 <- dat_final_sex0[outlier_ind, ]
      
      # stages
      lower_bound <- quantile(dat_final_stage0[, j], 0.01)
      upper_bound <- quantile(dat_final_stage0[, j], 0.99)
      outlier_ind <- which(dat_final_stage0[, j] > lower_bound&dat_final_stage0[, j] < upper_bound)
      dat_final_stage0 <- dat_final_stage0[outlier_ind, ]
    }
    
    # define var range interval, min and max var values, and N iterations for bootstrapping
    # sexes
    min_var_sexes <- min(dat_final_sex0[, j])
    max_var_sexes <- max(dat_final_sex0[, j])
    rango_var_sexes <- (max_var_sexes - min_var_sexes) / 500
    
    # stages
    min_var_stages <- min(dat_final_stage0[, j])
    max_var_stages <- max(dat_final_stage0[, j])
    rango_var_stages <- (max_var_stages - min_var_stages) / 500
    
    # P&S analysis, f(t), g(t) and Dmax calculation for sexes
    perry_females <- perry_D(AMB = dat_final_sex0[, j], CAP = dat_final_sex0$female, rango_t = rango_var_sexes, min_t = min_var_sexes, max_t = max_var_sexes)
    perry_males <- perry_D(AMB = dat_final_sex0[, j], CAP = dat_final_sex0$male, rango_t = rango_var_sexes, min_t = min_var_sexes, max_t = max_var_sexes)
    
    # P&S analysis, f(t), g(t) and Dmax calculation for stages
    perry_adults <- perry_D(AMB = dat_final_stage0[, j], CAP = dat_final_stage0$adult, rango_t = rango_var_stages, min_t = min_var_stages, max_t = max_var_stages)
    perry_juveniles <- perry_D(AMB = dat_final_stage0[, j], CAP = dat_final_stage0$juvenile, rango_t = rango_var_stages, min_t = min_var_stages, max_t = max_var_stages)
    
    # P&S analysis, bootstrapping of f(t) and g(t) for sexes
    boot_females <- perry_BOOT(AMB = dat_final_sex0[, j], CAP = dat_final_sex0$female, rango_t = rango_var_sexes, min_t = min_var_sexes, max_t = max_var_sexes, N_iter = N.iter)
    boot_males <- perry_BOOT(AMB = dat_final_sex0[, j], CAP = dat_final_sex0$male, rango_t = rango_var_sexes, min_t = min_var_sexes, max_t = max_var_sexes, N_iter = N.iter)
    
    # P&S analysis, bootstrapping of f(t) and g(t) for stages
    boot_adults <- perry_BOOT(AMB = dat_final_stage0[, j], CAP = dat_final_stage0$adult, rango_t = rango_var_stages, min_t = min_var_stages, max_t = max_var_stages, N_iter = N.iter)
    boot_juveniles <- perry_BOOT(AMB = dat_final_stage0[, j], CAP = dat_final_stage0$juvenile, rango_t = rango_var_stages, min_t = min_var_stages, max_t = max_var_stages, N_iter = N.iter)
    
    # 95% environmental range of f(t) and g(t)
    rangeAMB <- matrix(c(0.025, 0.5, 0.975), ncol = 1)
    rangeAMB_females <- data.frame(percentil = rangeAMB, ranFt = apply(rangeAMB, 1, rangeAMB.ft_females), ranGt = apply(rangeAMB, 1, rangeAMB.gt_females))
    rangeAMB_males <- data.frame(percentil = rangeAMB, ranFt = apply(rangeAMB, 1, rangeAMB.ft_males), ranGt = apply(rangeAMB, 1, rangeAMB.gt_males))
    rangeAMB_adults <- data.frame(percentil = rangeAMB, ranFt = apply(rangeAMB, 1, rangeAMB.ft_adults), ranGt = apply(rangeAMB, 1, rangeAMB.gt_adults))
    rangeAMB_juveniles <- data.frame(percentil = rangeAMB, ranFt = apply(rangeAMB, 1, rangeAMB.ft_juveniles), ranGt = apply(rangeAMB, 1, rangeAMB.gt_juveniles))
    
    # differences between female g(t) and male g(t) based on Sagarese et al. (2014)
    # maximum distance between female boot g(t) and male boot g(t)
    DmaxDIST_f_m <- apply(abs(boot_females$randomGt - boot_males$randomGt), 2, max)
    # maximum observed distance between female g(t) and male g(t)
    DmaxOBS_f_m <- max(abs(perry_females$Gt - perry_males$Gt))
    # p-value, how many boots Dmax are >= to Dmax observed
    p.value_f_m <- length(which(DmaxDIST_f_m >= DmaxOBS_f_m)) / N.iter
    
    # differences between adult g(t) and juvenile g(t) based on Sagarese et al. (2014)
    # maximum distance between adult boot g(t) and juvenile boot g(t)
    DmaxDIST_a_j <- apply(abs(boot_adults$randomGt - boot_juveniles$randomGt), 2, max)    
    # maximum observed distance between adult g(t) and juvenile g(t)
    DmaxOBS_a_j <- max(abs(perry_adults$Gt - perry_juveniles$Gt))
    # p-value, how many boots Dmax are >= to Dmax observed
    p.value_a_j <- length(which(DmaxDIST_a_j >= DmaxOBS_a_j)) / N.iter
    
    # create data frame for all categories
    pre_df = rbind(data.frame(Season = k, Var = j, category = 'female', Xambiente = perry_females$Xambiente,
                            Ft = perry_females$Ft, Gt = perry_females$Gt, Dmax = perry_females$Dmax, p.value = boot_females$p,
                            rangeAMBmin = rangeAMB_females$ranFt[1], rangeAMBmax = rangeAMB_females$ranFt[3],
                            rangeCAPmin = rangeAMB_females$ranGt[1], rangeCAPmax = rangeAMB_females$ranGt[3],
                            Dmax2 = DmaxOBS_f_m,p.value2 = p.value_f_m,randomGt = boot_females$randomGt),
                 data.frame(Season = k, Var = j, category = 'male', Xambiente = perry_males$Xambiente, 
                            Ft = perry_males$Ft, Gt = perry_males$Gt, Dmax = perry_males$Dmax, p.value = boot_males$p, 
                            rangeAMBmin = rangeAMB_males$ranFt[1], rangeAMBmax = rangeAMB_males$ranFt[3], 
                            rangeCAPmin = rangeAMB_males$ranGt[1], rangeCAPmax = rangeAMB_males$ranGt[3], 
                            Dmax2 = DmaxOBS_f_m, p.value2 = p.value_f_m, randomGt = boot_males$randomGt), 
                 data.frame(Season = k, Var = j, category = 'adult', Xambiente = perry_adults$Xambiente, 
                            Ft = perry_adults$Ft, Gt = perry_adults$Gt, Dmax = perry_adults$Dmax, p.value = boot_adults$p, 
                            rangeAMBmin = rangeAMB_adults$ranFt[1], rangeAMBmax = rangeAMB_adults$ranFt[3], 
                            rangeCAPmin = rangeAMB_adults$ranGt[1], rangeCAPmax = rangeAMB_adults$ranGt[3], 
                            Dmax2 = DmaxOBS_a_j, p.value2 = p.value_a_j, randomGt = boot_adults$randomGt), 
                 data.frame(Season = k, Var = j, category = 'juvenile', Xambiente = perry_juveniles$Xambiente, 
                            Ft = perry_juveniles$Ft, Gt = perry_juveniles$Gt, Dmax = perry_juveniles$Dmax, p.value = boot_juveniles$p, 
                            rangeAMBmin = rangeAMB_juveniles$ranFt[1], rangeAMBmax = rangeAMB_juveniles$ranFt[3], 
                            rangeCAPmin = rangeAMB_juveniles$ranGt[1], rangeCAPmax = rangeAMB_juveniles$ranGt[3], 
                            Dmax2 = DmaxOBS_a_j, p.value2 = p.value_a_j, randomGt = boot_juveniles$randomGt))
                 
    pre_df_final <- rbind(pre_df_final, pre_df)
    
  }
  
  df_final <- rbind(df_final, pre_df_final)
  
}

# Plot

detachAllPackages()
library(tidyverse)
library(viridis)

# Prepare data
dtaBOOT_all <- cbind(df_final[, 1:14], as.data.frame(t(apply(df_final[, -c(1:14)], 1, my.quantile))))
names(dtaBOOT_all)[15:17] <- c('ICinf', 'median', 'ICsup')

# Get p-values
test.df <- aggregate(cbind(Dmax,  p.value,  rangeAMBmin,  rangeAMBmax, rangeCAPmin, rangeCAPmax, Dmax2, p.value2) ~ 
                    Season + Var + category, data = dtaBOOT_all, FUN = mean)
test.df$Dmax <- round(test.df$Dmax, 3)
test.df$Dmax2 <- round(test.df$Dmax2, 3)
test.df$rangeAMBmin <- round(test.df$rangeAMBmin, 2)
test.df$rangeAMBmax <- round(test.df$rangeAMBmax, 2)
test.df$rangeCAPmin <- round(test.df$rangeCAPmin, 2)
test.df$rangeCAPmax <- round(test.df$rangeCAPmax, 2)

# Final plot, g(t) lines and f(t) lines
final_df <- transform(dtaBOOT_all, Var = factor(Var, levels = Vars), 
                   category = factor(category, levels = c('female', 'male', 'adult', 'juvenile')), 
                   Season = factor(Season, levels = season)) # ordering of factors

sub_sexes <- subset(final_df, category %in% c('female', 'male'))
sub_stages <- subset(final_df, category %in% c('adult', 'juvenile'))

cols <- c('#00204DFF', '#FFEA46FF')

# Labeller
var_names <- as_labeller(c(Latitude = 'Latitude~(ºS)', Bathymetry = 'Bathymetry~(m)', Turbidity = 'Kd490~(m^-1)',
                           Surface_temperature = 'SST~(ºC)', SST_fronts = 'SST~fronts~(ºC)', 
                           summer = 'Summer', autumn = 'Autumn', winter = 'Winter', spring = 'Spring'), default = label_parsed)

# Sexes
ggplot(sub_sexes, aes(x = Xambiente, color = category)) + ylab('Cumulative frequency') + xlab('Environmental variable') +
  geom_ribbon(aes(ymin = ICinf, ymax = ICsup, fill = category), alpha = 0.35, col = NA) +
  geom_line(aes(y = Gt, col = category), size = 0.6) +
  geom_line(aes(y = Ft), size = 0.6, col = 'black', linetype = 'dashed') +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme(strip.background = element_rect(fill = '#7C7B78FF'),
        panel.grid.minor = element_blank(), legend.position = 'none',
        panel.grid.major = element_line(size = 0.3),
        strip.text = element_text(size = 10)) +
  facet_grid(c('Season', 'Var'), labeller = var_names, scales = 'free') 
ggsave('Figure S7.pdf', dpi = 900, width = 19, height = 13, units = 'cm')

# Stages
ggplot(sub_stages, aes(x = Xambiente, color = category)) + ylab('Cumulative frequency') + xlab('Environmental variable') +
  geom_ribbon(aes(ymin = ICinf, ymax = ICsup, fill = category), alpha = 0.35, col = NA) +
  geom_line(aes(y = Gt, col = category), size = 0.6) +
  geom_line(aes(y = Ft), size = 0.6, col = 'black', linetype = 'dashed') +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme(strip.background = element_rect(fill = '#7C7B78FF'),
        panel.grid.minor = element_blank(), legend.position = 'none',
        panel.grid.major = element_line(size = 0.3),
        strip.text = element_text(size = 10)) +
  facet_grid(c('Season', 'Var'), labeller = var_names, scales = 'free') 
ggsave('Figure S8.pdf', dpi = 900, width = 19, height = 13, units = 'cm')

write.csv(test.df, 'Table_aim_2.csv', row.names = F)


#--------------------------------------------- END -------------------------------------
