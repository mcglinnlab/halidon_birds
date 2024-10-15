# Vignette Set Up
library(mobr)

# remember to add species filtering step in data processing
# check on which species has 31 and 12 for point count abundance

dat <- read.csv('./data/bird_data_year4.csv')
# small data cleaning
dat <- subset(dat, subset = dat$species != "")
dat <- subset( dat, subset = dat$prelim_data != 1)
dat$X25m <- ifelse(is.na(dat$X25m), 0, dat$X25m)
dat$X50m <- ifelse(is.na(dat$X50m), 0, dat$X50m)
dat$date <- as.Date(dat$date, format = "%m/%d/%Y")
head(dat)

# create unique sampling event id
dat$uni_id_date <- with(dat, paste(new.site.id, date, sep='_'))

# site by species index - total
comm <- with(dat, tapply(total, list(uni_id_date, species), sum))
comm <- ifelse(is.na(comm), 0, comm)
summary(comm)
table(dat$species)

######## mobr Beta Diversity Vignette START ######## 
# 1:_ where _ is the number of sites that it is comparing between
  # Do we need to have the specific number of sites? Will mobr automatically separate the values by site?
  # Or, is this better to have by the uni_id_date that was created above (of which we have 112)?
# Determining how many sites I have
library(dplyr)
length(unique(dat$new.site.id))
# I have 38 sites
# Calculate whitaker's beta
calc_comm_div(comm[1:38, ], 'S')
# Calculate beta for S_PIE
calc_comm_div(comm[1:38, ], 'S_PIE')
# Calculate beta for C (specific coverage)
calc_comm_div(comm[1:38, ], 'S_C')
# Calculate beta for S_n (rarefied richness)
calc_comm_div(comm[1:38, ], 'S_n', effort = 5)
calc_comm_div(comm[1:38, ], 'S_n', effort = 20)
# Calculate just beta diversity - not generally recommended without alpha and gamma
calc_beta_div(comm[1:38, ], 'S')

######## mobr Scale-Dependent Biodiversity Changes Vignette START ######## 
library(mobr)
library(dplyr)
library(ggplot2)
library(googlesheets4)

# Data for this vignette
dat <- read.csv('./data/bird_data_year4.csv')
# Remove wetland 43 because it was mulched in 2024
dat <- subset(dat, !(wetland_id %in% c('HH-43-W', 'HH-43-U')))
# remove wetland 5 because it was a control disturbed (only one replicate of this type)
dat <- subset(dat, !(wetland_id %in% c('HH-05-W', 'HH-05-U')))
# Removes empty identifications and preliminary data
dat <- subset(dat, subset = dat$species != "")
dat <- subset( dat, subset = dat$prelim_data != 1)
dat$X25m <- ifelse(is.na(dat$X25m), 0, dat$X25m)
dat$X50m <- ifelse(is.na(dat$X50m), 0, dat$X50m)
dat$date <- as.Date(dat$date, format = "%m/%d/%Y")
head(dat)
dat$uni_id_date <- with(dat, paste(new.site.id, date, sep='_'))

# check that total in dat is actual total
dat[(dat$X25m + dat$X50m) != dat$total, ]

# Create Community Matrix
# going to combine 25 and 50 m radius by using total (potentially change in future)
comm <- with(dat, tapply(total, list(wetland_id, species), sum))
comm <- ifelse(is.na(comm), 0, comm)
summary(comm)
table(dat$species)

dim(comm)
comm[1:5, 1:5]
comm_sites <- substr(row.names(comm),1,7)



# Site Attribute Table
hh_att <- read_sheet("https://docs.google.com/spreadsheets/d/1LsEtgikV3PtJJPZsW5IFkuCzZSZTZVg3Y2vQdgvS3V8/edit?usp=sharing", 
                         sheet = "attributes")
# check that old sites in hh_att match comm site id
# reorder the attribute table to match the order of the sites in the 
# community matrix
hh_att <- hh_att[match(row.names(comm), hh_att$site), ]


# Work Through
str(comm)
head(hh_att)

# Data Prep
hh_mob_in <- make_mob_in(comm, hh_att, coord_names = c('utm_easting', 'utm_northing'))
hh_mob_in

# Data Analysis
plot_rarefaction(hh_mob_in, 'treatment', ref_level = 'control-closed', 'sSBR', lwd = 4)

oldpar <- par(no.readonly = TRUE)
par(mfrow = c(1,2))
plot_rarefaction(hh_mob_in, 'treatment', 'cut-leave', 'IBR', 
                 leg_loc = 'bottomright')

par(oldpar)

oldpar <- par(no.readonly = TRUE)
par(mfrow = c(1, 2))
plot_abu(hh_mob_in, 'treatment', type = 'rad', scale = 'alpha', log = 'x')
plot_abu(hh_mob_in, 'treatment', type = 'rad', scale = 'gamma' , log = 'x')

par(oldpar)

# Two-Scale Analysis
indices <- c('N', 'S', 'S_n', 'S_PIE')
hh_div <- tibble(comm) %>%
  group_by(group = hh_att$treatment) %>%
  group_modify(~ calc_comm_div(.x, index = indices, effort = 5,
                               extrapolate = TRUE))
  # Having same issue as above with different plot values
head(hh_div)

# Plot Individual Diversity Metrics
plot_comm_div(hh_div, 'S')

plot_comm_div(hh_div, 'N')

plot_comm_div(hh_div, 'S_n')

plot_comm_div(hh_div, 'S_PIE')

# Compute beta diversity directly
calc_beta_div(comm, c('S', 'S_n', 'S_PIE'), effort = 5)

# Plots all the diversity metrics 
plot_comm_div(hh_div)

# Multi-Scale Analysis - examines the difference between sSBR, nsSBR, and IBR rarefaction 
# curves to tell how treatment influences ruchness on different community structure components
hh_mob_in <- make_mob_in(comm, hh_att,
                          coord_names = c('utm_easting', 'utm_northing'))
hh_deltaS <- get_delta_stats(hh_mob_in, env_var = 'treatment', ref_level='control-closed',
                              type='discrete', log_scale=TRUE, n_perm = 199)

# Make rarefaction curves to look at
plot(hh_deltaS, stat = 'b1', scale_by = 'indiv', display='S ~ effort')

# Considers effect sizes as a function of scale
plot(hh_deltaS, stat = 'b1', scale_by = 'indiv', display='stat ~ effort')


hh_mob_in2 <- subset(hh_mob_in, treatment %in% c('control-closed', 'cut-remove'))
hh_deltaS2 <- get_delta_stats(hh_mob_in2, env_var = 'treatment', ref_level='control-closed',
                             type='discrete', log_scale=TRUE, n_perm = 199)
plot(hh_deltaS2, stat = 'b1', scale_by = 'indiv', display='stat ~ effort')




