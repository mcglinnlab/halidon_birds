
######## mobr Scale-Dependent Biodiversity Changes Vignette START ######## 
library(mobr)
library(googlesheets4)

# Data for this vignette
dat <- read.csv('./data/bird_data_year4.csv')

sp_to_drop <- c("AMCR", "BBWD", "WOST", "NOBO", "CAGO", 
                "CONI", "FICR", "GREG", "MIKI", "RSHA",
                "RTHA", "LAGU")
dat <- subset(dat, !(species %in% sp_to_drop))

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
# there is a lot where total is incorrect
# recaculate total
dat$total <- dat$X25m + dat$X50m


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
# add recoding of treatments
hh_att$treatment2 <- hh_att$treatment
hh_att$treatment2 <- ifelse(hh_att$treatment2 == 'hack-squirt', 'control-closed',
                            hh_att$treatment2)
#hh_att$treatment2 <- ifelse(hh_att$treatment2 == 'cut-leave', 'control-closed',
#                            hh_att$treatment2)

# Work Through
str(comm)
head(hh_att)

# Data Prep
hh_mob_in <- make_mob_in(comm, hh_att, coord_names = c('utm_easting', 'utm_northing'))
hh_mob_in

# Data Analysis
par(mfrow=c(1,1))
plot_rarefaction(hh_mob_in, 'treatment', ref_level = 'control-closed', 'sSBR', lwd = 4)
par(mfrow=c(1,2))
plot_rarefaction(hh_mob_in, 'treatment', ref_level = 'control-closed', 'IBR', lwd = 4)
par(mfrow=c(1,2))
plot_rarefaction(hh_mob_in, 'treatment', ref_level = 'control-closed', 'IBR', lwd = 4,
                 scales = 'gamma')


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


comm_div <- calc_comm_div(hh_mob_in$comm, index = c('N', 'S', 'S_n', 'S_PIE', 'S_C'),
                          effort = 20, scale = 'alpha')

hh_dat <- data.frame(comm_div, hh_mob_in$env)

plot(value ~ dist_avg, subset = index == 'S' & treatment != 'upland',
     data = hh_dat)


boxplot(densiometer_avg ~ treatment, subset = index == 'S' & treatment != 'upland',
     data = hh_dat)


par(mfrow=c(4,1))
boxplot(value ~ site, subset = index == 'S' & treatment != 'upland',
     data = hh_dat, main = 'S')
boxplot(value ~ site, subset = index == 'S_n' & treatment != 'upland',
        data = hh_dat, main = 'S_n')
boxplot(value ~ site, subset = index == 'S_PIE' & treatment != 'upland',
        data = hh_dat, main = 'S_PIE')
boxplot(value ~ site, subset = index == 'S_C' & treatment != 'upland',
        data = hh_dat, main = 'S_C')

# Two-Scale Analysis

hh_stats <- get_mob_stats(hh_mob_in, "treatment")
hh_stats

plot(hh_stats, 'treatment')


table(hh_mob_in$env$treatment)

trts <- c('control-closed', 'hack-squirt', 'cut-leave','cut-remove', 'control-open')
hh_stats <- get_mob_stats(subset(hh_mob_in, treatment %in% trts),
                     "treatment", index = c('N', 'S', 'S_n', 'S_C', 'S_PIE'), )
plot(hh_stats, 'treatment', group_order = trts)


library(vegan)

mod <- rda(hh_mob_in$comm ~ hh_mob_in$env$treatment,
           subset = hh_mob_in$env$treatment %in% trts)
anova(mod, by = 'margin')
RsquareAdj(mod)
#$r.squared
#[1] 0.2685041
#
#adj.r.squared
#[1] 0.02467216

plot(mod)




# BACS
mean(hh_mob_in$comm$BACS)



treatments <- unique(hh_mob_in$env$treatment)
sample_list <- hh_mob_in$comm %>%
  group_by(hh_mob_in$env$treatment) %>%
  group_map(~ get_samples(.x, algo = 'boot', n_boot = 500))
names(sample_list) <- treatments

indices <- c('N', 'S', 'S_C', 'S_n', 'S_PIE')
effort <- 20
sample_div <- lapply(sample_list, calc_comm_div_ci, index = indices,
                     effort = effort)
sample_div <- bind_rows(sample_div, .id = 'id')
sample_div

C_min <- min(sample_div$gamma_coverage, na.rm =TRUE)
# rerun analysis for this target coverage - will only change S_C values
sample_div <- lapply(sample_list, calc_comm_div_ci, index = indices,
                     effort = effort, C_target_gamma = C_min)
sample_div <- bind_rows(sample_div, .id = 'id')
sample_div

# alpha and gamma scale comparison
subset(sample_div, scale != 'beta') |> 
  ggplot() + 
  geom_boxplot(aes(x = id, y = me, col=scale)) + 
  geom_errorbar(aes(x = id, y = me, ymin = lo, ymax = hi, col=scale),
                width = 0.2, position=position_dodge(0.8)) + 
  facet_wrap(vars(index), scales = 'free')

# beta diversity comparison
subset(sample_div, scale == 'beta') |> 
  ggplot() + 
  geom_boxplot(aes(x = id, y = me)) + 
  geom_errorbar(aes(x = id, y = me, ymin = lo, ymax = hi),
                width = 0.2, position=position_dodge(0.8)) + 
  facet_wrap(vars(index), scales = 'free')



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




