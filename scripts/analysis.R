library(mobr)
library(dplyr)
library(vegan)
#library(remotes)
#install_github("mobiodiv/mobr", ref = "dev")

# start at point count scale
comm_25p <- read.csv('./data/comm_25p.csv', row.names = 1)
comm_25p[1:5, 1:5]

hh_attp <- read.csv('./data/hh_attp.csv')

# group hack and squirt with control closed treatment
hh_attp <- hh_attp %>%
  mutate(treatment = ifelse(treatment == "hack-squirt", "control-closed", treatment)) %>%
  group_by(treatment)


birds_25p <- make_mob_in(comm_25p, hh_attp,
                       coord_names = c('utm_easting', 'utm_northing'))

# use mobr::calc_comm_div to compute diversity indices for each sampling event.
# should do this for all species in each site and also for just observations 
# within 25m. Jackson only used 25m observations


# spatial analysis of 2024 data only at point count scale
# remove uplands from analysis
stats_trt <- get_mob_stats(subset(birds_25p, year == 2024 & treatment != "upland"), 
                           group_var = 'treatment', 
                           index = c('N', 'S', 'S_n', 'S_PIE', 'S_C'),
                           ci_n_boot = 100)
# no apparent treatment effects
plot(stats_trt, group_var = 'treatment')


# todo: 
# recode treatments to be more informative
# look at temporal change


# temporal analysis of pre / post 
stats_pp <- get_mob_stats(birds_25p, 
                           group_var = 'pre_post', 
                           index = c('N', 'S', 'S_n', 'S_PIE', 'S_C'),
                           ci_n_boot = 100)

plot(stats_pp, 'pre_post')


bird_rda <- rda(comm_25p ~ hh_attp$pre_post + hh_attp$treatment)
anova(bird_rda)
anova(bird_rda, by='terms')
RsquareAdj(bird_rda)

plot(bird_rda, display ='species')
points(bird_rda, display ='bp', col='red')
text(bird_rda, display ='bp', col='red')


boxplot(comm_25p$BACS ~ hh_attp$year, subset = comm_25p$BACS > 0)



