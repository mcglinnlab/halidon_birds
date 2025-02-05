library(mobr)
library(dplyr)
library(vegan)
library(nlme)
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


# make treatments _pre or _post
birds_25p$env$treatment <- ifelse(birds_25p$env$pre_post == "pre", paste(birds_25p$env$treatment, "_pre", sep = ""), 
                              paste(birds_25p$env$treatment, "_post", sep = ""))
print(birds_25p)


# use mobr::calc_comm_div to compute diversity indices for each sampling event.
# should do this for all species in each site and also for just observations 
# within 25m. Jackson only used 25m observations

indices <- c('N', 'S', 'S_n', 'S_PIE', 'S_asymp')
sub_dat <- subset(birds_25p,  treatment != "upland")

sub_dat <- birds_25p$comm

sub_rda <- rda(sub_dat ~  treatment, data =  birds_25p$env, subset = treatment != "upland_pre" & treatment != "upland_post")
plot(sub_rda, display = c('sp', 'cn'))
anova(sub_rda, by = 'terms')

bggn_lm <- lm(sub_dat$BGGN ~ treatment, data =  birds_25p$env, subset = treatment != "upland")
summary(bggn_lm)

boxplot(sub_dat$BHNU ~ treatment, data =  birds_25p$env, subset = treatment != "upland")

stats_raw <- calc_comm_div(sub_dat$comm,
                           index = indices, effort = 5, scales = 'alpha')
tst <- cbind(stats_raw, sub_dat$env)

Nmod <- lm(value ~ treatment + densiometer_avg + dist_avg + dbh_avg + disked ,
            data = tst, subset = index == 'N')
summary(Nmod)
car::Anova(Nmod, type = 3)

Nmodre <- lme(value ~ treatment + disked,
              random = ~1 | year / site, data = tst, na.action = na.omit, 
              subset = index == "S")
summary(Nmodre)
car::Anova(Nmodre, type = 3)


# spatial analysis of 2024 data only at point , scale
# remove uplands from analysis
stats_trt <- get_mob_stats(subset(birds_25p, treatment != "upland_pre" & treatment != "upland_post"), 
                           group_var = 'treatment', 
                           index = c('N', 'S', 'S_n', 'S_PIE', 'S_asymp'),
                           ci_n_boot = 100)

# no apparent treatment effects
plot(stats_trt, group_var = 'treatment')

# plot _pre & _post treatments next to each other for comparison
# N index ----
N_data <- subset(stats_trt$comm_div, scale == 'alpha' & index == 'N')
# Create a bar plot
barplot_N <- barplot(N_data$value, 
        names.arg = N_data$treatment, 
        beside = TRUE, 
        col = "lightblue", 
        ylab = "N", 
        las = 2)  
arrows(barplot_N, N_data$lo_value, barplot_N, N_data$hi_value, 
       angle = 90, code = 3, length = 0.1, col = "black")
# S index ----
S_data <- subset(stats_trt$comm_div, scale == 'alpha' & index == 'S')
# Create a bar plot
barplot_S <- barplot(S_data$value, 
                     names.arg = S_data$treatment, 
                     beside = TRUE, 
                     col = "lightblue",
                     ylab = "S", 
                     las = 2)  
arrows(barplot_S, S_data$lo_value, barplot_S, S_data$hi_value, 
       angle = 90, code = 3, length = 0.1, col = "black")
# S_n index ----
S_n_data <- subset(stats_trt$comm_div, scale == 'alpha' & index == 'S_n')
# Create a bar plot
barplot_S_n <- barplot(S_n_data$value, 
                     names.arg = S_n_data$treatment, 
                     beside = TRUE, 
                     col = "lightblue", 
                     ylab = "S_n", 
                     las = 2)  
arrows(barplot_S_n, S_n_data$lo_value, barplot_S_n, S_n_data$hi_value, 
       angle = 90, code = 3, length = 0.1, col = "black")
# S_PIE index ----
S_PIE_data <- subset(stats_trt$comm_div, scale == 'alpha' & index == 'S_PIE')
# Create a bar plot
barplot_S_PIE <- barplot(S_PIE_data$value, 
                     names.arg = S_PIE_data$treatment, 
                     beside = TRUE, 
                     col = "lightblue", 
                     ylab = "S_PIE", 
                     las = 2)  
arrows(barplot_S_PIE, S_PIE_data$lo_value, barplot_S_PIE, S_PIE_data$hi_value, 
       angle = 90, code = 3, length = 0.1, col = "black")
# S_asymp index ----
S_asymp_data <- subset(stats_trt$comm_div, scale == 'alpha' & index == 'S_asymp')
# Create a bar plot
barplot_S_asymp <- barplot(S_asymp_data$value, 
                         names.arg = S_asymp_data$treatment, 
                         beside = TRUE, 
                         col = "lightblue",
                         ylab = "S_asymp", 
                         las = 2)  
arrows(barplot_S_asymp, S_asymp_data$lo_value, barplot_S_asymp, S_asymp_data$hi_value, 
       angle = 90, code = 3, length = 0.1, col = "black")
# ----

# todo: 
# recode treatments to be more informative
# look at temporal change

plot(1:10, 1:10, col = "#1462AE", pch =19)

# temporal analysis of pre / post 
stats_pp <- get_mob_stats(birds_25p, 
                           group_var = 'pre_post', 
                           index = c('N', 'S', 'S_n', 'S_PIE', 'S_C'),
                           ci_n_boot = 100)

svg("./figs/temp_S_plot.svg", width = 7*1.5, height = 5)
plot(stats_pp, group_var = 'pre_post', index = 'S')
dev.off()

svg("./figs/temp_N_plot.svg", width = (7.5*1.5)*0.66, height = 5)
plot(stats_pp, group_var = 'pre_post', index = 'N')
dev.off()

# rda
bird_rda <- rda(sub_dat ~  treatment, data =  birds_25p$env, subset = treatment != "upland_pre" & treatment != "upland_post")
plot(bird_rda, display = c('sp', 'cn'))
anova(bird_rda, by = 'terms')
RsquareAdj(bird_rda)

plot(bird_rda, display = 'species', type = 'n')
orditorp(bird_rda, display = 'species', )
points(bird_rda, display = 'bp', col = 'red') 
text(bird_rda, display = 'cn', col = 'red') 

# cca 
bird_cca<- cca(comm_25p ~ hh_attp$pre_post + hh_attp$treatment) 
anova(bird_cca) 
anova(bird_cca, by='terms')
RsquareAdj(bird_cca)

plot(bird_cca, display = 'species', type = 'n')
orditorp(bird_cca, display = 'species', )
points(bird_cca, display = 'bp', col = 'red') 
text(bird_cca, display = 'cn', col = 'red') 


boxplot(comm_25p$BACS ~ hh_attp$year, subset = comm_25p$BACS > 0)

# multiple regression modeling of diversity indices
subset(stats_trt$comm_div, index = "N")




