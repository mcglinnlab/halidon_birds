library(mobr)
library(dplyr)
library(vegan)
library(nlme)
library(ggplot2)
#library(remotes)
#install_github("mobiodiv/mobr", ref = "dev")

# start at point count scale
comm_25p <- read.csv('./data/comm_25p.csv', row.names = 1)
comm_25p[1:5, 1:5]

# drop the species that never occur
comm_25p <- comm_25p[, colSums(comm_25p) > 0]


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
sub_dat <- subset(birds_25p, 
                  !(treatment %in% c("upland_pre", "upland_post")))

#sub_dat <- birds_25p$comm

sub_rda <- rda(sub_dat$comm ~  treatment, data =  sub_dat$env)
RsquareAdj(sub_rda)
anova(sub_rda, by = 'terms')

par(mfrow=c(1,1))
plot(sub_rda, display = c('cn','sp'), type = 'n')
text(sub_rda, display = 'cn', col = 'red') 
orditorp(sub_rda, display = 'sp', priority = colSums(sub_dat$comm), air = 1.5)
#text(sub_rda, display = 'cn', col = cols)

# indirect ordination
#
sub_pca <- rda(sqrt(sub_dat$comm))
#sub_envfit <- envfit(sub_pca, sub_dat$env[ , c("treatment")], na.rm = TRUE)
plot(sub_pca, display = 'sp')
# indirect ordination does not differ strongly from the direct ordination using the
# rda method

# rarefaction analysis ----------------
# check if any samples have zero birds

table(rowSums(sub_dat$comm))

# 7 sites have no individuals - 
# this is a problem for the method
# drop these sites real quick

sub_dat2 <- subset(sub_dat, rowSums(sub_dat$comm) > 0)

# getting cryptic bug here
#tmp <- get_delta_stats(sub_dat2, 'treatment',
#                type = 'discrete', log_scale = FALSE, n_perm = 3,
#                overall_p = TRUE)

# go back to sites with no observations
S_op_pre <- specaccum(sub_dat$comm[sub_dat$env$treatment == 'control-open_pre', ])
S_op_post <- specaccum(sub_dat$comm[sub_dat$env$treatment == 'control-open_post', ])
S_cl_pre <- specaccum(sub_dat$comm[sub_dat$env$treatment == 'control-closed_pre', ])
S_cl_post <- specaccum(sub_dat$comm[sub_dat$env$treatment == 'control-closed_post', ])
S_CL_pre <- specaccum(sub_dat$comm[sub_dat$env$treatment == 'cut-leave_pre', ])
S_CL_post <- specaccum(sub_dat$comm[sub_dat$env$treatment == 'cut-leave_post', ])
S_CR_pre <- specaccum(sub_dat$comm[sub_dat$env$treatment == 'cut-remove_pre', ])
S_CR_post <- specaccum(sub_dat$comm[sub_dat$env$treatment == 'cut-remove_post', ])

# remove error bars
plot(S_op_pre, ylim = c(0, 30), xlim = c(0, 40), ylab = "Species Richness (S)",
     xlab = "Samples", ci = 0, lwd = 2, col= 'grey')
lines(S_op_post, col='black', ci = 0, lwd = 2)
lines(S_cl_pre, col='pink', ci = 0, lwd = 2)
lines(S_cl_post, col='red', ci = 0, lwd = 2)
lines(S_CL_pre, col='lightblue', ci = 0, lwd = 2)
lines(S_CL_post, col='blue', ci = 0, lwd = 2)
lines(S_CR_pre, col='lightgreen', ci = 0, lwd = 2)
lines(S_CR_post, col='darkgreen', ci = 0, lwd = 2)
legend("bottomright", legend = c("Control-Open Pre", "Control-Open Post", "Control-Closed Pre",
                              "Control-Closed Post", "Cut-Leave Pre", "Cut-Leave Post", 
                              "Cut-Remove Pre", "Cut-Remove Post"), 
       col = c("grey", "black", "pink", "red", "lightblue", "blue", "lightgreen", "darkgreen"), 
       lwd = 2, bty = "n")

# new rarefaction with less colors
plot(S_op_pre, ylim = c(0, 30), xlim = c(0, 40), ylab = "Species Richness (S)",
     xlab = "Samples", ci = 0, lwd = 2, col= 'blue', lty= 'dotted')
lines(S_op_post, col='blue', ci = 0, lwd = 2)
lines(S_cl_pre, col='orange', ci = 0, lwd = 2, lty= 'dotted')
lines(S_cl_post, col='orange', ci = 0, lwd = 2)
lines(S_CL_pre, col='green', ci = 0, lwd = 2, lty= 'dotted')
lines(S_CL_post, col='green', ci = 0, lwd = 2)
lines(S_CR_pre, col='purple', ci = 0, lwd = 2, lty= 'dotted')
lines(S_CR_post, col='purple', ci = 0, lwd = 2)
legend("bottomright", legend = c("Control-Open", "Control-Closed", 
                                 "Cut-Leave", "Cut-Remove", "Pre", "Post"), 
       col = c("blue", "orange", "green", "purple", "black", "black"), 
       lwd = 2, lty = c(1, 1, 1, 1, 2, 1), bty = "n")

#bggn_lm <- lm(sub_dat$BGGN ~ treatment, data =  birds_25p$env, subset = treatment != "upland")
#summary(bggn_lm)

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
# split out pre-post variable
stats_trt$comm_div$pre_post <- sapply(strsplit(stats_trt$comm_div$treatment, '_'),
                                      function(x) x[2])
stats_trt$comm_div$pre_post <- factor(stats_trt$comm_div$pre_post, 
                                      levels = c('pre', 'post'))
stats_trt$comm_div$trt <- sapply(strsplit(stats_trt$comm_div$treatment, '_'),
                                 function(x) x[1])

stats_trt$comm_div$treatment <- factor(stats_trt$comm_div$treatment, 
                  levels = c('control-closed_pre', 
                             'control-closed_post', 
                             'control-open_pre',
                             'control-open_post',
                             'cut-leave_pre',
                             'cut-leave_post',
                             'cut-remove_pre',
                             'cut-remove_post'))

indices <- unique(stats_trt$comm_div$index)
labs <- c('Abundance (N)' ,'Species Richness (S)', '', 'Rarefied Richness (S_n)', 
          '', 'Asymptotic Richness (S_asymp)', 'Evenness (S_PIE)')
p <- vector("list", length(indices))
names(p) <- indices

for(i in seq_along(indices)) {

  p[[i]] <- 
    subset(stats_trt$comm_div,
           subset = scale == 'alpha' & index == indices[i]) %>%
    ggplot(aes(x = treatment, y = value)) + 
            geom_bar(stat = "identity", aes(fill = pre_post)) +
            geom_errorbar(aes(x=treatment, ymin=lo_value, ymax=hi_value), width=0.15, 
                          colour="black", alpha=0.7, size=0.5) +
            ylab(labs[i]) + theme_bw() + xlab("Treatment") +
            scale_fill_manual(name = "Pre/Post" , 
                              values = alpha(c("lightblue", "darkblue"))) +
            theme(axis.title = element_text(face = "bold"))
}

p$N
p$S
p$S_n
p$S_PIE
p$S_asymp

# model year and observer differences (if any)
stats_obs <- get_mob_stats(subset(birds_25p, treatment != "upland_pre" & treatment != "upland_post"), 
                           group_var = 'pre_post', 
                           index = c('N', 'S', 'S_n', 'S_PIE', 'S_asymp'),
                           ci_n_boot = 100)

plot(stats_obs, group_var = "pre_post")

indices_obs <- unique(stats_obs$comm_div$index)
labs_obs <- c('Abundance (N)' ,'Species Richness (S)', '', 'Rarefied Richness (S_n)', 
          '', 'Asymptotic Richness (S_asymp)', 'Evenness (S_PIE)')
p_obs <- vector("list", length(indices_obs))
names(p_obs) <- indices

for(i in seq_along(indices_obs)) {
  
  p_obs[[i]] <- 
    subset(stats_obs$comm_div,
           subset = scale == 'alpha' & index == indices_obs[i]) %>%
    ggplot(aes(x = factor(pre_post, levels = rev(levels(factor(pre_post)))), 
               y = value, fill = pre_post)) +  # Ensuring colors apply
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin=lo_value, ymax=hi_value), width=0.15, 
                  colour="black", alpha=0.7, linewidth=0.5) +
    ylab(labs[i]) + theme_bw() + xlab("Treatment") +
    scale_x_discrete(labels = c("pre" = "Barratt Heitmann", "post" = "Barnes")) +
    scale_fill_manual(name = "Pre/Post", 
                      values = c("brown", "tan")) +  # Remove alpha() for proper color display
    theme(axis.title = element_text(face = "bold"))
}


p_obs$N
p_obs$S
p_obs$S_n

# observer S_PIE plot
obs_S_PIE <- 
  subset(stats_obs$comm_div,
         subset = scale == 'alpha' & index == indices_obs[7]) %>%
  ggplot(aes(x = factor(pre_post, levels = rev(levels(factor(pre_post)))), 
             y = value, fill = pre_post)) +  # Ensuring colors apply
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin=lo_value, ymax=hi_value), width=0.15, 
                colour="black", alpha=0.7, linewidth=0.5) +
  ylab("Evenness (S_PIE)") + theme_bw() + xlab("Treatment") +
  scale_x_discrete(labels = c("pre" = "Barratt Heitmann", "post" = "Barnes")) +
  scale_fill_manual(name = "Pre/Post", 
                    values = c("brown", "tan")) +  # Remove alpha() for proper color display
  theme(axis.title = element_text(face = "bold"))
obs_S_PIE

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
bird_rda <- rda(sub_dat ~  treatment, data =  birds_25p$env, 
                subset = treatment != "upland_pre" & treatment != "upland_post")
plot(bird_rda, display = c('sp', 'cn'))
anova(bird_rda, by = 'terms')
RsquareAdj(bird_rda)

svg("./figs/RDA_plot.svg", width = 7*1.5, height = 5)
plot(bird_rda, display = 'species', type = 'n')
orditorp(bird_rda, display = 'species', )
points(bird_rda, display = 'bp', col = 'red') 
text(bird_rda, display = 'cn', col = 'red') 
dev.off()

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

ggplot(birds_25p, aes(x = treatment, y = comm_25p$BACS))

# multiple regression modeling of diversity indices
subset(stats_trt$comm_div, index = "N")

# vegetation change between pre/post over time figure
hh_attp$pre_post <- factor(hh_attp$pre_post, levels = c("pre", "post"))

ggplot(hh_attp, aes(x = treatment, y = canopy_cover, fill = pre_post)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge") +
  theme_minimal() +
  labs(x = "Treatment", y = "Canopy Cover", title = "Average Canopy Cover Before and After Treatment") +
  scale_fill_manual(name = "Pre/Post", values = c("pre" = "lightblue", "post" = "darkblue")) +
  theme(axis.title = element_text(face = "bold"))

ggplot(hh_attp, aes(x = treatment, y = dist_avg, fill = pre_post)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge") +
  theme_minimal() +
  labs(x = "Treatment", y = "Tree Density^-1", title = "Average Tree Density^-1 Before and After Treatment") +
  scale_fill_manual(name = "Pre/Post", values = c("pre" = "lightblue", "post" = "darkblue")) +
  theme(axis.title = element_text(face = "bold"))

ggplot(hh_attp, aes(x = treatment, y = dbh_avg, fill = pre_post)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge") +
  theme_minimal() +
  labs(x = "Treatment", y = "Tree Basal Area", title = "Average Basal Area Before and After Treatment") +
  scale_fill_manual(name = "Pre/Post", values = c("pre" = "lightblue", "post" = "darkblue")) +
  theme(axis.title = element_text(face = "bold"))

## BACS means as a function of treatment
birds_25p$env$treatment2 <- with(birds_25p$env, paste(treatment, pre_post, sep='_'))

avgs <- tapply(birds_25p$comm$BACS, list(birds_25p$env$treatment2), mean)
se <- tapply(birds_25p$comm$BACS, list(birds_25p$env$treatment2), function(x) sd(x)/sqrt(length(x)))

avgs <- avgs[c(2,1, 4,3, 6,5, 8,7, 10,9)]
se <- se[c(2,1, 4,3, 6,5, 8,7, 10,9)]

dat <- data.frame(grp = names(avgs), abu = avgs, se,
                  pre_post = sapply(strsplit(names(avgs), '_'), function(x)x[2]))
dat$grp <- factor(dat$grp, levels = names(avgs)) # to reorder the group variable (again :()
                 

ggplot(dat, aes(x = grp, y = abu, fill = pre_post)) +
 geom_bar(stat = "identity", alpha = 0.7, color = "black") + 
  geom_errorbar(aes(ymin = abu - se, ymax = abu + se), width = 0.4, colour = "black",
    alpha = 0.9, size = 1) + scale_fill_manual(name = "Pre/Post", 
                                                 values = c("pre" = "lightblue", 
                                                            "post" = "darkblue")) +
  theme_bw() + labs( x = "Treatment", y = "BACS Abundance") 

## BHNU means as a function of treatment

avgs <- tapply(birds_25p$comm$BHNU, list(birds_25p$env$treatment2), mean)
se <- tapply(birds_25p$comm$BHNU, list(birds_25p$env$treatment2), function(x) sd(x)/sqrt(length(x)))

avgs <- avgs[c(2,1, 4,3, 6,5, 8,7, 10,9)]
se <- se[c(2,1, 4,3, 6,5, 8,7, 10,9)]

dat <- data.frame(grp = names(avgs), abu = avgs, se,
                  pre_post = sapply(strsplit(names(avgs), '_'), function(x)x[2]))
dat$grp <- factor(dat$grp, levels = names(avgs)) # to reorder the group variable (again :()


ggplot(dat, aes(x = grp, y = abu, fill = pre_post)) +
  geom_bar(stat = "identity", alpha = 0.7, color = "black") + 
  geom_errorbar(aes(ymin = abu - se, ymax = abu + se), width = 0.4, colour = "black",
                alpha = 0.9, size = 1) + scale_fill_manual(name = "Pre/Post", 
                                                           values = c("pre" = "lightblue", 
                                                                      "post" = "darkblue")) +
  theme_bw() + labs( x = "Treatment", y = "BHNU Abundance") 

## INBU means as a function of treatment

habitat <- ifelse(grepl('W', birds_25p$env$site), 'wetland', 'upland')

avgs <- tapply(birds_25p$comm$INBU, list(habitat), mean)
se <- tapply(birds_25p$comm$BHNU, list(habitat), function(x) sd(x)/sqrt(length(x)))

#avgs <- avgs[c(2,1, 4,3, 6,5, 8,7, 10,9)]
#se <- se[c(2,1, 4,3, 6,5, 8,7, 10,9)]

dat <- data.frame(grp = names(avgs), abu = avgs, se,
                  pre_post = sapply(strsplit(names(avgs), '_'), function(x)x[2]))
dat$grp <- factor(dat$grp, levels = names(avgs)) # to reorder the group variable (again :()


ggplot(dat, aes(x = grp, y = abu, fill = grp)) +
  geom_bar(stat = "identity", alpha = 0.7, color = "black") + 
  geom_errorbar(aes(ymin = abu - se, ymax = abu + se), width = 0.4, colour = "black",
                alpha = 0.9, size = 1) + scale_fill_manual(name = "Pre/Post", 
                                                           values = c("wetland" = "cornflowerblue", 
                                                                      "upland" = "seagreen")) +
  theme_bw() + labs( x = "Treatment", y = "INBU Abundance") 

