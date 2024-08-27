library(mobr)
library(dplyr)
#library(remotes)
#install_github("mobiodiv/mobr", ref = "dev")

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

# calculate div stats
div_stats <-  calc_comm_div(comm, index = c("N", "S", "S_n", "S_C", "S_PIE"), effort = 5)

# group by different treatments
grouped_data <- dat %>%
  group_by(treatment)


# use mobr::calc_comm_div to compute diversity indices for each sampling event.
# should do this for all species in each site and also for just observations 
# within 25m. Jackson only used 25m observations


