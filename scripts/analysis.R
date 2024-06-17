library(mobr)

dat <- read.csv('./data/bird_data_year4 - IB_data.csv')
# small data cleaning
dat$X25m <- ifelse(is.na(dat$X25m), 0, dat$X25m)
dat$X50m <- ifelse(is.na(dat$X50m), 0, dat$X50m)
dat$date <- as.Date(dat$date, format = "%m/%d/%Y")
head(dat)
# create unique sampling event id


# use mobr::calc_comm_div to compute diversity indices for each sampling event.
# should do this for all species in each site and also for just observations 
# within 25m. Jackson only used 25m observations


