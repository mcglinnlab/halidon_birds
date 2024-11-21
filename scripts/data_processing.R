library(googlesheets4)

# read in raw data
hh_abu <- read_sheet("https://docs.google.com/spreadsheets/d/1LsEtgikV3PtJJPZsW5IFkuCzZSZTZVg3Y2vQdgvS3V8/edit?usp=sharing", 
                                sheet = "all_data")
# Site Attribute Table
hh_att <- read_sheet("https://docs.google.com/spreadsheets/d/1LsEtgikV3PtJJPZsW5IFkuCzZSZTZVg3Y2vQdgvS3V8/edit?usp=sharing", 
                     sheet = "attributes")

head(hh_abu)
names(hh_abu)
# change NAs in count column to 0
hh_abu$`25m` <- ifelse(is.na(hh_abu$`25m`), 0, hh_abu$`25m`)
hh_abu$`50m` <- ifelse(is.na(hh_abu$`50m`), 0, hh_abu$`50m`)
# fix data formating
hh_abu$date <- as.Date(hh_abu$date, format = "%m/%d/%Y")

# species to drop
sp_to_drop <- c("AMCR", "BBWD", "WOST", "NOBO", "CAGO", 
                "CONI", "FICR", "GREG", "MIKI", "RSHA",
                "RTHA", "LAGU")
hh_abu <- subset(hh_abu, !(species %in% sp_to_drop))

# fix site ids
hh_abu$site <- paste("HH", substring(hh_abu$site, 3, 4), 
                     ifelse("UP" == substring(hh_abu$site, 1, 2), "U", "W"),
                     sep="-")

hh_abu$site_date <- with(hh_abu, paste(site, date, sep='_'))

table(hh_abu$site, hh_abu$site_type)
table(hh_abu$site, hh_abu$year)

# check hh_abu for data issues
which(hh_abu$total != (hh_abu$`25m` + hh_abu$`50m`))
sum(hh_abu$total == (hh_abu$`25m` + hh_abu$`50m`))
nrow(hh_abu)
sum(hh_abu$`25m`)
sum(hh_abu$`50m`)


# Remove wetland 43 because it was mulched in 2024
hh_abu <- subset(hh_abu, !(site %in% c('HH-43-W', 'HH-43-U')))
# remove wetland 5 because it was a control disturbed (only one replicate of this type)
hh_abu <- subset(hh_abu, !(site %in% c('HH-05-W', 'HH-05-U')))


# Create Community Matrices -----
# first create these at the site level
comm_25 <- with(hh_abu, tapply(`25m`, list(site, species), sum))
comm_25 <- ifelse(is.na(comm_25), 0, comm_25)

# going to combine 25 and 50 m radius by using total
comm_tot <- with(hh_abu, tapply(total, list(wetland_id, species), sum))
comm_tot <- ifelse(is.na(comm_tot), 0, comm_tot)

# now create them at the point count scale
comm_25p <- with(hh_abu, tapply(`25m`, list(site_date, species), sum))
comm_totp <- with(hh_abu, tapply(total, list(site_date, species), sum))

comm_25p <- ifelse(is.na(comm_25p), 0, comm_25p)
comm_totp <- ifelse(is.na(comm_totp), 0, comm_totp)


# create attribute table on the point count scale ----
head(hh_att)
sites <- substring(row.names(comm_25p), 1, 7)
dates <- as.Date(substring(row.names(comm_25p), 9 , 18))
years <- substring(row.names(comm_25p), 9, 12)
pre_post <- ifelse(years < 2024, 'pre', 'post')
table(pre_post)

tmp <- data.frame(site = sites, date = dates, year = years, 
                  pre_post)

hh_attp <- dplyr::left_join(tmp, hh_att, by = c('site', 'pre_post'))
all.equal(hh_attp$site, sites)
head(hh_attp)
# create canopy cover variable and 
#bring in lidar cover under a single variable

lidar <- read.csv('https://raw.githubusercontent.com/mcglinnlab/wetland_birds/refs/heads/main/data/lidar_tree_metrics.csv')
head(lidar)
names(lidar)[1] <- 'site_old'

tmp <- dplyr::left_join(hh_attp, lidar, by = 'site_old')

hh_attp$canopy_cover <- with(tmp, 
                             ifelse(is.na(densiometer_avg),
                             canopy_cover_sim,
                             100 - hh_attp$densiometer_avg))

# export files -----
write.csv(comm_25, file='./data/comm_25.csv', row.names = TRUE)
write.csv(comm_tot, file='./data/comm_tot.csv', row.names = TRUE)
write.csv(comm_25p, file='./data/comm_25p.csv', row.names = TRUE)
write.csv(comm_totp, file='./data/comm_totp.csv', row.names = TRUE)
write.csv(hh_att, file = './data/hh_att.csv', row.names = FALSE)
write.csv(hh_attp, file = './data/hh_attp.csv', row.names = FALSE)

