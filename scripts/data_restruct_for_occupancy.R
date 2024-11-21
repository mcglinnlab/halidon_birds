library(tidyr)

comm_25p <- read.csv('./data/comm_25p.csv', row.names = 1)
hh_attp <- read.csv('./data/hh_attp.csv')

# create a variable called visit
visit <- NULL
uni_sites <- unique(hh_attp$site)
uni_yrs <- unique(hh_attp$year)
for(i in seq_along(uni_sites)) {
    for(j in seq_along(uni_yrs)) {
       true <- hh_attp$site == uni_sites[i] & 
               hh_attp$year == uni_yrs[j]
       visit <- c(visit, order(hh_attp$date[true]))
  }
}

# within each year pull out the three surveys into 3 seperate columns
#chose a species
NOCA <- data.frame(y = comm_25p$NOCA, visit, hh_attp)

tmp <- pivot_wider(subset(NOCA, select = -date),
                   names_from = visit, values_from = y)

as.data.frame(tmp)
