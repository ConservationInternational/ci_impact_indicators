library(sf)
library(terra)
library(foreach)
library(units)
library(tidyverse)
library(lubridate)

data_folder <- 'data'

template <- rast(paste0(data_folder, "avoided_emissions/land_1km_eck4.tif"))

sites_2022 <- st_read("data/ci_sites/FY2022_Sites.shp") 
sites_2022$data_year <- 2022

sites <- sites_2022

sites_cea <- st_transform(sites, '+proj=cea')
sites_cea$area_cea <- st_area(sites_cea)
units(sites_cea$area_cea) <- 'hectares'

sites <- st_transform(sites_cea, crs(template))

table(sites_cea$area_cea < as_units(100, 'hectares'))

sites %>%
    select(CI_ID=ci_id,
           Data_Year=data_year,
           Area=area_name,
           CI_Start_Date=ci_start_d,
           CI_End_Date=ci_end_dat,
           CI_Division=ci_divisio,
           Restoration=restoratio,
           Area_ha=area_cea) -> sites
sites$CI_ID <- factor(sites$CI_ID)

sites$CI_Start_Date_clean <- as.character(sites$CI_Start_Date)
sites$CI_Start_Date_clean <- str_replace(sites$CI_Start_Date_clean, '^([0-9]{4})$', '1/1/\\1')
sites$CI_Start_Date_clean <- lubridate::ymd(sites$CI_Start_Date_clean)
sites %>%
    ggplot() +
    geom_histogram(aes(CI_Start_Date_clean))
# Set all start dates that are missing to 2016 (the median year)
sites$CI_Start_Date_clean[is.na(sites$CI_Start_Date_clean)] <- ymd('2016-01-01')
sites$CI_Start_Year <- year(sites$CI_Start_Date_clean)

sites$CI_End_Date_clean <- as.character(sites$CI_End_Date)
sites$CI_End_Date_clean <- str_replace(sites$CI_End_Date_clean, '^([0-9]{4})$', '1/1/\\1')
sites$CI_End_Date_clean <- ymd(sites$CI_End_Date_clean)
sites %>%
    ggplot() +
    geom_histogram(aes(CI_End_Date_clean))
# Set all end dates that are greater than 12/31/2019 to NA, so they are treated 
# as ongoing
sites$CI_End_Date_clean[sites$CI_End_Date_clean > ymd('2022-12-31')] <- NA
sites$CI_End_Year <- year(sites$CI_End_Date_clean)
table(is.na(sites$CI_End_Year))

sites$ID <- 1:nrow(sites)
write_csv(select(sites, CI_ID, ID), 'data/avoided_emissions/site_ID_key.csv')
saveRDS(sites, 'data/avoided_emissions/sites.RDS')

# Check for overlaps
# intersections <- foreach (year in c(2019, 2019)) %do% {
#     these_sites <- sites_cea[(!sites$area_cea_lt_100ha) & (sites$data_year == 2018), ]
#     return st_intersects(these_sites, these_sites)
# }
