library(sf)
library(terra)
library(foreach)
library(units)
library(tidyverse)
library(lubridate)

data_folder <- 'data/'

year = "2024"
template <- rast(paste0(data_folder, "avoided_emissions/land_1km_eck4.tif"))

sites_2024 <- readRDS(paste0("data/ci_sites/FY", year, "_Sites_Clean.rds")) 
sites_2024$Data_Year <- 2024

fiscalyear <- "FY2024"


#Load in 2024 as RDS
sites <- sites_2024


sites_cea <- st_transform(sites, '+proj=cea')
sites_cea$area_cea <- st_area(sites_cea)
units(sites_cea$area_cea) <- 'hectares'

sites <- st_transform(sites_cea, crs(template))

table(sites_cea$area_cea < as_units(100, 'hectares'))

sites <-  sites %>% dplyr::select("CI_ID", "Area_Name",
          "CI_Start_Date", "CI_End_Date",
          "Division_Primary_CI", "Intervention_Type", "area_cea", "Data_Year")

sites$CI_ID <- factor(sites$CI_ID)


sites$CI_End_Date_clean <- as.Date(sites$CI_End_Date , format = "%Y/%m/%d")
sites$CI_Start_Date_clean<- as.Date(sites$CI_Start_Date , format = "%Y/%m/%d")

sites %>%
    ggplot() +
    geom_histogram(aes(CI_Start_Date_clean))

median(sites$CI_Start_Date_clean, na.rm = TRUE)

sum(is.na(sites$CI_Start_Date_clean))

# Set all start dates that are missing to 2022 (the median year)
sites$CI_Start_Date_clean[is.na(sites$CI_Start_Date_clean)] <- ymd('2022-07-23')
sites$CI_Start_Year <- year(sites$CI_Start_Date_clean)

sites %>%
    ggplot() +
    geom_histogram(aes(CI_End_Date_clean))

# Set all end dates that are greater than 12/31/2024 to NA, so they are treated 
# as ongoing

sites$CI_End_Date_clean[sites$CI_End_Date_clean > ymd('2024-12-31')] <- NA
sites$CI_End_Year <- year(sites$CI_End_Date_clean)
table(is.na(sites$CI_End_Year))

sites$ID <- 1:nrow(sites)
write_csv(dplyr::select(sites, CI_ID, ID), 'data/avoided_emissions/site_ID_key_added_FY2024.csv')
saveRDS(sites, 'data/avoided_emissions/sites_FY2024.RDS')

