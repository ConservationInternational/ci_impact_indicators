# Data preparation for institutional indicators analysis
# Contact: Anna Ballasiotes
# Last updated: 04/7/2025
# Purpose: This script preps some of the analytical datasets
# as well as the Intervention Sites

library(tidyverse)
library(raster)
library(terra)
library(sf)
library(janitor)
library("aws.s3")
library(dplyr)

#########################
#### Analytical Data ####
#########################

## Outputs: 
# IC carbon
# Woody Carbon, Soil Carbon (later added for total carbon)
# Carbon Sequestration Potential
# Population
# 


###############
#### Sites ####
###############
## Inputs:
# CI Sites Feature Class from geodatabase

## Outputs: 
# FYxxxx_Overlaps.rds
# FYxxxx_Overlaps_Clean.rds
# FYxxxx_Sites_Clean.rds



# Update for Fiscal Year

year <- "2024"

#####

# I would recommend getting the newest version of IC. Contact Monica Noon.

my_crs <- crs(rast("data/irrecoverable_carbon/tstor_ic_ecosystem_stack.tif"))

#####
##### Carbon Stored
#####

## Woody biomass
# Global forest watch Aboveground Live Woody Biomass stock in megagrams/pixel (tonnes/pixel)

## THIS SHOULD BE DONE IN GEE. IT IS A WASTE OF PROCESSING POWER AND HARD DRIVE
## SPACE TO DO THIS IN R.
## CODE IS HERE: https://code.earthengine.google.com/7473ada262f98876b41782d6b15d26b1 

# mosaic GEE output together

rast_list <- c()
rast_list <- list.files(path="data/carbon_stored", 
                        pattern = "^biomass_carbon_gfw_bc_300m.*\\.tif$", full.names = TRUE)

# # Here is how you return just the files names
# rastlist.names <- c( unlist( lapply(str_sub(rast_list,65, -5), FUN=function(x) { x[1] })))   
# print(rastlist.names)
# 
# # Making rasters for each file, if you need to.
# for (i in 1:length(rast_list)){
#   assign(rastlist.names[i], rast(rast_list[i]))
# }

## Use terra to mosaic the carbon into a single data layer
rsrc <- terra::sprc(rast_list)
carbon_bio_mosaic <- mosaic(rsrc)

# I can't remember if I figured out why the carbon mosaic was multiplied by 0.1
# to calculate this sum. Probably a units thing. 
carbon_bio_mosaic <- carbon_bio_mosaic * 0.1
global(carbon_bio_mosaic, "sum", na.rm = TRUE)

plot(carbon_bio_mosaic)
print(carbon_bio_mosaic)

writeRaster(carbon_bio_mosaic, paste0("data/carbon_stored/biomass_prepped_", year, ".tif"))


## If reprojection needed:

# carbon_bio_mosaic_rp <- carbon_bio_mosaic %>% 
#   project(biomass_prepped_2022,
#           align = TRUE,
#           method = "bilinear",
#           threads = TRUE) %>% 
#   extend(biomass_prepped_2022) %>% 
#   crop(biomass_prepped_2022)
# global(carbon_bio_mosaic_rp, "sum", na.rm = TRUE)
# 
# 
# writeRaster(carbon_bio_mosaic_rp, paste0("data/carbon_stored/biomass_prepped_", year, ".tif"))




## Soil biomass

# # soil depth 1-3 for 30cm
# # ocstha = organic soil carbon tonnes/ha
# 
# sd1 <- rast("data/carbon_stored/OCSTHA_M_sd1_1km_soc.tif")
# sd2 <- rast("data/carbon_stored/OCSTHA_M_sd2_1km_soc.tif")
# sd3 <- rast("data/carbon_stored/OCSTHA_M_sd3_1km_soc.tif")
# # crs(sd3) == my_crs # TRUE
# 
# soil <- (sd1 + sd2 + sd3) * 100
# 
# writeRaster(soil, "data/carbon_stored/OCSTHA_30cm_1km.tif")


#####

##### Irrecoverable Carbon
#####

# # If you pull in new IC data, you should run this. 

# # These are the Dinerstein ecoregions
#eco <- rast("data/ecosystems/ecosystems.tif") 
# 
# # total tonnes of carbon per hectare
# total_ic <- rast("data/irrecoverable_carbon/Irrecoverable_C_Total_2018.tif") %>% 
#   resample(eco, method = "bilinear")
# # checked multiple methods of resampling + combining with area and this was the most accurate
# # at achieving 139B tonnes
# 
# # convert to tonnes per pixel
# area <- raster::area(raster(total_ic)) * 100 #km2 to ha
# 
# tonnes_ic <- raster(total_ic) * area
# 
# global(rast(tonnes_ic), "sum", na.rm = TRUE)
# 
# writeRaster(tonnes_ic, "tonnes_ic_bilmethod.tif")
# 
# total_ic <- rast("tonnes_ic_bilmethod.tif") 
# 
# # ic total with vals below 25 tonnes/ha as NA
# high_ic <- rast("data/irrecoverable_carbon/ic_above_25.tif") %>% 
#   resample(eco, method = "bilinear")
# 
# # ic_total w/ vals below 25 as NA
# 
# # make raster of 'any' ic (above 0.01)
# any_ic <- total_ic
# any_ic[any_ic > 0.01] <- 1
# any_ic[any_ic <= 0.01] <- NA
# 
# 
# # change ic_high to binary
# high_ic[high_ic > 0] <- 1
# 
# # need 'high' and 'any' to reflect has/pixel
# # change ic_high and ic_any to represent area
# 
# high_area <- cellSize(high_ic, unit = "ha")
# 
# any_area <- cellSize(any_ic, unit = "ha")
# 
# 
# # create stack
# ic_stack <- c(total_ic, high_area, any_area)  
# 
# names(ic_stack) <- c("tstor_ic", "ha_high_ic", "ha_ic")
# 
# writeRaster(ic_stack, "data/irrecoverable_carbon/ic_stack_prepped.tif",
#             overwrite = TRUE)
# 
# ic_stack <- rast("data/irrecoverable_carbon/ic_stack_prepped.tif")
# 
# # check total to confirm it matches with our 140B tonnes total in Noon et al
# global(ic_stack$tstor_ic, "sum", na.rm = TRUE)
# 
# # divvy up irrecoverable carbon by ecosystem
# 
# ### ecosystems to mask irr carbon to
# 
# eco <- rast("data/ecosystems/ecosystems.tif") 
# eco_codes <- data.frame(cats(eco$class)) %>%
#   rename('ecosystem_name' = class) %>% 
#   drop_na()
# 
# eco_seg <- eco %>%
#   terra::segregate()
# 
# eco_seg[eco_seg == 0] <- NA
# 
# writeRaster(eco_seg, "eco_seg_nas.tif", overwrite = TRUE)
# 
# eco_seg <- rast("eco_seg_nas.tif")
# 
# # for each ic layer
# for (i in 1:nlyr(ic_stack)){
#   
#   ic_lyr <- ic_stack[[i]]
#   og_name <- names(ic_lyr)
#   ic_eco_stack <- ic_lyr
#   
#   # mask irr carbon to each ecosystem layer then add to stack
#   for(e in 1:nlyr(eco_seg)){
#     eco_mask <- eco_seg[[e]]
#     class <- eco_codes$ecosystem_name[e] %>% 
#       make_clean_names()
#     ic_masked <- ic_lyr[[1]] %>%
#       terra::mask(mask = eco_mask)
#     names(ic_masked) <- paste0(og_name, "_", class)
#     add(ic_eco_stack) <- ic_masked
#     remove(eco_mask)
#     remove(ic_masked)
#     gc()
#   }
#   
#   # then save stack for that ic 
#   writeRaster(
#     ic_eco_stack, 
#     filename = paste0("data/irrecoverable_carbon/", og_name, "_ecosystem_stack.tif"),
#     overwrite = TRUE)
# }




#####

##### Population
#####

# Global mosaic of worldpop population count for 2020
# no prep needed



#####

##### Carbon Sequestration
#####

# ### Carbon sequestration potential
# # These are essentially Bernal coefficients spatially mapped by
# # Mariano. I was intending to update carbon sequestration potential
# # this year using Cook-Patton along with other datasets. 
# # Probably just use the carbon seq potl stack that exists already
# # in the data. However, note that Peru is randomly missing. 

# seq_files <- list.files("data/carbon_sequestration/",
#                         full.names = TRUE)
# 
# # define variables based on file names
# vars <- data.frame(vars = str_sub(
#   seq_files, start = 27, end = -30
# )) %>% 
#   distinct() %>% 
#   pluck("vars")
# 
# # loop through to mosaic per variable and add to stack
# for (i in seq_along(vars)){
#   var <- vars[i]
#   nos <- which(str_detect(seq_files, var))
#   mosaic <- mosaic(
#     rast(seq_files[nos[1]]),
#     rast(seq_files[nos[2]])
#   )
#   
#   if(i == 1){
#     seq_stack <- mosaic
#   } else {
#     add(seq_stack) <- mosaic
#   }
# }
# 
# names(seq_stack) <- vars
# 
# writeRaster(seq_stack, 
#             "data/carbon_sequestration/carbon_sequestration_potl_stack.tif", overwrite = TRUE)
# 

####################################################################

#####
##### Define Sites & Site Intersections/Overlaps
#####

### Edits for Feature Class ###

# The input file geodatabase
# Edit to reflect where your data is. 
fgdb <- "C:/Users/aballasiotes/Dev/Projects/ci_impact_indicators/data/ci_sites/CISites.gdb"

gdb_layers <- vect(fgdb)

fc <- vect(fgdb, layer = "FY2024_Draft_R_102524")

# Determine the FC extent, projection, and attribute information
summary(fc)

# View the feature class
plot(fc)

# Step 2: Convert SpatVector to sf object
sf_object <- st_as_sf(fc)

# Step 3: Check and fix invalid geometries
if (!all(st_is_valid(sf_object))) {
  sf_object <- st_make_valid(sf_object)
}

# Optional: Disable s2 geometry processing if errors persist
sf_use_s2(FALSE)

#########################
##      ALL DATA       ##
#########################
# clean, calculate area
# Remove historic sites!!
# Note, if pattern for Site Status changes for FY25,
# Will need to edit the filter here. 

fc_clean <- sf_object %>% 
  st_transform("epsg:5070") %>%
  st_zm() %>%  
  st_make_valid() %>%
  sf::st_buffer(dist = 0) %>% 
  lwgeom::st_snap_to_grid(50) %>% #Snaps the geometries to a grid with a cell size of 50 units.
  st_set_precision(50) %>% #Sets the precision of the geometries to 50 units.
  st_make_valid() %>%  #Makes the geometries valid again after snapping to the grid.
  sf::st_buffer(dist = 0) %>%  # Another buffer operation with a distance of 0, likely to further ensure validity.
  filter(Site_Status == "New this year" |
           Site_Status == "Continued this year" |
           Site_Status_Sequestration == "Continued this year" |
           Site_Status_Sequestration == "New this year") %>%
  filter(!is.na(CI_Start_Date)) %>%
  filter(Intervention_Type != "Not Applicable") %>%
  mutate(undr_rest = (Intervention_Category == 'Restoration Areas')) %>%
  rowwise() %>% 
  mutate("rest_area" = case_when(
    undr_rest == TRUE ~ Area_ha,
    T ~ 0
  )) %>% # Creates a new column rest_area based on conditions specified
  ungroup() %>%
  mutate(origin = row_number()) %>% # Adds a new column origin with the row numbers.
  dplyr::select(!c(
    GlobalID, CreationDate, Creator, EditDate, Editor, Shape_Area, Shape_Length
  )) %>%  #Drops selected columns from the dataset.
  filter(!st_is_empty(.)) #Removes rows with empty geometries.

shp_crs <- fc_clean %>%
  st_transform(my_crs) %>%
  st_make_valid() 

# Save to Shapefile so you can see how the sites were cleaned / changed
# etc. in GIS. 
write_sf(shp_crs, paste0("data/ci_sites/FY", year, "_Sites_Clean.shp"))

# Save to RDS so you can use it in R
saveRDS(fc_clean, paste0("data/ci_sites/FY", year, "_Sites_Clean.rds"))

# Save to a CSV to perform other edits more quickly that
# don't require spatial information.

fc_clean_df <- fc_clean %>% st_drop_geometry()
write_csv(fc_clean_df, paste0("data/ci_sites/FY", year, "_Sites_Clean.csv"))  


# # These sites were giving me issues, so I removed them for FY24 (and fixed them later)
# fc_clean_fy24 <- FY2024_Sites_Clean %>%
#   filter(!CI_ID %in% c("PHL1105", "PHL1106", "PHL1107", "PHL1108", "PHL1109", "PHL1110", "BNA1001", "BNA1002", "BNA1003", "BNA1017", "BNA1018", "BNA1019",
#                        "BNA1074", "BNA1075", "BNA1076", "BNA1077", "BNA1078", "BNA1079",
#                        "BNA1080", "BNA1081"))
# 
# 
# fc_clean <- rbind(fc_clean_fy24, fc_clean_added)


#########################
##      OVERLAPS       ##
#########################

## This is how we will tie overlaps to the original site - the "origin" column
fc_clean <- fc_clean %>%
  mutate(origin = row_number())

# calculate intersections
# start w/ cleaned fc
# I think Cam called these intersections
# because they possibly could've included self intersections
  
intersections <- fc_clean %>%  
  dplyr::select(origin, geometry) %>% 
  st_intersection() %>%  
  st_transform(my_crs) %>%
  st_make_valid()

saveRDS(intersections, paste0("data/ci_sites/FY", year, "_intersections.rds"))

# Keep the valid intersections
good_intersections <- intersections[st_is_valid(intersections),, drop = FALSE]

# From now on, we call them overlaps. I think Cam
# made this distinction originally to highlight that we are
# only taking good intersections where there is more than one site
# with overlaps. 

overlaps <- good_intersections %>%
  filter(n.overlaps > 1) %>%
  dplyr::select(origin, n.overlaps, origins, geometry) %>%
  rowid_to_column()

saveRDS(overlaps,paste0("data/ci_sites/FY", year, "_Overlaps.rds"))

#########################################################

# # Read in the sites if you're re-doing things.
# overlaps <- readRDS(paste0("data/ci_sites/FY", year, "_Overlaps.rds"))
# fc_clean <- readRDS(paste0("data/ci_sites/FY", year, "_Sites_Clean.rds"))

#########################################################


# link intersections to CI_ID

origin_ref <- fc_clean %>% 
  mutate(origin = row_number()) %>%
  sf::st_as_sf() %>%  # Convert to sf object
  st_transform(my_crs) %>% 
  st_zm() %>% 
  st_drop_geometry() %>% 
  dplyr::select(origin, CI_ID, Country, Division_Primary_CI, Division_Secondary_CI,
                CI_SLS, CI_Role, Intervention_Category, Geographic_Priority, Initial_Reporting_Year)


origin_ref$origin <- seq_len(nrow(fc_clean))

# At one point, I was having trouble with the overlaps being read properly.
# This step may not be necessary. 

st_agr(overlaps)
st_agr(overlaps) <- "constant"

overlaps <- overlaps %>% 
  rename(originOG = origin)

# This looks complicated, but it's essentially
# creating a table where every CI ID is matched (up to 6 overlaps...
# hopefully there will not be more than 6 overlaps...)


int_fc_clean <- overlaps %>%
  st_transform(my_crs) %>% 
  st_zm() %>% 
  filter(n.overlaps > 1) %>% #this is an unnecessary line
  add_column("geom_type" = st_geometry_type(.)) %>% 
  dplyr::filter(geom_type == "POLYGON" | geom_type == "MULTIPOLYGON") %>% 
  add_column("Area_ha_R" = as.numeric(st_area(.))/10000) %>%
  filter(Area_ha_R > 0.25) %>%
  rename(origin = origins) %>%
  dplyr::select(c(
     n.overlaps, origin, Area_ha_R,
   )) %>%
  rowwise() %>%
  rename(origins = origin) %>%
  mutate("origin" = origins[1]) %>%
  left_join(origin_ref, by = "origin")  %>%
  rename(CI_ID_1 = CI_ID) %>%
  mutate("origin" = origins[2])  %>%
  left_join(origin_ref, by = "origin", suffix=c("",".y"))%>%
  dplyr::select(-ends_with(".y"))  %>%
  rename(CI_ID_2 = CI_ID) %>%
  mutate("origin" = case_when(length(origins) > 2 ~ origins[3])) %>%
  left_join(origin_ref, by = "origin", suffix=c("",".y"))%>%
  dplyr::select(-ends_with(".y")) %>%
  rename(CI_ID_3 = CI_ID) %>%
  mutate("origin" = case_when(length(origins) > 3 ~ origins[4])) %>%
  left_join(origin_ref, by = "origin", suffix=c("",".y"))%>%
  dplyr::select(-ends_with(".y")) %>%
  rename(CI_ID_4 = CI_ID) %>%
  mutate("origin" = case_when(length(origins) > 4 ~ origins[5])) %>%
  left_join(origin_ref, by = "origin", suffix=c("",".y"))%>%
  dplyr::select(-ends_with(".y")) %>%
  rename(CI_ID_5 = CI_ID) %>%
  mutate("origin" = case_when(length(origins) > 5 ~ origins[5])) %>%
  left_join(origin_ref, by = "origin", suffix=c("",".y"))%>%
  dplyr::select(-ends_with(".y")) %>%
  rename(CI_ID_6 = CI_ID) 
  # dplyr::select(!origin) %>%
  # ungroup()

saveRDS(int_fc_clean, paste0("data/ci_sites/FY", year, "_Overlaps_Clean.rds"))


