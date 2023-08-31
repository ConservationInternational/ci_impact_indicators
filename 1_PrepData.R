# Data preparation for impact indicators analysis
# Contact: Cameryn Brock
# Last updated: 3/19/2022

library(tidyverse)
library(raster)
library(terra)
library(sf)
library(janitor)

year <- "2022"

#####

##### Carbon Stored
#####

## Woody biomass
# Global forest watch Aboveground Live Woody Biomass stock in megagrams/pixel (tonnes/pixel)
gfw_files <- list.files("misc/woody_tiles_30m",
                        full.names = TRUE)
tiles <- str_sub(gfw_files, start = 22, end = -5)

hansen_dest30 <- "misc/hansen_tiles_30m/"

# read in each tile for hansen and gfw, calculate bgb, mask biomass to hansen, aggregate, save

rast_list <- c()


for (i in seq_along(tiles)){
  tile <- tiles[i]
  
  # above-ground biomass from gfw and bgb calculated from that
  gfw_file <- gfw_files[i]
  gfw_agb <- rast(gfw_file)
  bgb <- 0.489 * gfw_agb^(0.89) # calculate below-ground biomass based on Mokany formula
  bio <- sum(gfw_agb, bgb, na.rm = TRUE) 
  # add above-ground and below-ground for total biomass)
  # convert to carbon equilavent (*0.5) to get Total Carbon 
  ## NOTE instead of * 0.5 I divde by 2 in the B: extract data script
  
  # forest loss from hansen
  hansen_url <- paste0(case_when(
    year == "2022" ~
      "https://storage.googleapis.com/earthenginepartners-hansen/GFC-2021-v1.9/Hansen_GFC-2021-v1.9_lossyear_",
    year == "2021" ~
      "https://storage.googleapis.com/earthenginepartners-hansen/GFC-2020-v1.8/Hansen_GFC-2020-v1.8_lossyear_"), 
    tile, ".tif")
  download.file(
    url = hansen_url, 
    destfile = paste0(hansen_dest30, "hansen_", tile, ".tif"),
    mode = "wb")
  
  hansen <- rast(paste0(hansen_dest30, "hansen_", tile, ".tif"))
  
  # mask out hansen forest loss from biomass
  bio_masked <- bio %>% 
    terra::mask(hansen, inverse = TRUE, maskvalues = c(NA, 0))
  
  # aggregate to 300m
  bio_300 <- bio_masked %>% 
    terra::aggregate(fact = 10, fun = "sum", na.rm = TRUE)
  
  # add to list
  rast_list[[i]] <- bio_300
  
  tmpFiles(remove = TRUE)
}


# mosaic all together
rsrc <- terra::src(rast_list)
carbon_bio_mosaic <- mosaic(rsrc)

writeRaster(carbon_bio_mosaic, paste0("data/carbon_stored/biomass_prepped_", year, ".tif"))

# need to test with sites overlaid

## Soil biomass

# soil depth 1-3 for 30cm
# ocstha = organic soil carbon tonnes/ha

sd1 <- rast("data/carbon_stored/OCSTHA_M_sd1_1km_soc.tif")
sd2 <- rast("data/carbon_stored/OCSTHA_M_sd2_1km_soc.tif")
sd3 <- rast("data/carbon_stored/OCSTHA_M_sd3_1km_soc.tif")
# crs(sd3) == my_crs # TRUE

soil <- (sd1 + sd2 + sd3) * 100

writeRaster(soil, "data/carbon_stored/OCSTHA_30cm_1km.tif")





#####

##### Irrecoverable Carbon
#####

eco <- rast("data/ecosystems/ecosystems.tif") 

# total tonnes of carbon per hectare
total_ic <- rast("data/irrecoverable_carbon/Irrecoverable_C_Total_2018.tif") %>% 
  resample(eco, method = "bilinear")
# checked multiple methods of resampling + combining with area and this was the most accurate
# at achieving 139B tonnes

# convert to tonnes per pixel
area <- raster::area(raster(total_ic)) * 100 #km2 to ha

tonnes_ic <- raster(total_ic) * area

global(rast(tonnes_ic), "sum", na.rm = TRUE)

writeRaster(tonnes_ic, "tonnes_ic_bilmethod.tif")

total_ic <- rast("tonnes_ic_bilmethod.tif") 

# ic total with vals below 25 tonnes/ha as NA
high_ic <- rast("data/irrecoverable_carbon/ic_above_25.tif") %>% 
  resample(eco, method = "bilinear")

# ic_total w/ vals below 25 as NA

# make raster of 'any' ic (above 0.01)
any_ic <- total_ic
any_ic[any_ic > 0.01] <- 1
any_ic[any_ic <= 0.01] <- NA


# change ic_high to binary
high_ic[high_ic > 0] <- 1

# need 'high' and 'any' to reflect has/pixel
# change ic_high and ic_any to represent area

high_area <- cellSize(high_ic, unit = "ha")

any_area <- cellSize(any_ic, unit = "ha")


# create stack
ic_stack <- c(total_ic, high_area, any_area)  

names(ic_stack) <- c("tstor_ic", "ha_high_ic", "ha_ic")

writeRaster(ic_stack, "data/irrecoverable_carbon/ic_stack_prepped.tif",
            overwrite = TRUE)

ic_stack <- rast("data/irrecoverable_carbon/ic_stack_prepped.tif")

# check total to confirm it matches with our 140B tonnes total in Noon et al
global(ic_stack$tstor_ic, "sum", na.rm = TRUE)

# divvy up irrecoverable carbon by ecosystem

### ecosystems to mask irr carbon to

eco <- rast("data/ecosystems/ecosystems.tif") 
eco_codes <- data.frame(cats(eco$class)) %>%
  rename('ecosystem_name' = class) %>% 
  drop_na()

eco_seg <- eco %>%
  terra::segregate()

eco_seg[eco_seg == 0] <- NA

writeRaster(eco_seg, "eco_seg_nas.tif", overwrite = TRUE)

eco_seg <- rast("eco_seg_nas.tif")

# for each ic layer
for (i in 1:nlyr(ic_stack)){
  
  ic_lyr <- ic_stack[[i]]
  og_name <- names(ic_lyr)
  ic_eco_stack <- ic_lyr
  
  # mask irr carbon to each ecosystem layer then add to stack
  for(e in 1:nlyr(eco_seg)){
    eco_mask <- eco_seg[[e]]
    class <- eco_codes$ecosystem_name[e] %>% 
      make_clean_names()
    ic_masked <- ic_lyr[[1]] %>%
      terra::mask(mask = eco_mask)
    names(ic_masked) <- paste0(og_name, "_", class)
    add(ic_eco_stack) <- ic_masked
    remove(eco_mask)
    remove(ic_masked)
    gc()
  }
  
  # then save stack for that ic 
  writeRaster(
    ic_eco_stack, 
    filename = paste0("data/irrecoverable_carbon/", og_name, "_ecosystem_stack.tif"),
    overwrite = TRUE)
}




#####

##### Population
#####

# Global mosaic of worldpop population count for 2020
# no prep needed



#####

##### Carbon Sequestration
#####

### Carbon sequestration potential
seq_files <- list.files("data/carbon_sequestration/",
                        full.names = TRUE)

# define variables based on file names
vars <- data.frame(vars = str_sub(
  seq_files, start = 27, end = -30
)) %>% 
  distinct() %>% 
  pluck("vars")

# loop through to mosaic per variable and add to stack
for (i in seq_along(vars)){
  var <- vars[i]
  nos <- which(str_detect(seq_files, var))
  mosaic <- mosaic(
    rast(seq_files[nos[1]]),
    rast(seq_files[nos[2]])
  )
  
  if(i == 1){
    seq_stack <- mosaic
  } else {
    add(seq_stack) <- mosaic
  }
}

names(seq_stack) <- vars

writeRaster(seq_stack, 
            "data/carbon_sequestration/carbon_sequestration_potl_stack.tif", overwrite = TRUE)


#####

##### Define Sites & Site Intersections/Overlaps
#####

my_crs <- crs(rast("data/irrecoverable_carbon/tstor_ic_ecosystem_stack.tif"))

sf_use_s2(FALSE)

# clean, calculate area, 
# remove wwf, bna, swl, and proposed sites per request
shp <- read_sf(
  dsn = "data/ci_sites",
  layer = case_when(year == "2022" ~ "FY22_Vetting")) %>% 
  clean_names() %>% 
  st_transform("epsg:5070") %>%
  st_zm() %>% 
  st_make_valid() %>% 
  sf::st_buffer(dist = 0) %>%
  lwgeom::st_snap_to_grid(50) %>% 
  st_set_precision(50) %>% 
  st_make_valid() %>%  
  sf::st_buffer(dist = 0) %>%
  filter(!str_detect(ci_id, "WWF"),
         !str_detect(ci_id, "BNA"),
         !str_detect(ci_id, "SLW")) %>%
  filter(!interventi == "Protected Area - Proposed (National or Regional)") %>% 
  mutate(undr_rest = (restoratio != ' ' & 
                        restoratio != 'Not Applicable' & 
                        !is.na(restoratio)))  %>% 
  mutate(area_ha = as.numeric(round(area_ha, digits = 2))) %>% 
  rowwise() %>% 
  mutate("rest_area" = case_when(
    undr_rest == TRUE ~ area_ha,
    T ~ 0
  )) %>% 
  ungroup() %>%
  mutate(origin = row_number()) %>% 
  dplyr::select(!c(
    global_id, creation_da, creator, edit_date, editor,
    internal_e, shape_area, shape_leng, undr_rest
  )) %>% 
  filter(!st_is_empty(.))

shp_crs <- shp %>% 
  st_transform(my_crs) %>% 
  st_make_valid() 

write_sf(shp_crs, paste0("data/ci_sites/FY", year, "_Sites.shp"))

# calculate intersections
# start w/ original shapefiles with sites requested removed

intersections <- shp %>%  
  dplyr::select(origin, geometry) %>% 
  st_intersection() %>%  
  st_transform(my_crs) %>% 
  st_make_valid()  

int_shp <- intersections %>% 
  filter(n.overlaps > 1) %>% 
  dplyr::select(!origins) %>% 
  st_collection_extract("POLYGON")

write_sf(int_shp, paste0("data/ci_sites/FY", year, "_Intersections.shp"))

write_sf(dplyr::select(intersections, !origins), 
         paste0("data/ci_sites/FY", year, "_Intersections.shp"))

int <- intersections %>% 
  filter(n.overlaps > 1) %>% 
  dplyr::select(n.overlaps, origins, geometry) %>% 
  mutate(row_number = row_number())

saveRDS(int,paste0("data/ci_sites/FY", year, "_Overlaps.rds"))

# link intersections to ci_id
origin_ref <- read_sf(
  dsn = "data/ci_sites",
  layer = paste0("FY", year, "_Sites")
) %>% 
  mutate(origin = row_number()) %>% 
  st_transform(my_crs) %>% 
  st_zm() %>% 
  clean_names() %>%
  st_drop_geometry() %>% 
  dplyr::select(origin, ci_id, country, ci_divisio, ci_divis_1, ci_sls_1, ci_sls_2)
# want country, division, sls, site

shp_df <- shp %>% st_drop_geometry()

int <- readRDS(paste0("data/ci_sites/FY", year, "_Overlaps.rds")) %>%
  clean_names() %>% 
  st_transform(my_crs) %>% 
  st_zm() %>% 
  filter(n_overlaps > 1) %>% 
  mutate("geom_type" = st_geometry_type(.)) %>% 
  filter(geom_type == "POLYGON" | geom_type == "MULTIPOLYGON") %>% 
  mutate("area_ha" = as.numeric(st_area(.))/10000) %>% 
  filter(area_ha > 0.25) %>%
  rename(origin = origins) %>% 
  dplyr::select(c(
    n_overlaps, origin, area_ha, 
  )) %>% 
  rename(origins = origin) %>% 
  rowwise() %>% 
  mutate(origin = origins[1]) %>% 
  left_join(origin_ref, by = "origin") %>% 
  rename(ci_id_1 = ci_id) %>% 
  mutate(origin = origins[2]) %>% 
  left_join(origin_ref, by = "origin") %>% 
  rename(ci_id_2 = ci_id) %>% 
  mutate(origin = case_when(length(origins) > 2 ~ origins[3])) %>% 
  left_join(origin_ref, by = "origin") %>% 
  rename(ci_id_3 = ci_id) %>% 
  mutate(origin = case_when(length(origins) > 3 ~ origins[4])) %>% 
  left_join(origin_ref, by = "origin") %>% 
  rename(ci_id_4 = ci_id) %>% 
  mutate(origin = case_when(length(origins) > 4 ~ origins[5])) %>% 
  left_join(origin_ref, by = "origin") %>% 
  rename(ci_id_5 = ci_id) %>% 
  dplyr::select(!origin) %>% 
  ungroup()

saveRDS(int, paste0("data/ci_sites/FY", year, "_Overlaps_Clean.rds"))
