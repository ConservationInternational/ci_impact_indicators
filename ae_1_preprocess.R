library(tidyverse)
library(raster)
library(terra)
library(sf)
library(geodata)

template <- rast("~/Conservation International/Data/land_1km_eck4.tif")

# SRTM30+ Global 1-km Digital Elevation Model (DEM): Version 11: Land Surface
# from https://pae-paha.pacioos.hawaii.edu/erddap/griddap/srtm30plus_v11_land.html

# need to mosaic tiles
dem_files <- list.files("data/avoided_emissions", pattern = "srtm30",
                        full.names = TRUE)
for (i in seq_along(dem_files)){
  rast <- rast(dem_files[i])
  if(i == 1){dem <- rast} else {dem <- mosaic(dem, rast)}
}

writeRaster(
  dem,
  "data/avoided_emissions/SRTM_DEM.tif")

# project to be able to calculate slope

dem_proj <- dem %>% 
  project(template, 
          method = "bilinear",
          threads = TRUE,
          mask = TRUE)

writeRaster(
  dem_proj,
  "data/avoided_emissions/covariate_dem.tif")

# create slope
slope <- terrain(dem_proj, "slope")

writeRaster(slope,
            "data/avoided_emissions/covariate_slope.tif")



# suitability probability for cropland
# for now, mosaic from Chen

crop_suit <- rast(
  "data/avoided_emissions/croplandSuitabilityProbability_Chen2022.tif") %>% 
  project(template, 
          method = "bilinear",
          threads = TRUE,
          mask = TRUE)

writeRaster(
  crop_suit,
  "covariate_cropland_suitability.tif")

precip <- rast("data/avoided_emissions/chirps-v2.0.1981-2020.40yrs.tif") %>% 
  project(template, 
          method = "bilinear",
          threads = TRUE) %>% 
  extend(template) %>% 
  crop(template) %>% 
  mask(slope)

writeRaster(
  precip,
  "covariate_precipitation.tif",
  overwrite = TRUE)

# temperature
# have mean per month... take average of all

temp_files <- list.files("data/avoided_emissions/wc2.1_30s_tavg/", 
                         pattern = ".tif", full.names = TRUE)

for (i in seq_along(temp_files)){
  rast <- rast(temp_files[i])
  if(i == 1) {temp_stack <- rast} else {temp_stack <- c(temp_stack, rast)}
}

temp <- mean(temp_stack) %>% 
  project(template, 
          method = "bilinear",
          threads = TRUE,
          mask = TRUE)

writeRaster(
  temp,
  "data/avoided_emissions/covariate_temp.tif")

# cropland
crop <- rast("data/avoided_emissions/Global_cropland_3km_2019.tif") %>% 
  project(template, 
          method = "bilinear",
          threads = TRUE,
          mask = TRUE)

writeRaster(
  crop,
  "data/avoided_emissions/covariate_cropland.tif")


# accessibility to cities from Weiss et al 
cities <- rast("data/avoided_emissions/accessibility_to_cities_2015_v1.0.tif") %>% 
  project(template, 
          method = "bilinear",
          threads = TRUE,
          mask = TRUE)

writeRaster(
  cities,
  "data/avoided_emissions/covariate_cities.tif")


wdpa_poly <- read_sf("data/avoided_emissions/global-2023-04-10.gpkg") %>%
  janitor::clean_names() %>% 
  st_transform(crs(template)) %>% 
  filter(!marine == "marine") %>% 
  vect() 

wdpa <- terra::rasterize(
  x = wdpa_poly,
  y = template,
  cover = FALSE)

writeRaster(
  wdpa,
  "data/avoided_emissions/covariate_wdpa.tif")

# distance to roads
# from https://www.globio.info/download-grip-dataset
# then calculated eucdistance with template as template in ArcGIS

roads <- rast("data/avoided_emissions/GRIP4roads_eucdistance.tif") %>% 
  project(template, 
          method = "bilinear",
          threads = TRUE,
          mask = TRUE)

writeRaster(
  roads, 
  "data/avoided_emissions/covariate_distroads.tif")

# administrative unit
# per GADM 1st sub-national level

gadm <- read_sf("data/avoided_emissions/gadm_410.gpkg") %>% 
  dplyr::select(NAME_1) %>% 
  st_transform(crs = crs(template)) %>% 
  vect() %>% 
  rasterize(y = template,
            field = "NAME_1")

writeRaster(
  gadm, 
  "data/avoided_emissions/covariate_gadm.tif")

# ecoregion

eco <- read_sf(
  dsn = "data/avoided_emissions",
  layer = "Ecoregions2017") %>% 
  dplyr::select(ECO_NAME) %>% 
  st_transform(crs(template)) %>% 
  vect() %>% 
  rasterize(y = template,
            field = "ECO_NAME")

writeRaster(
  eco,
  "data/avoided_emissions/covariate_ecoregions.tif")

biome <- read_sf(
  dsn = "data/avoided_emissions",
  layer = "Ecoregions2017") %>% 
  dplyr::select(BIOME_NAME) %>% 
  st_transform(crs(template)) %>% 
  vect() %>% 
  rasterize(y = template,
            field = "BIOME_NAME")

writeRaster(
  biome,
  "data/avoided_emissions/covariate_biome.tif")


# Population for 2000, 2005, 2010, 2015 2020 from 
# https://hub.worldpop.org/geodata/listing?id=64

mask_rast <- rast("data/avoided_emissions/covariate1_slope.tif")

for(i in 1:5){
  rast <- rast(paste0("data/population/ppp_", 
                      seq(from = 2000, to = 2020, by = 5)[i],
                      "_1km_Aggregated.tif")) %>%   
    project(template, 
            method = "bilinear",
            threads = TRUE) %>% 
    mask(mask_rast)
  
  writeRaster(rast, 
              paste0("data/population/ppp_", 
                     seq(from = 2000, to = 2020, by = 5)[i],
                     "_1km_Aggregated_proj.tif"),
              overwrite = TRUE)
  
  if(i == 1) {population <- rast } else {population <- c(population, rast)}
}

population <- raster::stack(population)

names(population) <-
  c('pop_2000', 'pop_2005', 'pop_2010', 'pop_2015', 'pop_2020')
NAvalue(population) <- -32768

population_growth <-
  overlay(
    population$pop_2000,
    population$pop_2020,
    fun = function(pop_2000, pop_2020) {
      r = ((pop_2020 / pop_2000) ^ (1 / 15) - 1) * 100
      # Set growth to zero when 2000 and 2020
      # populations are equal
      r[pop_2000 == pop_2020] <- 0
      # Set growth rate to maximum of other
      # pixels for areas that grew from zero
      # population in 2000.
      isinfs <- is.infinite(r)
      r[isinfs] <- 0
      r[isinfs] <- max(r, na.rm = TRUE)
      r
    }
  )

writeRaster(
  population_growth,
  filename = 'data/avoided_emissions/covariate_population_growth.tif',
  overwrite = TRUE,
  options = "COMPRESS=LZW",
  datatype = "INT2S"
)


# biomass

biomass <- rast('data/carbon_stored/biomass_prepped_2022.tif') %>%   
  project(template, 
          method = "bilinear",
          threads = TRUE) %>% 
  mask(mask_rast)

writeRaster(
  biomass,
  "data/avoided_emissions/covariate_biomass.tif")


# land classes
# from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8349160/ 1km agg
# downloaded from https://figshare.com/articles/dataset/Open_Synergy_Land_Cover_Series/14329025/1?file=27322868
# information here https://data.apps.fao.org/map/catalog/srv/eng/catalog.search#/metadata/b68a0adb-06da-44bb-a5aa-7d0ca95b654e
# chose because closest to land cover classes needed, and already agg to 1km

lc <- rast("data/avoided_emissions/OS2015.tif")

names(lc) <- c(
  "cropland",
  "forest",
  "grassland",
  "shrubland",
  "wetland",
  "water",
  "tundra",
  "impervious surface",
  "bareland")

list <- as.list(lc)

forest <- list[[2]]
grassland <- list[[3]]
agriculture <- list[[1]]
wetlands <- list[[5]]
artificial <- list[[8]]
other <- max(list[[4]], list[[7]], list[[9]], na.rm = TRUE)
water <- list[[6]]

lc_cov <- c(forest, 
            grassland,
            agriculture,
            wetlands,
            artificial, 
            other,
            water) %>% 
  project(template, 
          method = "bilinear",
          threads = TRUE,
          mask = TRUE)

writeRaster(
  lc_cov,
  "data/avoided_emissions/covariate_lc2015.tif")



lc <- rast("data/avoided_emissions/OS2000.tif")

names(lc) <- c(
  "cropland",
  "forest",
  "grassland",
  "shrubland",
  "wetland",
  "water",
  "tundra",
  "impervious surface",
  "bareland")

list <- as.list(lc)

forest <- list[[2]]
grassland <- list[[3]]
agriculture <- list[[1]]
wetlands <- list[[5]]
artificial <- list[[8]]
other <- max(list[[4]], list[[7]], list[[9]], na.rm = TRUE)
water <- list[[6]]

lc_cov <- c(forest, 
            grassland,
            agriculture,
            wetlands,
            artificial, 
            other,
            water) %>% 
  project(template, 
          method = "bilinear",
          threads = TRUE,
          mask = TRUE)

writeRaster(
  lc_cov,
  "data/avoided_emissions/covariate_lc2000.tif")


# hansen forest cover prepped in hansen_prep