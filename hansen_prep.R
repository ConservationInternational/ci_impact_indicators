# Hansen data download and prep
# Note this script takes ~4 days to run

library(tidyverse)
library(raster)
library(terra)
library(sf)
library(gfcanalysis)

template <- rast("~/Conservation International/Data/land_1km_eck4.tif")
hansen_folder <- "misc/hansen_fc_tiles/"


# only have to do once
# sites_poly <- read_sf(
#   dsn = "data/ci_sites",
#   layer = "FY2022_Sites") %>%
#   filter(biome != "Marine")
# hansen_tiles <- calc_gfc_tiles(sites_poly)
# download_tiles(
#   tiles = hansen_tiles,
#   output_folder = hansen_folder,
#   images = c("treecover2000", "lossyear"),
#   dataset = "GFC-2022-v1.10")
# # get rid of tiles that are just ocean
# tile_names <- list.files(hansen_folder,
#                          pattern = "lossyear") %>% 
#   str_sub(., start = 32, end = -5) 
# 
# for(i in seq_along(tile_names)){
#   
#   print(paste0(i, "/", length(tile_names)))
#   
#   tile_name <- tile_names[[i]]
#   
#   tc_2000 <- rast(paste0(
#     hansen_folder, "Hansen_GFC-2022-v1.10_treecover2000_", tile_name, ".tif")) 
#   
#   # remove tile if no values in it
#   if(global(tc_2000, "max", na.rm = TRUE)$max == 0) { 
#     file.remove(paste0(
#       hansen_folder, "Hansen_GFC-2022-v1.10_treecover2000_", tile_name, ".tif")) 
#     file.remove(paste0(
#       hansen_folder, "Hansen_GFC-2022-v1.10_lossyear_", tile_name, ".tif"))
#     
#     print(paste0("...removed"))
#     }
# }

tile_names <- list.files(hansen_folder,
                         pattern = "lossyear") %>% 
  str_sub(., start = 32, end = -5) 

start_time <- Sys.time() # get start time

# years 2009 to 2022

for(y in 1:22){
  
  year <- case_when(y < 10 ~ paste0("200", y), 
                    T ~ paste0("20", y))
  
  rast_list <- list()
  
  for(i in seq_along(tile_names)){
    
    print(paste0("Year ", year, "... ", i, "/", length(tile_names)))
    
    tile_name <- tile_names[[i]]
    
    tc_2000 <- rast(paste0(
      hansen_folder, "Hansen_GFC-2022-v1.10_treecover2000_", tile_name, ".tif")) %>%
      aggregate(fact = 30, 
                fun = "mean",
                cores = 8,
                na.rm = TRUE)
    
    ly <- rast(paste0(
      hansen_folder, "Hansen_GFC-2022-v1.10_lossyear_", tile_name, ".tif")) %>% 
      classify(
        matrix(data = c(
          0, y, 1,
          y, 23, 0),
          ncol = 3,
          byrow = TRUE),
        include.lowest = FALSE) %>% 
      aggregate(fact = 30, 
                fun = "mean",
                cores = 8,
                na.rm = TRUE) * 100
    
    # remove forest from 2000 treecover dataset from loss up to year y
    tc_y <- tc_2000 - ly 
    tc_y[tc_y < 0] <- 0
    
    rast_list[i] <- tc_y
    
    # if(i == 1){
    #   ty_mosaic <- tc_y
    #   
    # } else {
    #   ty_mosaic <- merge(ty_mosaic, tc_y)
    # }
    
    gc()
    tmpFiles(remove = TRUE)
    
  }
  
  mosaic <- merge(sprc(rast_list))
  
  mosaic_proj <- mosaic %>% 
    project(template,
            align = TRUE,
            method = "bilinear",
            threads = TRUE) %>% 
    extend(template) %>% 
    crop(template)
  
  names(mosaic_proj) <- paste0("fc_", str_sub(year, start = 3))
  
  writeRaster(
    mosaic_proj, 
    paste0("avoided_emissions/data/covariate_forest_cover_", year, ".tif"),
    overwrite = TRUE)
  
  gc()
  
  # print elapsed time
  print(paste0("Start time: ", start_time))
  print(paste0("Current time: ", Sys.time()))
}




# now get forest change from those outputs

for(y in 1:22){
  
  year <- as.numeric(case_when(y < 10 ~ paste0("200", y), 
                    T ~ paste0("20", y)))
  
  print(year)
  
  year_before <- year - 1
  
  mosaic <-  rast( 
    paste0("avoided_emissions/data/covariate_forest_cover_", year, ".tif"))
  
  mosaic_before <- rast( 
    paste0("avoided_emissions/data/covariate_forest_cover_", year_before, ".tif"))
  
  change <- mosaic - mosaic_before
  
  
  writeRaster(
    change,
    paste0("avoided_emissions/data/covariate_forest_cover_change_", year, ".tif"),
    overwrite = TRUE)
}
