library(sf)
library(tidyverse)
library(raster)
library(terra)
library(fasterize)
library(foreach)

data_folder <- 'data/avoided_emissions'
template <- rast(paste0(data_folder, "/land_1km_eck4.tif"))


load_as_vrt <- function(folder, pattern, band=FALSE, raster=TRUE) {
    vrt_file <- tempfile(fileext='.vrt')
    files <- list.files(folder, pattern=pattern, full.names = TRUE)
    if (length(files) == 0) {
        stop('No files found')
    }
    if (band) {
        gdalbuildvrt(paste0(folder, '/', files), vrt_file, b=band)
        r <- raster(vrt_file)
    } else {
        gdalbuildvrt(paste0(folder, '/', files), vrt_file)
        r <- raster::stack(vrt_file)
    }
    if (raster) {
        return(r)
    } else {
        return(vrt_file)
    }
}

# Function used to get IDs from a rasterized set of polygons (to determine 
# which polygons were lost due to rasterization (very small polygons drop out)
get_unique <- function(x) {
    bs <- raster::blockSize(x)
    n_blocks <- bs$n
    for (block_num in 1:n_blocks) {
        these_vals <- unique(raster::getValues(x,
                                               row=bs$row[block_num], 
                                               nrows=bs$nrows[block_num]))
        if (block_num == 1) {
            out <- these_vals
        } else {
            out <- unique(c(out, these_vals))
        }
    }
    return (na.omit(out))
}

###############################################################################
fp_files <- list.files(data_folder, pattern="*.tif$", full.names = TRUE)

# Load covariates
covariates_1 <- raster::stack(
    list.files(file.path(data_folder), pattern = "^covariate1.*\\.tif$", full.names = TRUE))


names(covariates_1) <- c('biome',
                         'elev',
                         'ecoregion',
                         'precip',
                         'slope',
                         'temp')
NAvalue(covariates_1) <- -32768

covariates_2 <- raster::stack(
    list.files(file.path(data_folder), pattern = "^covariate2.*\\.tif$", 
               full.names = TRUE))
names(covariates_2) <- c('dist_cities',
                         'crop_suitability',
                         'dist_roads',
                         'pa')
NAvalue(covariates_2) <- -32768

population <- raster::stack(
    list.files("data/population", pattern = "_proj.tif", full.names = TRUE))
names(population) <-
    c('pop_2000', 'pop_2005', 'pop_2010', 'pop_2015', 'pop_2020')
NAvalue(population) <- -32768

biomass <- raster('data/avoided_emissions/covariate_biomass_2024.tif')
names(biomass) <- c('total_biomass')

population_growth <- raster(
    'data/avoided_emissions/covariate_population_growth.tif')
names(population_growth) <- c('pop_growth')

covariates <- stack(covariates_1, covariates_2, biomass, population, 
                    population_growth)
writeRaster(
    covariates, 
    filename='data/avoided_emissions/covariates_covariates_2024.tif', 
    overwrite=TRUE, options="COMPRESS=LZW", datatype="INT2S")

write_csv(data.frame(names=names(covariates)), 'data/avoided_emissions/covariates_covariates.csv')

lc_2000 <- raster::stack("data/avoided_emissions/covariate_lc2000.tif")
names(lc_2000) <- c('lc_2000_forest',
                    'lc_2000_grassland',
                    'lc_2000_agriculture',
                    'lc_2000_wetlands',
                    'lc_2000_artificial',
                    'lc_2000_other',
                    'lc_2000_water')
writeRaster(lc_2000, filename='data/covariates_lc_2000.tif', 
            overwrite=TRUE, options="COMPRESS=LZW", datatype="INT2S")

lc_2015 <- raster::stack("data/avoided_emissions/covariate_lc2015.tif")
names(lc_2015) <- c('lc_2015_forest',
                    'lc_2015_grassland',
                    'lc_2015_agriculture',
                    'lc_2015_wetlands',
                    'lc_2015_artificial', 
                    'lc_2015_other',
                    'lc_2015_water')
writeRaster(lc_2015, filename='data/avoided_emissions/covariates_lc_2015.tif', 
            overwrite=TRUE, options="COMPRESS=LZW", datatype="INT2S")
write_csv(data.frame(names=names(lc_2015)), 
          'data/avoided_emissions/covariates_lc_2015.csv')


# prep Forest Cover & Forest Cover Change using GEE script below:
# https://code.earthengine.google.com/81280861c10c98478d9e3fa7b2d6e65c
# if in later scripts you are having projection issues, there are some
# code snippets below that take care of that. 

# forest cover

fc <- raster::stack(list.files(data_folder, pattern = "^covariate_forest_cover_fc.*\\.tif$",,
                               full.names = TRUE))
NAvalue(fc) <- -32768

fcnames_list <- c("fc2000", "fc2001", "fc2002", "fc2003", "fc2004", "fc2005", "fc2006", "fc2007", "fc2008", "fc2009", "fc2010", "fc2011", "fc2012", "fc2013", "fc2014", "fc2015", "fc2016", "fc2017", "fc2018", "fc2019", "fc2020", "fc2021", "fc2022", "fc2023")
names(fc) <- fcnames_list
print(fc)


## SAVE IT
writeRaster(fc, filename='data/avoided_emissions/covariates_fc.tif', 
            overwrite=TRUE, options="COMPRESS=LZW", datatype="INT2S")

## WRITE NAMES TO CSV
write_csv(data.frame(names=names(fc)), 'data/avoided_emissions/covariates_fc.csv')


# forest cover change
fcc <- raster::stack(list.files(data_folder, pattern = "^covariate_forest_cover_change_fcc_20.*\\.tif$",,
                               full.names = TRUE))
NAvalue(fcc) <- -32768

fcc_names_list <- c("fcc_2001", "fcc_2002", "fcc_2003", "fcc_2004", "fcc_2005", "fcc_2006", "fcc_2007", "fcc_2008", "fcc_2009", "fcc_2010", "fcc_2011", "fcc_2012", "fcc_2013", "fcc_2014", "fcc_2015", "fcc_2016", "fcc_2017", "fcc_2018", "fcc_2019", "fcc_2020", "fcc_2021", "fcc_2022", "fcc_2023")
print(fcc_names_list)
names(fcc) <- fcc_names_list
print(fcc)

writeRaster(fcc, filename='data/avoided_emissions/covariates_fc_change.tif', 
            overwrite=TRUE, options="COMPRESS=LZW", datatype="INT2S")

write_csv(data.frame(names=names(fcc)), 'data/avoided_emissions/covariates_fc_change.csv')

fc <- rast('data/avoided_emissions/covariates_forcov.tif')
fcc <- rast('data/avoided_emissions/covariates_fc_change.tif')

###############################################################################
### Load GADM boundaries
regions <- st_read("data/avoided_emissions/gadm_410-levels.gpkg"
                   , layer = "ADM_1") %>% 
    st_transform(crs(template))
regions$level0_ID <- as.numeric(factor(regions$GID_0))
regions$level1_ID <- as.numeric(factor(regions$GID_1))
regions_rast <- fasterize(regions, raster(template),
                          field = 'level1_ID')
region_IDs_after_rasterization <- get_unique(regions_rast)
regions <- regions[regions$level1_ID %in% region_IDs_after_rasterization, ]

regions$level0_ID <- as.numeric(factor(as.character(regions$GID_0)))
regions$level1_ID <- as.numeric(factor(as.character(regions$GID_1)))
# Now re-rasterize boundaries (with ID's that will disappear dropped) to ensure
# that all IDs are sequential and that they match between the data.frame and 
# the raster.
regions_rast <- fasterize(regions, raster(template),
                          field='level1_ID')

names(regions_rast) <- 'region'
region_IDs_after_rasterization <- get_unique(regions_rast)
stopifnot(sort(region_IDs_after_rasterization) == sort(regions$level1_ID))
saveRDS(regions, file='data/avoided_emissions/regions.RDS')


##############
### TAKE GEE FILES & REPROJECT
###########

fc_files <- list.files('data/avoided_emissions/', pattern="^covariate_forest_cover_fc.*\\.tif$", full.names = TRUE)
print(fc_files)

fcc_files <- list.files('data/avoided_emissions/', pattern="^covariate_forest_cover_change_fcc.*\\.tif$", full.names = TRUE)
print(fcc_files)


##
# Step 2: Loop over the list of files and reproject each
for (file in fcc_files[1:2]) {
  
  # Read the raster data
  raster_data <- rast(file)
  
  reprojected_data <- raster_data %>% 
    project(template,
            align = TRUE,
            method = "bilinear",
            threads = TRUE) %>% 
    extend(template) %>% 
    crop(template)
  
  # Step 3: Generate a new file name (e.g., append "_reprojected" to the original name)
  new_file <- gsub("\\.tif$", "_p.tif", file)
  
  # Save the reprojected file
  writeRaster(reprojected_data, new_file, overwrite = TRUE)
  
  # Optionally, print the new file path to track progress
  print(paste("Reprojected and saved:", new_file))
}


