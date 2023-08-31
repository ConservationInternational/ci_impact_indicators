library(sf)
library(tidyverse)
library(raster)
library(terra)
library(fasterize)
library(gdalUtils)
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

# Load covariates
covariates_1 <- raster::stack(
    list.files(file.path(data_folder), pattern = "covariate1", full.names = TRUE))

names(covariates_1) <- c('biome',
                         'elev',
                         'ecoregion',
                         'precip',
                         'slope',
                         'temp')
NAvalue(covariates_1) <- -32768

covariates_2 <- raster::stack(
    list.files(file.path(data_folder), pattern = "covariate2", 
               full.names = TRUE))
names(covariates_2) <- c('dist_cities',
                         'crop_suitability',
                         'dist_roads',
                         'pa')
NAvalue(covariates_2) <- -32768

population <- raster::stack(
    list.files("data/population", pattern = "_proj", full.names = TRUE))
names(population) <-
    c('pop_2000', 'pop_2005', 'pop_2010', 'pop_2015', 'pop_2020')
NAvalue(population) <- -32768

biomass <- raster('data/avoided_emissions/covariate_biomass.tif')
names(biomass) <- c('total_biomass')

population_growth <- raster(
    'data/avoided_emissions/covariate_population_growth.tif')
names(population_growth) <- c('pop_growth')

covariates <- stack(covariates_1, covariates_2, biomass, population, 
                    population_growth)
writeRaster(
    covariates, 
    filename='data/avoided_emissions/covariates_covariates.tif', 
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



# forest cover
fc <- raster::stack(list.files(data_folder, pattern = "forest_cover_20",
                               full.names = TRUE))
NAvalue(fc) <- -32768
writeRaster(fc, filename='data/avoided_emissions/covariates_fc.tif', 
            overwrite=TRUE, options="COMPRESS=LZW", datatype="INT2S")
write_csv(data.frame(names=names(fc)), 'data/avoided_emissions/covariates_fc.csv')


# forest cover change
fcc <- raster::stack(list.files(data_folder, pattern = "forest_cover_change_20",
                               full.names = TRUE))
NAvalue(fcc) <- -32768
writeRaster(fcc, filename='data/avoided_emissions/covariates_fc_change.tif', 
            overwrite=TRUE, options="COMPRESS=LZW", datatype="INT2S")
write_csv(data.frame(names=names(fcc)), 'data/avoided_emissions/covariates_fc_change.csv')

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