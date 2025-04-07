library(raster)
library(fasterize)
library(sf)
library(tidyverse)
library(units)
library(foreach)
library(mapview)
library(exactextractr)
library(raster)
library(tictoc)

options("optmatch_max_problem_size"=Inf)

data_folder <- "data/avoided_emissions/"

# Basic function to extract variable names from a formula object
get_names <- function(f) {
    f <- paste0(as.character(f), collapse=' ')
    v <- strsplit(f, split='[+ ~]')[[1]]
    v <- v[v != '']
    gsub('strata\\(([a-zA-Z_]*)\\)', '\\1', v)
}

###############################################################################
### Final data setup

f <- treatment ~ lc_2015_agriculture + precip + temp + elev + slope + 
    dist_cities + dist_roads + crop_suitability + pop_2020 + pop_growth + 
    total_biomass

#readRDS('ae_output/formula.RDS')
saveRDS(f, 'ae_output/formula.RDS')


covariates <- brick(paste0(data_folder, 'covariates_covariates_2024.tif'))
names(covariates) <- read_csv(paste0(data_folder, 'covariates_covariates.csv'))$names
lc_2015 <- brick(paste0(data_folder, 'covariates_lc_2015.tif'))
names(lc_2015) <- read_csv(paste0(data_folder, 'covariates_lc_2015.csv'))$names

## adb note: probably stack covariates from the 2nd script rather than loading tif - it doesn't
## store very well
fc <- brick(paste0(data_folder, 'covariates_fc.tif'))
names(fc) <- read_csv(paste0(data_folder, 'covariates_fc.csv'))$names
fcc <- brick(paste0(data_folder, 'covariates_fc_change.tif'))
names(fcc) <- read_csv(paste0(data_folder, 'covariates_fc_change.csv'))$names

d <- stack(covariates, lc_2015, fc, fcc)

# Ensure only layers in the formula are included (so extra data isn't being 
# passed around)
d <- d[[c(get_names(f),
          'region',
          'ecoregion',
          'pa',
           paste0('fc200', seq(0, 9)),
           paste0('fc20', seq(10, 23)),
           paste0('fcc_200', seq(1, 9)),
           paste0('fcc_20', seq(10, 23))
          )]]
write_csv(data.frame(names=names(d)), file='data/avoided_emissions/all_covariates_names.csv')

###############################################################################
###  Load sites and covariates

sites <- readRDS(paste0(data_folder, 'sites_FY2024.RDS'))
dim(sites)

## ADB this is a variable we could maybe change... leaves us with fewer sites
# Filter to only sites over 100 ha
sites <- sites[!sites$area_cea < as_units(100, 'hectares'), ]
dim(sites)

# Select only sites that are not Rangeland Restoration or
# 'Rangeland Restoration - Alien/Native Bush Thinning
# because these we don't want to calculate AE for
# Would consider filtering other AE per table of queries
# and planned updates for FY25. Consider what sites we want
# to be performing AE on. 

sites$Rangeland <- FALSE
sites$Rangeland[sites$Intervention_Type == 'Rangeland Restoration - Planned Grazing' | sites$Intervention_Type == 'Rangeland Restoration - Alien/Native Bush Thinning'] <- TRUE
table(sites$Rangeland)
sites <- sites[!sites$Rangeland, ]
dim(sites)

regions <- readRDS(paste0(data_folder, 'regions.RDS'))
regions_rast <- fasterize(regions, raster(d[[1]]), field='level1_ID')
names(regions_rast) <- 'region'

## SAVE IT
writeRaster(regions_rast, filename='data/avoided_emissions/regions.tif', 
            overwrite=TRUE, options="COMPRESS=LZW", datatype="INT2S")

d <- stack(d, regions_rast)

###############################################################################
###  Load sites and covariates

# Run extractions of treatment points individually to make catching any polygon 
# errors easier
treatment_key <- foreach(n=1:nrow(sites), .combine=rbind) %do% {
    print(n)
    exact_extract(d$region, sites[n, ], include_cell=TRUE, 
                  include_cols=c('CI_ID', 'CI_Start_Year', 'CI_End_Year'))[[1]]
}

## This renames the value column to region drops any NA regions, and saves it
treatment_key %>%
    rename(region=value) -> treatment_key
treatment_key <- treatment_key[!is.na(treatment_key$region), ]
saveRDS(treatment_key, 'ae_output/treatment_cell_key.RDS')

# Run extraction of control and treatment data by region to make the problem 
# tractable in-memory

out <- foreach(this_region_ID=unique(treatment_key$region), 
        .packages=c('exactextractr', 'sf')) %do% {
    if (file.exists(paste0('ae_output/treatments_and_controls_', this_region_ID, 
                           '.RDS'))) {
        print(paste0('Skipping ', this_region_ID, '. Already processed.'))
    } else {
        this_region <- regions[regions$level1_ID == this_region_ID, ]
        vals <- exact_extract(d, this_region, include_cell=TRUE)[[1]]
        saveRDS(vals, paste0('ae_output/treatments_and_controls_', 
                             this_region_ID, '.RDS'))
    }
        }

