### Make maps of covariates

library(tidyverse)
library(terra)

data_folder <- 'avoided_emissions/data'

covariates_1 <- rast(
  list.files(file.path(data_folder), pattern = "covariate1", full.names = TRUE))

names(covariates_1) <- c('biome',
                         'elev',
                         'ecoregion',
                         'precip',
                         'slope',
                         'temp')

covariates_2 <- rast(
  list.files(file.path(data_folder), pattern = "covariate2", 
             full.names = TRUE))
names(covariates_2) <- c('dist_cities',
                         'crop_suitability',
                         'dist_roads',
                         'pa')

population <- rast(
  list.files("data/population", pattern = "_proj", full.names = TRUE))[[4]]
names(population) <- 'pop_2020'

biomass <- rast('avoided_emissions/data/covariate_biomass.tif')
names(biomass) <- c('total_biomass')

population_growth <- rast(
  'avoided_emissions/data/covariate_population_growth.tif')
names(population_growth) <- c('pop_growth')

covariates <- c(covariates_1, covariates_2, biomass, population, 
                population_growth)

covariates



for(i in 1:nlyr(covariates)){
  
  rast <- rast("avoided_emissions/data/covariate_cropland.tif")
  rast <- covariates[[i]] 
  title <- names(rast)
  
  png(paste0("avoided_emissions/plot_covariates/plots/", title, ".png"),
      width = 1500,
      height = 1000)
  
  print({
    plot(
      rast,
      axes = FALSE,
      legend = FALSE,
      maxcell = 1e6,
      main = str_to_title(title),
      cex.main = 3 # Title text size
      )
  })
  
  dev.off()
  
}
