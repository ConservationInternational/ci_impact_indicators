# Data Acquisition
# Contact: Anna Ballasiotes
# Last updated: 04/7/2025

library(tidyverse)
library(raster)
library(terra)
library(sf)
library(janitor)
library("aws.s3")
library(dplyr)


# Acquire "data" folder
# Link to sharepoint below directs you to "ci_sites"
# and "previous" folders -- two folders that should not
# be uploaded to GitHub but that I think are useful to have.
# I also included a couple of other misc folders that may be 
# helpful at some point but also should not be uploaded.
# Only someone at CI should be able to access this link:
# https://conservation.sharepoint.com/:f:/s/CIImpactMetrics/EoaTN2DnDQZMnRM217VOjZcBCgKY45srIxkudS2Y1nYXIQ?e=wYsGee


privatekey <- #Check Private Key file
privatesecret <- #Check Private Secret File

ii_files <- get_bucket_df(
  bucket = "impact-indicators-data",
  region = "us-west-1",
  key = privatekey,
  secret = privatesecret,
) %>%
  as_tibble()

save_impact_data <- function(file) {

  # if object doesn't exist in bucket, return NA
  ok <- object_exists(
    object = file,
    bucket = "impact-indicators-data",
    region = "us-west-1",
    key = privatekey,
    secret = privatesecret,
  )
  if(!ok) return(NA_character_)
  if(endsWith(file, "/")) return("It's a folder!")

  # if object exists, save it and return file path
  save_object(
    object = file,
    bucket = "impact-indicators-data",
    region = "us-west-1",
    key = privatekey,
    secret = privatesecret,
    file = file
  )
}

#Map over function to get data

ii_files$Key %>%
  map_chr(save_impact_data)

#####
