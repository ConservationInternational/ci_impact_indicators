# Workflow for summarizing species data


library(tidyverse)
library("ritis")

year <- "2024"

# Use python code linked below to run / query species data 
# from GEE and create the "species_all_sites_CR-EN-VU" data
# that is required for this script
## https://colab.research.google.com/drive/16-e2GIzhPXAuFjvrLqFFg_S48RPcclxs?usp=sharing 
# You'll also need to upload an asset of the CI Sites to your GEE account

species_all_sites_CR_EN_VU <- read_csv("data/ecosystems/fy24species_all_sites_CR-EN-VU.csv")

sp_raw <- species_all_sites_CR_EN_VU 
species_merge <- species_all_sites_CR_EN_VU

sp_raw %>%
  dplyr::select(-site_id) %>%
  distinct() %>%
  dplyr::mutate(id = 1:n()) %>%
  relocate(id) -> species
saveRDS(species, file.path("data/ecosystems", 'species_fy24.rds'))


length(unique(sp_raw$sci_name))

species_sum_combos <- sp_raw %>%
  group_by(site_id, category, class) %>%
  summarize(combo_count = n(), .groups = 'drop')

sites <- read_csv(paste0(
  "results/FY", str_sub(year, start = 3), "_ImpactIndicators_Primary_Sites_ALL.csv"))


overlaps <- readRDS(
  paste0("results/FY", str_sub(year, start = 3), "_ImpactIndicators_Primary_Overlaps.rds"))

overlaps <- overlaps %>%
  rename(Area_ha_list = Area_ha) %>%
  rename(Area_ha = Area_ha_R)

relevant_site_columns <- sites %>%
  dplyr::select(CI_ID, Biome, Country, Area_Name, Division_Primary_CI, Geographic_Priority, CI_Portfolio)

species_merge <- sp_raw %>%
    rename(CI_ID = site_id)


species_merge_fields <- species_merge %>%
  left_join(relevant_site_columns, by = "CI_ID")
############# 

species_sum_sites <- species_merge_fields %>%
  group_by(CI_ID) %>%
  summarize(Unique_species_count = n_distinct(sci_name), .groups = 'drop')


species_sum_sites <- species_sum_sites %>%
  left_join(relevant_site_columns, by = "CI_ID")

# Species tables
species_sum_Country <- species_merge_fields %>%
  group_by(Country, category, class) %>%
  summarize(Species_count = n_distinct(sci_name), .groups = 'drop')

species_sum_Country_all <- species_merge_fields %>%
  group_by(Country) %>%
  summarize(Species_count = n_distinct(sci_name), .groups = 'drop')

species_sum_Division <- species_merge_fields %>%
  group_by(Division_Primary_CI, category, class) %>%
  summarize(Species_count = n_distinct(sci_name), .groups = 'drop')

species_sum_Portfolio <- species_merge_fields %>%
  dplyr::filter(CI_Portfolio != "<Null>") %>%
  group_by(CI_Portfolio, category, class) %>%
  summarize(Species_count = n_distinct(sci_name), .groups = 'drop')

species_sum_all_Division <- species_merge_fields %>% 
  dplyr::filter(Division_Primary_CI != "Not Applicable") %>%
  group_by(Division_Primary_CI) %>%
  summarize(species_count = n_distinct(sci_name))

species_sum_all_Portfolio <- species_merge_fields %>% 
  dplyr::filter(NPE != "<Null>") %>%
  summarize(species_count = n_distinct(sci_name))


species_sum_all <- species_merge_fields %>% 
  summarize(species_count = n_distinct(sci_name))

write_csv(species_sum_Country, 
          paste0("data/ecosystems/FY", str_sub(year, start = 3), 
                 "_species_sum_Country.csv"))

write_csv(species_sum_all_Division, 
          paste0("data/ecosystems/FY", str_sub(year, start = 3), 
                 "_species_sum_Division.csv"))

write_csv(species_sum_Country_all, 
          paste0("data/ecosystems/FY", str_sub(year, start = 3), 
                 "_species_sum_Country_all.csv"))

write_csv(species_sum_SLS, 
          paste0("data/ecosystems/FY", str_sub(year, start = 3), 
                 "_species_sum_SLS.csv"))

write_csv(species_sum_all, 
          paste0("data/ecosystems/FY", str_sub(year, start = 3), 
                 "_species_total.csv"))


write_csv(species_sum_sites, 
          paste0("data/ecosystems/FY", str_sub(year, start = 3), 
                 "_species_sum_sites.csv"))


write_csv(species_sum_Portfolio, 
          paste0("data/ecosystems/FY", str_sub(year, start = 3), 
                 "_species_sum_Portfolio.csv"))

species_sciNames <- read_csv("data/ecosystems/scientificnames.csv")

get_common_name <- function(scientific_name) {
  tryCatch({
    # Search ITIS for the scientific name
    results <- ritis::search_scientific(scientific_name)
    if (nrow(results) > 0) {
      # Extract the first common name if available
      common_name <- results$commonName[1]
      return(common_name)
    } else {
      return(NA)  # Return NA if no common name is found
    }
  }, error = function(e) {
    message(paste("Error with:", scientific_name, "->", e))
    return(NA)
  })
}


species_sciNames$Common_Name <- sapply(species_sciNames$scientificName, get_common_name)

# View the updated DataFrame
head(df)
