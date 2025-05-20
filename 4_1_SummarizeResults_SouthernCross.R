# Workflow for creating 
# Total CI Footprint - Direct Contribution Report (TF)
# Southern Cross - Direct Contribution Report (SC)

# Contact: Anna Ballasiotes
# Last updated: 04/07/2025

# Much of this will likely need to change for SC to 2030
# I started a 'grouping' / 'query' table for Nathalie
# The bones / methodology should be here, but things will likely
# need to be re-grouped, filtered, changed, etc.
# There were a lot of changes / lot of last minute requests
# So I can't guarantee everything will work smoothly
# Would recommend working in chunks to catch any errors

library(tidyverse)
library(lubridate)
year <- "2024"


#####
##### Read in data
#####


# irrecoverable carbon broken up by ecosystem
sites_ic <- read_csv(
  paste0("results/FY", str_sub(year, start = 3), "_ImpactIndicators_IrrecoverableCarbon_Sites.csv")) %>%
  dplyr::select(!tonnes_ha_ic) # can't sum

overlaps_ic <- readRDS(
  paste0("results/FY", str_sub(year, start = 3), "_ImpactIndicators_IrrecoverableCarbon_Overlaps.rds")) %>%
  dplyr::select(!tonnes_ha_ic) # can't sum

#read in results tables as defined in script 2

sites_primary <- read_csv(
  paste0("results/FY", str_sub(year, start = 3), "_ImpactIndicators_Primary_Sites_ALL.csv"))

overlaps_primary <- readRDS(
  paste0("results/FY", str_sub(year, start = 3), "_ImpactIndicators_Primary_Overlaps.rds"))

### CAST AS DATE BC WHY

sites_primary$CI_End_Date <- parse_date_time(sites_primary$CI_End_Date, orders = c("dmy", "ymd", "mdy"))
sites_primary$CI_Start_Date <- parse_date_time(sites_primary$CI_Start_Date, orders = c("dmy", "ymd", "mdy"))
sites_primary$Gazettement_Date <- parse_date_time(sites_primary$Gazettement_Date, orders = c("dmy", "ymd", "mdy"))
sites_primary$Tree_Planting_Date <- parse_date_time(sites_primary$Tree_Planting_Date, orders = c("dmy", "ymd", "mdy"))

### avoided emissions - by site, already prepped
ae <- readRDS("data/avoided_emissions/ae_by_site_by_year_2024.rds") %>% 
  dplyr::select(CI_ID, emissions_avoided_MgCO2e_2023)

#####
##### Correct small items
#####

#Remove these CI'IDs from the main dataset
#CI_ID IN ('BNA1001', 'BNA1002', 'BNA1003', 'BNA1017', 'BNA1018', 'BNA1019', 'BNA1074', 'BNA1075', 'BNA1076', 'BNA1077', 'BNA1078', 'BNA1079', 'BNA1080', 'BNA1081')

overlaps_ic <- overlaps_ic %>%
  rename(Area_ha_list = Area_ha)  %>%
  filter(!map_chr(CI_ID, ~ .x[1]) %in% c("BNA1022", "BNA1007", "BNA1057", "BNA1073", "BNA1024", "BNA1052")) %>%
  filter(!map_chr(CI_ID, ~ .x[2]) %in% c("BNA1018", "PER1044")) %>%
  filter(!map_lgl(CI_ID, ~ any(.x %in% c("PHL1105", "PHL1106", "PHL1107", "PHL1108", "PHL1109", "PHL1110"))))

sites_ic <- sites_ic %>%
  filter(!CI_ID %in% c("PHL1105", "PHL1106", "PHL1107", "PHL1108", "PHL1109", "PHL1110"))

sites_primary <- sites_primary %>%
  filter(!CI_ID %in% c("PHL1105", "PHL1106", "PHL1107", "PHL1108", "PHL1109", "PHL1110"))

overlaps_primary <- overlaps_primary %>%
  filter(!map_chr(CI_ID, ~ .x[1]) %in% c("BNA1022", "BNA1007", "BNA1057", "BNA1073", "BNA1024", "BNA1052")) %>%
  filter(!map_chr(CI_ID, ~ .x[2]) %in% c("BNA1018", "PER1044")) %>%
  filter(!map_lgl(CI_ID, ~ any(.x %in% c("PHL1105", "PHL1106", "PHL1107", "PHL1108", "PHL1109", "PHL1110"))))

overlaps_primary <- overlaps_primary %>%
  rename(Area_ha_list = Area_ha) %>%
  rename(Area_ha = Area_ha_R)

#################################
#################################
#################################
## For Total Footprint Report
#################################
#################################
#################################

## ALL TOTAL FOOTPRINT
sites_tf <- sites_primary %>%
  filter(Site_Status %in% c("New this year", "Continued this year"))
overlaps_tf <- overlaps_primary %>%
  filter(map_lgl(CI_ID, ~all(.x %in% sites_tf$CI_ID)))


sites_tf_division <- sites_primary %>%
  filter(Site_Status %in% c("New this year", "Continued this year"))

overlaps_tf_division <- overlaps_tf %>%
                            filter(!map_lgl(CI_ID, ~ any(.x %in% c("BNA1071"))))
                                      
sites_ic_tf <- sites_ic %>%
  filter(Site_Status %in% c("New this year", "Continued this year"))
overlaps_ic_tf <- overlaps_ic %>%
  filter(map_chr(Site_Status, ~ .x[1]) %in% c("New this year", "Continued this year"))

## NEW
sites_tf_new <- sites_tf %>%
  filter(Site_Status %in% c("New this year"))
overlaps_tf_new <- overlaps_tf %>%
  filter(map_lgl(CI_ID, ~all(.x %in% sites_tf_new$CI_ID)))

## CONTD
sites_tf_contd <- sites_tf %>%
  filter(Site_Status %in% c("Continued this year"))
overlaps_tf_contd <- overlaps_tf %>%
  filter(map_lgl(CI_ID, ~all(.x %in% sites_tf_contd$CI_ID)))

## ALL RESTORATION
sites_restor <- sites_primary %>%
  filter(Intervention_Category == "Restoration Areas")
overlaps_restor <- overlaps_primary %>%
  filter(map_lgl(CI_ID, ~all(.x %in% sites_restor$CI_ID)))

## NEW RESTORATION
sites_restor_new <- sites_restor %>%
  filter(Site_Status == "New this year")
overlaps_restor_new <- overlaps_restor %>%
  filter(map_lgl(CI_ID, ~all(.x %in% sites_restor_new$CI_ID)))




all_identical <- function(x) {
  length(unique(x)) == 1
}

overlaps_tf_contd_division <- overlaps_tf_contd %>%
  filter(map_lgl(Division_Primary_CI, all_identical))

sites_tf_contd_division <- sites_tf_contd


########################################################
###############################
###############################
## For Southern Cross Reporting
###############################
###############################
###############################


## ALL OCEAN

sites_ocean <- sites_primary %>%
  filter(Site_Status %in% c("New this year", "Continued this year")) %>%
  filter(!CI_Start_Date < "2019-01-01") %>%
  filter(!grepl("Proposed", Intervention_Type)) %>%
  filter(Biome %in% c("Marine","Terrestrial + Marine") |
           Division_Primary_CI %in% c('Blue Nature Alliance (C4O - BNA)', 'Center for Oceans (C4O)') |
           Division_Secondary_CI %in% c('Blue Nature Alliance (C4O - BNA)', 'Center for Oceans (C4O)'))

overlaps_ocean <- overlaps_primary %>%
  filter(map_lgl(CI_ID, ~all(.x %in% sites_ocean$CI_ID)))

write_csv(
  sites_ocean,
  file = paste0("results/summaries_", year, "/reports/FY24_sites_ocean.csv")
)

sites_ocean_new <- sites_ocean %>%
  filter(Intervention_Category == "Contractual Conservation" & Site_Status == "New this year" |
           Intervention_Type == "Ocean Conservation Area (National or Regional)" & Site_Status == "New this year")
overlaps_ocean_new <- overlaps_primary %>%
  filter(map_lgl(CI_ID, ~all(.x %in% sites_ocean_new$CI_ID)))



sites_ocean_strengthened <- sites_ocean %>%
  filter(Intervention_Category %in% c("Protected Areas", "Contractual Conservation", "Areas Governed by IPLCs"))
overlaps_ocean_strengthened <- overlaps_primary %>%
  filter(map_lgl(CI_ID, ~all(.x %in% sites_ocean_strengthened$CI_ID)))


## Climate Star

sites_climate <- sites_primary %>%
  filter(Site_Status %in% c("New this year", "Continued this year")) %>%
  filter(!CI_Start_Date < "2019-01-01") %>%
  filter(!grepl("Proposed", Intervention_Type)) %>%
  filter(Biome %in% c("Terrestrial", "Terrestrial + Marine",
                      "Terrestrial + Freshwater"))
overlaps_climate <- overlaps_primary %>%
  filter(map_lgl(CI_ID, ~all(.x %in% sites_climate$CI_ID)))

write_csv(
  sites_climate,
  file = paste0("results/summaries_", year, "/reports/FY24_sites_climate.csv")
)

sites_climate_new <- sites_climate %>%
  filter(Intervention_Category == "Protected Areas" & Gazettement_Date > "2023-07-01" |
           Intervention_Category == "Contractual Conservation Areas" & Site_Status == "New this year")
overlaps_climate_new <- overlaps_primary %>%
  filter(map_lgl(CI_ID, ~all(.x %in% sites_climate_new$CI_ID)))


sites_climate_strengthened <- sites_climate %>%
  filter(grepl("Yes", Improved_Management)) %>%
  filter(Intervention_Category == "Protected Areas" & Gazettement_Date < "2023-07-01" | 
           Intervention_Category == "Contractual Conservation Areas" & Site_Status == "Continued this year" |
          Intervention_Category == "Areas Governed by IPLCs")
overlaps_climate_strengthened <- overlaps_primary %>%
  filter(map_lgl(CI_ID, ~all(.x %in% sites_climate_strengthened$CI_ID)))


sites_climate_star <- sites_climate_new %>%
  rbind(sites_climate_strengthened)
overlaps_climate_star <- overlaps_primary %>%
  filter(map_lgl(CI_ID, ~all(.x %in% sites_climate_star$CI_ID)))


## SEQUESTRATION

sites_climate_restor <- sites_restor %>%
  filter(!CI_Start_Date < "2019-01-01") %>%
  filter(!Site_Status_Sequestration %in% c("Unknown")) %>%
  filter(Intervention_Type %in% c("Restoration - Natural Regeneration",
                                  "Rangeland Restoration - Planned Grazing",
                                  "Restoration - Agroforestry",
                                  "Mangrove Shrub Restoration",
                                  "Restoration - Mangrove Tree",
                                  "Restoration - Plantations & Woodlots (Eucalyptus)",
                                  "Restoration - Plantations & Woodlots (Mixed 50/50)",
                                  "Restoration - Enrichment Planting/Assisted Natural Regeneration",
                                  "Restoration - Seed Dispersal",
                                  "Restoration - Silvopasture",
                                  "Restoration - Plantations & Woodlots  (Other Broadleaf)" ,
                                  "Restoration - Plantations & Woodlots (Pine)",
                                  "Restoration - Plantations & Woodlots (Other Conifer)"))

overlaps_climate_restor <- overlaps_primary %>%
  filter(map_lgl(CI_ID, ~all(.x %in% sites_climate_restor$CI_ID)))


sites_climate_restor_new <- sites_climate_restor %>%
  filter(Site_Status == "New this year")
overlaps_climate_restor_new <- overlaps_primary %>%
  filter(map_lgl(CI_ID, ~all(.x %in% sites_climate_restor_new$CI_ID)))


sites_southern_cross <- sites_climate %>%
  rbind(sites_ocean_strengthened)
write_csv(
  sites_southern_cross,
  file = paste0("results/summaries_", year, "/reports/FY24_sites_southern_cross.csv")
)

duplicates_sc <- sites_southern_cross$CI_ID[duplicated((sites_southern_cross$CI_ID))]
oceans_coast <- sites_ocean_strengthened %>% filter(Biome == "Terrestrial + Marine")


### Total Footprint

groups <- c(
  "tf",
  "tf_new",
  "tf_contd",
  "restor", 
  "restor_new",
  "seq")

vars_1 <- c(
  "Country",
  "CI_Portfolio")

groups <- c(
  "climate_restor",
  "climate_restor_new")

groups <- c(
  "tf_division",
  "tf_contd_division")

groups <- c(
  "tf_new",
  "restor", 
  "restor_new",
  "seq")

vars_1 <- c(
  "Division_Primary_CI")


## Southern Cross
groups <- c(
  "climate_new",
  "climate_strengthened",
  "climate_star",
  "climate_restor",
  "climate_restor_new",
  "ocean_new",
  "ocean_strengthened")

groups <- c(
  "climate_restor_new")

groups <- c(
  "ocean_new",
  "ocean_strengthened")

vars_1 <- c(
  "Country",
  "CI_Portfolio",
  "Division_Primary_CI")

vars_1 <- c(
  "Division_Primary_CI")

###############################
###############################
###############################
## Primary Variables
###############################
###############################
###############################



primary_summarize <- function(df){
  
  if(any(str_detect(colnames(df), 'emissions_avoided_MgCO2e'))) {
    
    df %>% 
      summarize(across(
        .cols = c(Area_ha, rest_area, population, 
                  tstor_woody, tstor_soil, tstor_total, 
                  carbon_seq_potl, emissions_avoided_MgCO2e),
        sum, na.rm = TRUE),
        .groups = "keep") 
  } else {
    
    df %>% 
      summarize(across(
        .cols = c(Area_ha, rest_area, population, 
                  tstor_woody, tstor_soil, tstor_total, 
                  carbon_seq_potl),
        sum, na.rm = TRUE),
        .groups = "keep") 
  }
}

# Loop through vars

for (n in seq_along(groups)) {
  group_name <- paste0("sites_", groups[n])
  overlaps_name <- paste0("overlaps_", groups[n])
  
  sites_df <- get(group_name)
  overlaps_df <- get(overlaps_name)

  for (i in seq_along(vars_1)) {
    
    user_group_1 <- vars_1[i]
    
    user_groups <- c(user_group_1)
    
    sites_summary <- sites_df %>% 
      group_by_at(user_groups) %>% 
      primary_summarize() %>% 
      ungroup() 
    
    # summarize overlaps
    overlaps_summary_wo_ae <- overlaps_df %>% 
      dplyr::select((!!as.symbol(user_group_1)), 
                    Area_ha, rest_area, population, tstor_woody, 
                    tstor_soil, tstor_total, carbon_seq_potl, n.overlaps) %>% 
      rowwise() %>% 
      # do again for second grouping if needed
      mutate(duplicates = list((!!as.symbol(user_group_1))[duplicated((!!as.symbol(user_group_1)))])) %>% 
      ungroup() %>% 
      filter(!duplicates == "character(0)") %>% 
      rowwise() %>%
      mutate({{ user_group_1 }} := unique(duplicates)) %>% 
      unnest({{ user_group_1 }}) %>% 
      ungroup() %>% 
      rowwise() %>% 
      mutate(duplicate_n = sum((!!as.symbol(user_group_1)) == duplicates, na.rm = TRUE)) %>% 
      mutate(across(
        .cols = c(Area_ha, rest_area, population, tstor_woody, 
                  tstor_soil, tstor_total, carbon_seq_potl),
        ~ .x * -duplicate_n))  %>% 
      
      # summarize as you did with sites
      group_by_at(user_groups) %>% 
      primary_summarize() %>% 
      ungroup() 
    
    ##
    # overlaps_df_corrected <- overlaps_df %>% 
    #   mutate(across(
    #     .cols = c(Area_ha, rest_area, population, 
    #               tstor_woody, tstor_soil, tstor_total, 
    #               carbon_seq_potl),
    #     ~ as.numeric(.x) * - (n.overlaps - 1))) %>% 
    #   dplyr::select(!n.overlaps)
    
    # bring in emissions avoided to overlap
    # break up ae by % area
    
    if (nrow(overlaps_summary_wo_ae)== 0) {
      
      write_csv(
        sites_summary,
        file = paste0("results/summaries_", year, "/reports/FY24", "_", str_to_title(group_name), "_", str_to_title(user_group_1), ".csv")
      ) 
      }
    else {
      total_aes <- sites_summary %>% 
        dplyr::select((!!as.symbol(user_group_1)), 
                      Area_ha, emissions_avoided_MgCO2e) %>% 
        rename(total_area = Area_ha)
      
      
      ae_ref <- overlaps_summary_wo_ae %>% 
        dplyr::select((!!as.symbol(user_group_1)), Area_ha) %>% 
        left_join(total_aes, by = c(user_group_1)) %>% 
        mutate(pct_area = Area_ha/total_area) %>% 
        mutate(emissions_avoided_MgCO2e = emissions_avoided_MgCO2e * pct_area) %>% 
        dplyr::select((!!as.symbol(user_group_1)), emissions_avoided_MgCO2e)
      
      overlaps_summary <- overlaps_summary_wo_ae %>% 
        left_join(ae_ref, by = c(user_group_1))
  
      # combine the two data frames and summarize to subtract the overlaps
      corrected_summary <- sites_summary %>% 
        bind_rows(overlaps_summary) %>% 
        group_by_at(user_groups) %>% 
        primary_summarize() %>% 
        ungroup() 
      write_csv(
        corrected_summary,
        file = paste0("results/summaries_", year, "/reports/FY24", "_", str_to_title(group_name), "_", str_to_title(user_group_1), ".csv")
      ) 
    }
  }
}



###############################
###############################
###############################
## Irrecoverable Carbon ######
###############################
###############################
###############################

sites_ic_climate_strengthened <- sites_ic %>%
  filter(CI_ID %in% sites_climate_strengthened$CI_ID)
overlaps_ic_climate_strengthened <- overlaps_ic_corrected %>%
  filter(map_lgl(CI_ID, ~all(.x %in% sites_ic_climate_strengthened$CI_ID)))

sites_ic_climate_new <- sites_ic %>%
  filter(CI_ID %in% sites_climate_new$CI_ID)
overlaps_ic_climate_new <- overlaps_ic_corrected %>%
  filter(map_lgl(CI_ID, ~all(.x %in% sites_ic_climate_new$CI_ID)))


sites_ic_climate <- sites_ic %>%
  filter(CI_ID %in% sites_climate$CI_ID)
overlaps_ic_climate <- overlaps_ic_corrected %>%
  filter(map_lgl(CI_ID, ~all(.x %in% sites_ic_climate$CI_ID)))

sites_ic_climate_star <- sites_ic %>%
  filter(CI_ID %in% sites_climate_star$CI_ID)
overlaps_ic_climate_star <- overlaps_ic_corrected %>%
  filter(map_lgl(CI_ID, ~all(.x %in% sites_ic_climate_star$CI_ID)))

groups_ic <- c(
  "climate_star",
  "climate_new",
  "climate_strengthened")


vars_1_ic <- c(
  "Country",
  "Division_Primary_CI",
  "CI_Portfolio")


# define functions
ic_summarize <- function(df){
  df %>% 
    summarize(across(
      .cols = contains("_ic"),
      sum, na.rm = TRUE),
      .groups = "keep") 
}


for (n in seq_along(groups_ic)) {
  group_name <- paste0("sites_ic_", groups_ic[n])
  overlaps_name <- paste0("overlaps_ic_", groups_ic[n])
  
  sites_ic_df <- get(group_name)
  overlaps_ic_df <- get(overlaps_name)
  
  for (i in seq_along(vars_1_ic)) {
    
    user_group_1 <- vars_1_ic[i]
    
    user_groups <- c(user_group_1, "ecosystem")
    
    sites_ic_summary <- sites_ic_df %>% 
      group_by_at(user_groups) %>% 
      ic_summarize() %>% 
      ungroup() 
    
    # summarize overlaps
    overlaps_ic_summary <- overlaps_ic_df %>% 
      dplyr::select((!!as.symbol(user_group_1)), 
                    ecosystem,
                    contains("_ic")) %>% 
      rowwise() %>%
      # do again for second grouping if needed
      mutate(duplicates = list((!!as.symbol(user_group_1))[duplicated((!!as.symbol(user_group_1)))])) %>% 
      ungroup() %>% 
      filter(!duplicates == "character(0)") %>% 
      rowwise() %>% 
      mutate({{ user_group_1 }} := unique(duplicates)) %>% 
      #unnest({{ user_group_1 }}) %>% 
      ungroup() %>% 
      mutate(duplicate_n = sum((!!as.symbol(user_group_1)) == duplicates, na.rm = TRUE)) %>% 
      
      # mutate(across(
      #   .cols = contains("_ic"),
      #   ~ .x * -duplicate_n))  %>%
      
      #summarize as you did with sites_ic
      group_by_at(user_groups) %>% 
      ic_summarize() %>% 
      ungroup() 
    
    # combine the two data frames and summarize to subtract the overlaps_ic
    corrected_summary <- sites_ic_summary %>% 
      bind_rows(overlaps_ic_summary) %>% 
      group_by_at(user_groups) %>% 
      ic_summarize() %>% 
      ungroup() 
    
    write_csv(
      corrected_summary,
      file = paste0("results/summaries_", year, "/reports/FY24_IC_", str_to_title(group_name), "_", str_to_title(user_group_1),".csv")
    )
  }
}



## COUNTRY + AVOIDED EMISSIONS BY COUNTRY

primary_overlap_summarize <- function(df){
  
  if(any(str_detect(colnames(df), 'emissions_avoided_MgCO2e'))) {
    
    df %>% 
      summarize(across(
        .cols = c(Area_ha, rest_area, population, 
                  tstor_woody, tstor_soil, tstor_total, 
                  carbon_seq_potl, emissions_avoided_MgCO2e),
        sum, na.rm = TRUE),
        .groups = "keep") 
  } else {
    
    df %>% 
      summarize(across(
        .cols = c(Area_ha, rest_area, population, 
                  tstor_woody, tstor_soil, tstor_total, 
                  carbon_seq_potl),
        sum, na.rm = TRUE),
        .groups = "keep") 
  }
}

# primary indicators

country_sites_primary_avoid <- sites_primary %>% 
  group_by(Country, AvoidedOrAdded) %>% 
  primary_summarize() %>% 
  ungroup() 

country_overlaps_primary_wo_ae_avoid <- overlaps_primary_corrected %>%
  rowwise() %>% 
  mutate(Country = (unique(Country))) %>% # all countries within polygons are the same
  ungroup() %>% 
  group_by(Country) %>% 
  primary_overlap_summarize() %>% 
  ungroup() 

# bring in emissions avoided to overlap
# break up ae by % area
total_aes <- country_sites_primary_avoid %>% 
  dplyr::select(Country, Area_ha, emissions_avoided_MgCO2e, AvoidedOrAdded) %>% 
  rename(total_area = Area_ha)

ae_ref <- country_overlaps_primary_wo_ae_avoid %>% 
  dplyr::select(Country, Area_ha) %>% 
  left_join(total_aes, by = "Country") %>% 
  mutate(pct_area = Area_ha/total_area) %>% 
  mutate(emissions_avoided_MgCO2e = emissions_avoided_MgCO2e * pct_area) %>% 
  dplyr::select(Country, emissions_avoided_MgCO2e)

ae_ref$AvoidedOrAdded <- ifelse(ae_ref$emissions_avoided_MgCO2e >= 0, "Added", "Avoided")


country_overlaps_primary_avoid <- country_overlaps_primary_wo_ae_avoid %>% 
  left_join(ae_ref, by = 'Country')

#Each row in `x` is expected to match at most 1 row in `y`.
# Row 1 of `x` matches multiple rows.
# If multiple matches are expected, set `multiple = "all"` to silence this warnin

country_corrected_primary_avoid <- country_sites_primary_avoid %>% 
  bind_rows(country_overlaps_primary_avoid) %>% 
  group_by(Country, AvoidedOrAdded) %>% 
  primary_summarize() %>%  # can sum bc overlap values are negative
  ungroup() 

country_corrected_primary_bycountry <- country_sites_primary_avoid %>% 
  bind_rows(country_overlaps_primary_avoid) %>% 
  group_by(Country) %>% 
  primary_summarize() %>%  # can sum bc overlap values are negative
  ungroup() 

write_csv(country_corrected_primary_avoid, 
          paste0("results/summaries_", year, "/reports/FY24_CountryEmissions_PrimaryIndicators_AvoidAdd.csv"))

write_csv(country_corrected_primary_bycountry, 
          paste0("results/summaries_", year, "/ImpactIndicators_CountryEmissions_PrimaryIndicators.csv"))



# climate star indicators
### MAKE SURE TO REMOVE THE EXTRACTIVES & REDD+

country_sites_climate_star_avoid <- sites_climate_star %>% 
  filter(Biome %in% c("Terrestrial",
                      "Terrestrial + Freshwater")) %>%
  group_by(Country, AvoidedOrAdded) %>% 
  primary_summarize() %>% 
  ungroup() 

## Correct to make negative
overlaps_climate_star_corrected <- overlaps_climate_star %>% 
  mutate(across(
    .cols = c(Area_ha, rest_area, population, 
              tstor_woody, tstor_soil, tstor_total, 
              carbon_seq_potl),
    ~ as.numeric(.x) * - (n.overlaps - 1))) %>% 
  dplyr::select(!n.overlaps)

country_overlaps_climate_star_wo_ae_avoid <- overlaps_climate_star_corrected %>%
  rowwise() %>% 
  mutate(Country = (unique(Country))) %>% # all countries within polygons are the same
  ungroup() %>% 
  group_by(Country) %>% 
  primary_overlap_summarize() %>% 
  ungroup() 

# bring in emissions avoided to overlap
# break up ae by % area
total_aes <- country_sites_climate_star_avoid %>% 
  dplyr::select(Country, Area_ha, emissions_avoided_MgCO2e, AvoidedOrAdded) %>% 
  rename(total_area = Area_ha)

ae_ref <- country_overlaps_climate_star_wo_ae_avoid %>% 
  dplyr::select(Country, Area_ha) %>% 
  left_join(total_aes, by = "Country") %>% 
  mutate(pct_area = Area_ha/total_area) %>% 
  mutate(emissions_avoided_MgCO2e = emissions_avoided_MgCO2e * pct_area) %>% 
  dplyr::select(Country, emissions_avoided_MgCO2e)

ae_ref$AvoidedOrAdded <- ifelse(ae_ref$emissions_avoided_MgCO2e >= 0, "Added", "Avoided")


country_overlaps_climate_star_avoid <- country_overlaps_climate_star_wo_ae_avoid %>% 
  left_join(ae_ref, by = 'Country')

#Each row in `x` is expected to match at most 1 row in `y`.
# Row 1 of `x` matches multiple rows.
# If multiple matches are expected, set `multiple = "all"` to silence this warnin

country_corrected_climate_star_avoid <- country_sites_climate_star_avoid %>% 
  bind_rows(country_overlaps_climate_star_avoid) %>% 
  group_by(Country, AvoidedOrAdded) %>% 
  primary_summarize() %>%  # can sum bc overlap values are negative
  ungroup() 

country_corrected_primary_bycountry <- country_sites_climate_star_avoid %>% 
  bind_rows(country_overlaps_climate_star_avoid) %>% 
  group_by(Country) %>% 
  primary_summarize() %>%  # can sum bc overlap values are negative
  ungroup() 

write_csv(country_corrected_climate_star_avoid, 
          paste0("results/summaries_", year, "/reports/FY24_CountryEmissions_ClimateStar_AvoidAdd_TerrestrialOnly.csv"))


################ 
ic_summary <- sites_ic %>%
  group_by(CI_ID) %>%
  summarise(across(ends_with("_ic"), 
                   list(sum = sum), 
                   na.rm = TRUE, 
                   .names = "{.col}_{.fn}"))

write_csv(ic_summary, 
          paste0("results/summaries_", year, "/reports/FY24_IC_summarized_by_site.csv"))
