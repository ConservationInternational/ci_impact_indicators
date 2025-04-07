# Workflow for aggregating site data for SC and TF reports
# for instances where 2 variables are requested

# Contact: Anna Ballasiotes
# Last updated: 04/07/2025

library(tidyverse)
library(lubridate)
year <- "2024"


#####
##### Read in data if you have not yet
#####


# irrecoverable carbon broken up by ecosystem
sites_ic <- read_csv(
  paste0("results/FY", str_sub(year, start = 3), "_ImpactIndicators_IrrecoverableCarbon_Sites.csv")) %>%
  dplyr::select(!tonnes_ha_ic) # can't sum

overlaps_ic <- readRDS(
  paste0("results/FY", str_sub(year, start = 3), "_ImpactIndicators_IrrecoverableCarbon_Overlaps.rds")) %>%
  dplyr::select(!tonnes_ha_ic) # can't sum

#read in results tables as defined in script 2

sites_primary<- read_csv(
  paste0("results/FY", str_sub(year, start = 3), "_ImpactIndicators_Primary_Sites_ALL.csv"))

overlaps_primary <- readRDS(
  paste0("results/FY", str_sub(year, start = 3), "_ImpactIndicators_Primary_Overlaps.rds"))

### avoided emissions - by site, already prepped
ae <- readRDS("data/avoided_emissions/ae_by_site_by_year_2024.rds") %>% 
  dplyr::select(CI_ID, emissions_avoided_MgCO2e_2023)

#####
##### Correct any last-min edits
#####

overlaps_ic <- overlaps_ic %>%
  rename(Area_ha_list = Area_ha)  %>%
  filter(!map_chr(CI_ID, ~ .x[1]) %in% c("BNA1022", "BNA1007", "BNA1057", "BNA1073", "BNA1024", "BNA1052")) %>%
  filter(!map_chr(CI_ID, ~ .x[2]) %in% c("BNA1018", "PER1044")) %>%
  filter(!map_lgl(CI_ID, ~ any(.x %in% c("PHL1105", "PHL1106", "PHL1107", "PHL1108", "PHL1109", "PHL1110"))))

sites_ic <- sites_ic %>%
  filter(!CI_ID %in% c("PHL1105", "PHL1106", "PHL1107", "PHL1108", "PHL1109", "PHL1110"))

sites <- sites %>%
  filter(!CI_ID %in% c("PHL1105", "PHL1106", "PHL1107", "PHL1108", "PHL1109", "PHL1110"))

overlaps <- overlaps %>%
  filter(!map_chr(CI_ID, ~ .x[1]) %in% c("BNA1022", "BNA1007", "BNA1057", "BNA1073", "BNA1024", "BNA1052")) %>%
  filter(!map_chr(CI_ID, ~ .x[2]) %in% c("BNA1018", "PER1044")) %>%
  filter(!map_lgl(CI_ID, ~ any(.x %in% c("PHL1105", "PHL1106", "PHL1107", "PHL1108", "PHL1109", "PHL1110"))))


overlaps <- overlaps %>%
  rename(Area_ha_list = Area_ha) %>%
  rename(Area_ha = Area_ha_R)

sites_tf_new_biome <- sites_tf_new
sites_tf_contd_biome <- sites_tf_contd

overlaps_tf_new_biome <- overlaps_tf_new %>%
  filter(map_lgl(Biome, all_identical))

overlaps_tf_contd_biome <- overlaps_tf_contd %>%
  filter(map_lgl(Biome, all_identical))

groups <- c(
  "tf_new_biome",
  "tf_contd_biome")

vars_1 <- c(
  "Country")

vars_2 <- c(
  "Biome")


# Primary indicators -----

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
  
  for (i in seq_along(vars_1)) {
    
    user_group_1 <- vars_1[i]
    
    for (g in seq_along(vars_2)) {
      
      user_group_2 <- vars_2[g]
      
      user_groups <- c(user_group_1, user_group_2)
      
      sites_summary <- sites_df %>% 
        group_by_at(user_groups) %>% 
        primary_summarize() %>% 
        ungroup() 
      
      # summarize overlaps
      overlaps_summary_wo_ae <- overlaps_df %>% 
        dplyr::select((!!as.symbol(user_group_1)), (!!as.symbol(user_group_2)), 
                      Area_ha, rest_area, population, tstor_woody, 
                      tstor_soil, tstor_total, carbon_seq_potl) %>% 
        rowwise() %>% 
        # get a list of duplicated values
        mutate(duplicates = list((!!as.symbol(user_group_2))[duplicated((!!as.symbol(user_group_2)))])) %>% 
        ungroup() %>% 
        # remove those with no duplicates
        filter(!duplicates == "character(0)") %>% 
        # multiply values of interest by the number of times its duplicated & make negative 
        rowwise() %>%
        mutate({{ user_group_2 }} := list(unique(duplicates))) %>% 
        unnest({{ user_group_2 }}) %>% 
        ungroup() %>% 
        rowwise() %>% 
        mutate(duplicate_n = sum((!!as.symbol(user_group_2)) == duplicates, na.rm = TRUE)) %>% 
        mutate(across(
          .cols = c(Area_ha, rest_area, population, tstor_woody, 
                    tstor_soil, tstor_total, carbon_seq_potl),
          ~ .x * -duplicate_n))  %>% 
        
        # do again for second grouping if needed
        mutate(duplicates = list((!!as.symbol(user_group_1))[duplicated((!!as.symbol(user_group_1)))])) %>% 
        ungroup() %>% 
        filter(!duplicates == "character(0)") %>% 
        rowwise() %>%
        mutate({{ user_group_1 }} := unique(duplicates)) %>% 
        #      unnest({{ user_group_1 }}) %>% 
        ungroup() %>% 
        # summarize as you did with sites
        group_by_at(user_groups) %>% 
        other_summarize() %>% 
        ungroup() 
      
      #
      
      # bring in emissions avoided to overlap
      # break up ae by % area
      total_aes <- sites_summary %>% 
        dplyr::select((!!as.symbol(user_group_1)), (!!as.symbol(user_group_2)), 
                      Area_ha, emissions_avoided_MgCO2e) %>% 
        rename(total_area = Area_ha)
      
      
      ae_ref <- overlaps_summary_wo_ae %>% 
        dplyr::select((!!as.symbol(user_group_1)), (!!as.symbol(user_group_2)), Area_ha) %>% 
        left_join(total_aes, by = c(user_group_1, user_group_2)) %>% 
        mutate(pct_area = Area_ha/total_area) %>% 
        mutate(emissions_avoided_MgCO2e = emissions_avoided_MgCO2e * pct_area) %>% 
        dplyr::select((!!as.symbol(user_group_1)), (!!as.symbol(user_group_2)), emissions_avoided_MgCO2e)
      
      overlaps_summary <- overlaps_summary_wo_ae %>% 
        left_join(ae_ref, by = c(user_group_1, user_group_2))
      
      #
      
      # combine the two data frames and summarize to subtract the overlaps
      corrected_summary <- sites_summary %>% 
        bind_rows(overlaps_summary) %>% 
        group_by_at(user_groups) %>% 
        primary_summarize() %>% 
        ungroup() 
      
      write_csv(
        corrected_summary,
        file = paste0("results/summaries_", year, "/reports/FY24", "_", str_to_title(group_name), "_", str_to_title(user_group_1),
                      "_", str_to_title(user_group_2), ".csv")
      )
    }
  }
}


###############################
###############################
###############################
## Irrecoverable Carbon
###############################
###############################
###############################


groups <- c(
  "ic_climate_star")

vars_1 <- c(
  "Country")

vars_2 <- c(
  "Biome")

# define functions
ic_summarize <- function(df){
  df %>% 
    summarize(across(
      .cols = contains("_ic"),
      sum, na.rm = TRUE),
      .groups = "keep") 
}



for (i in seq_along(vars_1)) {
  
  user_group_1 <- vars_1[i]
  
  for (g in seq_along(vars_2)) {
    
    user_group_2 <- vars_2[g]
    
    user_groups <- c(user_group_1, user_group_2, "ecosystem")
    
    sites_ic_summary <- sites_ic_climate_star %>% 
      group_by_at(user_groups) %>% 
      ic_summarize() %>% 
      ungroup() 
    
    # summarize overlaps
    overlaps_ic_summary <- overlaps_ic_climate_star %>% 
      dplyr::select((!!as.symbol(user_group_1)), (!!as.symbol(user_group_2)), 
                    ecosystem,
                    contains("_ic")) %>% 
      rowwise() %>% 
      # get a list of duplicated values
      mutate(duplicates = list((!!as.symbol(user_group_2))[duplicated((!!as.symbol(user_group_2)))])) %>% 
      ungroup() %>% 
      # remove those with no duplicates
      filter(!duplicates == "character(0)") %>% 
      # multiply values of interest by the number of times its duplicated & make negative 
      rowwise() %>%
      mutate({{ user_group_2 }} := list(unique(duplicates))) %>% 
      unnest({{ user_group_2 }}) %>% 
      ungroup() %>% 
      rowwise() %>% 
      
      mutate(duplicate_n = sum((!!as.symbol(user_group_2)) == duplicates, na.rm = TRUE)) %>% 
      mutate(across(
        .cols = contains("_ic"),
        ~ .x * -duplicate_n))  %>% 
      
      # do again for second grouping if needed
      mutate(duplicates = list((!!as.symbol(user_group_1))[duplicated((!!as.symbol(user_group_1)))])) %>% 
      ungroup() %>% 
      filter(!duplicates == "character(0)") %>% 
      rowwise() %>%
      mutate({{ user_group_1 }} := unique(duplicates)) %>% 
      #      unnest({{ user_group_1 }}) %>% 
      ungroup() %>% 
      
      # summarize as you did with sites_ic
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
      file = paste0("results/summaries_", year, "/reports/FY24_Check_IC_", str_to_title(user_group_1),
                    "_", str_to_title(user_group_2), ".csv")
    )
  }
  
}



