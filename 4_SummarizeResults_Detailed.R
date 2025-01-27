# Workflow for summarizing impact indicators data

library(tidyverse)

year <- "2022"

# irrecoverable carbon broken up by ecosystem
sites_ic <- read_csv(
  paste0("results/FY", str_sub(year, start = 3), "_ImpactIndicators_IrrecoverableCarbon_Sites.csv")) %>%
  dplyr::select(!tonnes_ha_ic) # can't sum

overlaps_ic <- readRDS(
  paste0("results/FY", str_sub(year, start = 3), "_ImpactIndicators_IrrecoverableCarbon_Overlaps.rds")) %>%
  dplyr::select(!tonnes_ha_ic) # can't sum

# read in data for 'other indicators' for all sites and overlaps
sites <- read_csv(paste0(
  "results/FY", str_sub(year, start = 3), "_ImpactIndicators_Other_Sites.csv"))

overlaps <- readRDS(
  paste0("results/FY", str_sub(year, start = 3), "_ImpactIndicators_Other_Overlaps.rds"))

vars_1 <- c(
  "country",
  "ci_divisio",
  "ci_sls_2")

vars_2 <- c(
  "new_or_c_1",
  "inter_cat",
  "interventi",
  "restoratio")


# Irrecoverable carbon -----

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
    
    sites_ic_summary <- sites_ic %>% 
      group_by_at(user_groups) %>% 
      ic_summarize() %>% 
      ungroup() 
    
    # summarize overlaps
    overlaps_ic_summary <- overlaps_ic %>% 
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
      file = paste0("results/summaries_", year, "/reports/ImpactIndicators_IC_", str_to_title(user_group_1),
                    "_", str_to_title(user_group_2), ".csv")
    )
  }
  
}







# Other indicators -----

# Define summarizing function
# other_summarize <- function(df){
#   df %>% 
#     summarize(across(
#       .cols = c(area_ha, rest_area, population, 
#                 tstor_woody, tstor_soil, tstor_total, 
#                 carbon_seq_potl, emissions_avoided_mgco2e),
#       sum, na.rm = TRUE),
#       .groups = "keep") 
# }

other_summarize <- function(df){
  
  if(any(str_detect(colnames(df), 'emissions_avoided_mgco2e'))) {
    
    df %>% 
      summarize(across(
        .cols = c(area_ha, rest_area, population, 
                  tstor_woody, tstor_soil, tstor_total, 
                  carbon_seq_potl, emissions_avoided_mgco2e),
        sum, na.rm = TRUE),
        .groups = "keep") 
  } else {
    
    df %>% 
      summarize(across(
        .cols = c(area_ha, rest_area, population, 
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
    
    sites_summary <- sites %>% 
      group_by_at(user_groups) %>% 
      other_summarize() %>% 
      ungroup() 
    
    # summarize overlaps
    overlaps_summary_wo_ae <- overlaps %>% 
      dplyr::select((!!as.symbol(user_group_1)), (!!as.symbol(user_group_2)), 
                    area_ha, rest_area, population, tstor_woody, 
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
        .cols = c(area_ha, rest_area, population, tstor_woody, 
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
                    area_ha, emissions_avoided_mgco2e) %>% 
      rename(total_area = area_ha)
    
    
    ae_ref <- overlaps_summary_wo_ae %>% 
      dplyr::select((!!as.symbol(user_group_1)), (!!as.symbol(user_group_2)), area_ha) %>% 
      left_join(total_aes, by = c(user_group_1, user_group_2)) %>% 
      mutate(pct_area = area_ha/total_area) %>% 
      mutate(emissions_avoided_mgco2e = emissions_avoided_mgco2e * pct_area) %>% 
      dplyr::select((!!as.symbol(user_group_1)), (!!as.symbol(user_group_2)), emissions_avoided_mgco2e)
    
    overlaps_summary <- overlaps_summary_wo_ae %>% 
      left_join(ae_ref, by = c(user_group_1, user_group_2))
    
    #
    
    # combine the two data frames and summarize to subtract the overlaps
    corrected_summary <- sites_summary %>% 
      bind_rows(overlaps_summary) %>% 
      group_by_at(user_groups) %>% 
      other_summarize() %>% 
      ungroup() 
    
    write_csv(
      corrected_summary,
      file = paste0("results/summaries_", year, "/reports/ImpactIndicators", "_", str_to_title(user_group_1),
                    "_", str_to_title(user_group_2), ".csv")
    )
  }
  
}

