# Workflow for summarizing impact indicators data

library(tidyverse)

# irrecoverable carbon broken up by ecosystem
sites_ic <- read_csv("results/FY21_ImpactIndicators_IrrecoverableCarbon_Sites.csv") %>%
  dplyr::select(!tonnes_ha_ic) # can't sum
overlaps_ic <- readRDS("results/FY21_ImpactIndicators_IrrecoverableCarbon_Overlaps.rds") %>%
  dplyr::select(!tonnes_ha_ic) # can't sum

# read in data for 'other indicators' for all sites and overlaps
sites <- read_csv("results/FY21_ImpactIndicators_Other_Sites.csv")
overlaps <- readRDS("results/FY21_ImpactIndicators_Other_Overlaps.rds")

vars_1 <- c(
  "country",
  "ci_divisio",
  "ci_sls_2")

vars_2 <- c(
  "new_or_c_1",
  "inter_cat",
  "interventi",
  "restoratio")


# Irrecoverable carbon 
# -----

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
      file = paste0("results/summaries/reports/ImpactIndicators_IC_", "_", str_to_title(user_group_1),
                    "_", str_to_title(user_group_2), ".csv")
    )
  }
  
}







# Other indicators
# -----

# Define summarizing function
other_summarize <- function(df){
  df %>% 
    summarize(across(
      .cols = c(area_ha, rest_area, population, 
                tstor_woody, tstor_soil, tstor_total, 
                carbon_seq_potl),
      sum, na.rm = TRUE),
      .groups = "keep") 
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
    overlaps_summary <- overlaps %>% 
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
    
    # combine the two data frames and summarize to subtract the overlaps
    corrected_summary <- sites_summary %>% 
      bind_rows(overlaps_summary) %>% 
      group_by_at(user_groups) %>% 
      other_summarize() %>% 
      ungroup() 
    
    write_csv(
      corrected_summary,
      file = paste0("results/summaries/reports/ImpactIndicators", "_", str_to_title(user_group_1),
                    "_", str_to_title(user_group_2), ".csv")
    )
  }
  
}
