#7_Southern Cross Refresh


library(sf)
library(terra)
library(exactextractr)
library(janitor)
library(googledrive)
library(tidyverse)
library(purrr)
library(data.table)

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

other_summarize <- function(df){
  
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

########SUMMARY###########
# Protect - secured (new only), (improved)
# Protect - reduced (new + cont'd), (improved)

# Manage - reduced (new + cont'd)

# Manage - removed (new + cont'd)
  ## Agroforestry, Silvopasture, rangeland restoration

# Restore - removed 
###########################

setwd("C:/Users/aballasiotes/Dev/Projects/ci_impact_indicators/")

#######################################
########### LOAD DATA #################
#######################################

FY21_ImpactIndicators_IrrecoverableCarbon_Sites <- read_csv("results/previous/datafiles_FY21/source_data/FY21_ImpactIndicators_IrrecoverableCarbon_Sites.csv")
FY21_ImpactIndicators_Other_Sites <- read_csv("results/previous/datafiles_FY21/source_data/FY21_ImpactIndicators_Other_Sites.csv")

FY22_SouthAfrica_Rangeland_FY21Proxy <- read_csv("results/previous/datafiles_FY22/source_data/FY22_SouthAfrica_Rangeland_FY21Proxy.csv")
FY22_SouthAfrica_Rangeland_FY23Proxy <- read_csv("results/previous/datafiles_FY22/source_data/FY22_SouthAfrica_Rangeland_FY23Proxy.csv")

FY22_ImpactIndicators_IrrecoverableCarbon_Sites <- read_csv("results/previous/datafiles_FY22/source_data/FY22_ImpactIndicators_IrrecoverableCarbon_Sites.csv")
FY22_ImpactIndicators_Other_Sites <- read_csv("results/previous/datafiles_FY22/source_data/FY22_ImpactIndicators_Other_Sites.csv")

FY23_ImpactIndicators_IrrecoverableCarbon_Sites <- read_csv("results/previous/datafiles_FY23/source_data/FY23_ImpactIndicators_IrrecoverableCarbon_Sites_ALL.csv")
FY23_ImpactIndicators_Other_Sites <- read_csv("results/previous/datafiles_FY23/source_data/FY23_ImpactIndicators_OtherMetrics_Sites_ALL.csv")

FY24_ImpactIndicators_Other_Sites <- read_csv("results/previous/datafiles_FY24/source_data/FY24_ImpactIndicators_Primary_Sites_ALL.csv")
FY24_ImpactIndicators_IrrecoverableCarbon_Sites <- read_csv("results/previous/datafiles_FY24/source_data/FY24_ImpactIndicators_IrrecoverableCarbon_Sites.csv")


#######################################
########### CLEAN NAMES ###############
#######################################

setnames(FY21_ImpactIndicators_Other_Sites, old = c("emissions_avoided_mgco2e", "area_ha", "country", "interventi", "restoratio"), new = c("emissions_avoided_MgCO2e", "Area_ha", "Country", "Intervention_Type", "Restoration_Type"), skip_absent = TRUE)
setnames(FY22_ImpactIndicators_Other_Sites, old = c("emissions_avoided_mgco2e", "area_ha", "country", "interventi", "restoratio"), new = c("emissions_avoided_MgCO2e", "Area_ha", "Country", "Intervention_Type", "Restoration_Type"), skip_absent = TRUE)
setnames(FY23_ImpactIndicators_Other_Sites, old = c("emissions_avoided_mgco2e", "area_ha", "country", "Intervention_Type_Primary", "restoratio"), new = c("emissions_avoided_MgCO2e", "Area_ha", "Country", "Intervention_Type", "Restoration_Type"), skip_absent = TRUE)

setnames(FY22_SouthAfrica_Rangeland_FY21Proxy, old = c("emissions_avoided_mgco2e", "area_ha", "country", "interventi", "restoratio"), new = c("emissions_avoided_MgCO2e", "Area_ha", "Country", "Intervention_Type", "Restoration_Type"), skip_absent = TRUE)
setnames(FY22_SouthAfrica_Rangeland_FY23Proxy, old = c("emissions_avoided_mgco2e", "area_ha", "country", "Intervention_Type_Primary", "restoratio"), new = c("emissions_avoided_MgCO2e", "Area_ha", "Country", "Intervention_Type", "Restoration_Type"), skip_absent = TRUE)



setnames(FY21_ImpactIndicators_Other_Overlaps, c("area_ha", "country", "n_overlaps", "interventi", "restoratio"), c("Area_ha", "Country", "n.overlaps", "Intervention_Type", "Restoration_Type" ), skip_absent = TRUE)
setnames(FY22_ImpactIndicators_Other_Overlaps, c("area_ha", "country", "n_overlaps", "interventi", "restoratio"), c("Area_ha", "Country", "n.overlaps", "Intervention_Type", "Restoration_Type"), skip_absent = TRUE)
setnames(FY23_ImpactIndicators_Other_Overlaps, c("Area_ha", "area_ha", "country", "Intervention_Type_Primary", "restoratio"), c("Area_ha_overlaps", "Area_ha", "Country", "Intervention_Type", "Restoration_Type"), skip_absent = TRUE)
setnames(FY24_ImpactIndicators_Other_Overlaps, c("Area_ha", "Area_ha_R"), c("Area_ha_overlaps", "Area_ha"), skip_absent = TRUE)

setnames(FY21_ImpactIndicators_IrrecoverableCarbon_Sites, old = c("emissions_avoided_mgco2e", "area_ha", "country", "interventi"), new = c("emissions_avoided_MgCO2e", "Area_ha", "Country", "Intervention_Type"), skip_absent = TRUE)
setnames(FY22_ImpactIndicators_IrrecoverableCarbon_Sites, old = c("emissions_avoided_mgco2e", "area_ha", "country", "interventi"), new = c("emissions_avoided_MgCO2e", "Area_ha", "Country", "Intervention_Type"), skip_absent = TRUE)
setnames(FY23_ImpactIndicators_IrrecoverableCarbon_Sites, old = c("emissions_avoided_mgco2e", "area_ha", "country", "Intervention_Type_Primary"), new = c("emissions_avoided_MgCO2e", "Area_ha", "Country", "Intervention_Type"), skip_absent = TRUE)

setnames(FY21_ImpactIndicators_IrrecoverableCarbon_Overlaps, c("area_ha", "country", "n_overlaps", "interventi"), c("Area_ha", "Country", "n.overlaps", "Intervention_Type"), skip_absent = TRUE)
setnames(FY22_ImpactIndicators_IrrecoverableCarbon_Overlaps, c("area_ha", "country", "n_overlaps", "interventi"), c("Area_ha", "Country", "n.overlaps", "Intervention_Type"), skip_absent = TRUE)
setnames(FY23_ImpactIndicators_IrrecoverableCarbon_Overlaps, c("Area_ha", "area_ha", "country", "Intervention_Type_Primary"), c("Area_ha_overlaps", "Area_ha", "Country", "Intervention_Type"), skip_absent = TRUE)
setnames(FY24_ImpactIndicators_IrrecoverableCarbon_Overlaps, c("Area_ha", "Area_ha_R"), c("Area_ha_overlaps", "Area_ha"), skip_absent = TRUE)

#######################################
########### JOIN PRO/MA ###############
#######################################

FY21_ImpactIndicators_Other_Sites <- FY21_ImpactIndicators_Other_Sites %>% 
  left_join(ProtectManage, by = c("Intervention_Type"))
FY22_ImpactIndicators_Other_Sites <- FY22_ImpactIndicators_Other_Sites %>% 
  left_join(ProtectManage, by = c("Intervention_Type"))
FY23_ImpactIndicators_Other_Sites <- FY23_ImpactIndicators_Other_Sites %>% 
  left_join(ProtectManage, by = c("Intervention_Type"))
FY24_ImpactIndicators_Other_Sites <- FY24_ImpactIndicators_Other_Sites %>% 
  left_join(ProtectManage, by = c("Intervention_Type"))

FY21_ImpactIndicators_IrrecoverableCarbon_Sites <- FY21_ImpactIndicators_IrrecoverableCarbon_Sites %>% 
  left_join(ProtectManage, by = c("Intervention_Type"))
FY22_ImpactIndicators_IrrecoverableCarbon_Sites <- FY22_ImpactIndicators_IrrecoverableCarbon_Sites %>% 
  left_join(ProtectManage, by = c("Intervention_Type"))
FY23_ImpactIndicators_IrrecoverableCarbon_Sites <- FY23_ImpactIndicators_IrrecoverableCarbon_Sites %>% 
  left_join(ProtectManage, by = c("Intervention_Type"))
FY24_ImpactIndicators_IrrecoverableCarbon_Sites <- FY24_ImpactIndicators_IrrecoverableCarbon_Sites %>% 
  left_join(ProtectManage, by = c("Intervention_Type"))

##############################################
##############################################
############## SECURE CARBON #################
##############################################
##############################################

## PROTECT - SECURE: 

## filter for new <- apply IC analysis; By Intervention Type & Country

## FIX FY21 --> new + continued & improved
FY21_protect_secure <- FY21_ImpactIndicators_IrrecoverableCarbon_Sites %>%
  filter(new_or_c_1 == "New this year" | (new_or_c_1 == "Continued this year" &
           grepl("Yes", improved_m))) %>%
  filter(Protect_Manage == "Protect")

FY22_protect_secure <- FY22_ImpactIndicators_IrrecoverableCarbon_Sites %>%
  filter(new_or_c_1 == "New this year") %>%
  filter(Protect_Manage == "Protect")

FY23_protect_secure <- FY23_ImpactIndicators_IrrecoverableCarbon_Sites %>%
  filter(New_or_Continued_1 == "New this year") %>%
  filter(Protect_Manage == "Protect")
  
FY24_protect_secure <- FY24_ImpactIndicators_IrrecoverableCarbon_Sites %>%
  filter(Site_Status == "New this year") %>%
  filter(Protect_Manage == "Protect")

## PROTECT - SECURE OVERLAPS:
FY21_protect_secure_overlaps <- FY21_ImpactIndicators_IrrecoverableCarbon_Overlaps %>%
  filter(map_lgl(ci_id, ~all(.x %in% FY21_protect_secure$ci_id))) %>%
  filter(map_lgl(Country, ~ length(unique(.x)) == 1))

FY22_protect_secure_overlaps <- FY22_ImpactIndicators_IrrecoverableCarbon_Overlaps %>%
  filter(map_lgl(ci_id, ~all(.x %in% FY22_protect_secure$ci_id))) %>%
  filter(map_lgl(Country, ~ length(unique(.x)) == 1))

FY23_protect_secure_overlaps <- FY23_ImpactIndicators_IrrecoverableCarbon_Overlaps %>%
  filter(map_lgl(CI_ID, ~all(.x %in% FY23_protect_secure$CI_ID))) %>%
  filter(map_lgl(Country, ~ length(unique(.x)) == 1))

FY24_protect_secure_overlaps <- FY24_ImpactIndicators_IrrecoverableCarbon_Overlaps %>%
  filter(map_lgl(CI_ID, ~all(.x %in% FY24_protect_secure$CI_ID))) %>%
  filter(map_lgl(Country, ~ length(unique(.x)) == 1))


# define functions
ic_summarize <- function(df){
  df %>% 
    summarize(across(
      .cols = contains("_ic"),
      sum, na.rm = TRUE),
      .groups = "keep") 
}


years <- c("FY21", "FY22", "FY23", "FY24")

groups <- c("protect_secure")

vars_1 <- c(
  "Country")

vars_2 <- c(
  "Intervention_Type")

for (y in seq_along(years)) { 
  
  for (n in seq_along(groups)) {
    group_name <- paste0(years[y], "_", groups[n])
    overlaps_name <- paste0(years[y], "_", groups[n], "_overlaps")
    
    sites_ic_df <- get(group_name)
    overlaps_ic_df <- get(overlaps_name)
    
    for (i in seq_along(vars_1)) {
      
      user_group_1 <- vars_1[i]
      
      for (g in seq_along(vars_2)) {
        
        user_group_2 <- vars_2[g]
        
        user_groups <- c(user_group_1, user_group_2)
        
        sites_ic_summary <- sites_ic_df %>% 
          group_by_at(user_groups) %>% 
          ic_summarize() %>% 
          ungroup()
        
      # summarize overlaps
      overlaps_ic_summary <- overlaps_ic_df %>% 
        dplyr::select((!!as.symbol(user_group_1)), (!!as.symbol(user_group_2)),
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
        # summarize as you did with sites
        group_by_at(user_groups) %>% 
        ic_summarize() %>% 
        ungroup() 
      
     
      # combine the two data frames and summarize to subtract the overlaps_ic
      corrected_ic_summary <- sites_ic_summary %>% 
        bind_rows(overlaps_ic_summary) %>% 
        group_by_at(user_groups) %>% 
        ic_summarize() %>% 
        ungroup() %>%
        left_join(CountryRegion, by = "Country")
      
      write_csv(
        corrected_ic_summary,
        file = paste0("results/SCto2030/", str_to_title(group_name), "_", str_to_title(user_group_1), ".csv")
      )
    }
    
  }
  }
}

write_csv(
  sites_ic_summary,
  file = paste0("results/SCto2030/FY21_", str_to_title(user_group_1), ".csv")
)
##############################################
##############################################
############## REDUCE CARBON #################
##############################################
##############################################

#######################
## PROTECT - REDUCE: ##
#######################

## CHANGE so that it is NEW OR Continued & Improved ##

FY21_protect_reduce <- FY21_ImpactIndicators_Other_Sites %>%
  filter(new_or_c_1 == "New this year" | new_or_c_1 == "Continued this year" & grepl("Yes",improved_m)) %>%
  filter(Protect_Manage == "Protect")

FY22_protect_reduce <- FY22_ImpactIndicators_Other_Sites %>%
  filter(new_or_c_1 == "New this year" | new_or_c_1 == "Continued this year" & grepl("Yes",improved_m)) %>%
  filter(Protect_Manage == "Protect")

FY23_protect_reduce <- FY23_ImpactIndicators_Other_Sites %>%
  filter(New_or_Continued_1 == "New this year" | New_or_Continued_1 == "Continued this year" & grepl("Yes",Improved_Management)) %>%
  filter(Protect_Manage == "Protect")

FY24_protect_reduce <- FY24_ImpactIndicators_Other_Sites %>%
  filter(Site_Status == "New this year" | Site_Status == "Continued this year" & grepl("Yes",Improved_Management)) %>%
  filter(Protect_Manage == "Protect")


## PROTECT - REDUCE OVERLAPS:
FY21_protect_reduce_overlaps <- FY21_ImpactIndicators_Other_Overlaps %>%
  filter(map_lgl(ci_id, ~all(.x %in% FY21_protect_reduce$ci_id)))

FY22_protect_reduce_overlaps <- FY22_ImpactIndicators_Other_Overlaps %>%
  filter(map_lgl(ci_id, ~all(.x %in% FY22_protect_reduce$ci_id)))

FY23_protect_reduce_overlaps <- FY23_ImpactIndicators_Other_Overlaps %>%
  filter(map_lgl(CI_ID, ~all(.x %in% FY23_protect_reduce$CI_ID)))

FY24_protect_reduce_overlaps <- FY24_ImpactIndicators_Other_Overlaps %>%
  filter(map_lgl(CI_ID, ~all(.x %in% FY24_protect_reduce$CI_ID)))

#######################
## MANAGE - REDUCE: ##
#######################

## filter for New & Continued <-- apply Other analysis; By Intervention Type, Country, & AvoidAdd

FY21_manage_reduce <- FY21_ImpactIndicators_Other_Sites %>%
  filter(new_or_c_1 == "New this year" | new_or_c_1 == "Continued this year") %>%
  filter(Protect_Manage == "Manage")

FY22_manage_reduce <- FY22_ImpactIndicators_Other_Sites %>%
  filter(new_or_c_1 == "New this year" | new_or_c_1 == "Continued this year") %>%
  filter(Protect_Manage == "Manage")

FY23_manage_reduce <- FY23_ImpactIndicators_Other_Sites %>%
  filter(New_or_Continued_1 == "New this year" | New_or_Continued_1 == "Continued this year") %>%
  filter(Protect_Manage == "Manage")

FY24_manage_reduce <- FY24_ImpactIndicators_Other_Sites %>%
  filter(Site_Status == "New this year" | Site_Status == "Continued this year") %>%
  filter(Protect_Manage == "Manage")


# CORRECT OVERLAP INTERVENTIONS

FY21_manage_reduce_overlaps <- FY21_ImpactIndicators_Other_Overlaps %>%
  filter(map_lgl(Intervention_Type, ~ length(unique(.x)) == 1)) %>%
  filter(map_lgl(ci_id, ~all(.x %in% FY21_manage_reduce$ci_id)))

FY22_manage_reduce_overlaps <- FY22_ImpactIndicators_Other_Overlaps %>%
  filter(map_lgl(Intervention_Type, ~ length(unique(.x)) == 1)) %>%
  filter(map_lgl(ci_id, ~all(.x %in% FY22_manage_reduce$ci_id)))

FY23_manage_reduce_overlaps <- FY23_ImpactIndicators_Other_Overlaps %>%
  filter(map_lgl(Intervention_Type, ~ length(unique(.x)) == 1)) %>%
  filter(map_lgl(CI_ID, ~all(.x %in% FY23_manage_reduce$CI_ID)))

FY24_manage_reduce_overlaps <- FY24_ImpactIndicators_Other_Overlaps %>%
  filter(map_lgl(Intervention_Type, ~ length(unique(.x)) == 1)) %>%
  filter(map_lgl(CI_ID, ~all(.x %in% FY24_manage_reduce$CI_ID)))

years <- c("FY21", "FY22", "FY23", "FY24")

groups <- c("protect_reduce")

vars_1 <- c(
  "Country")

vars_2 <- c(
  "Intervention_Type")

for (y in seq_along(years)) { 
  
  for (n in seq_along(groups)) {
    group_name <- paste0(years[y], "_", groups[n])
    overlaps_name <- paste0(years[y], "_", groups[n], "_overlaps")
    
    sites_df <- get(group_name)
    overlaps_df <- get(overlaps_name)
    
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
          primary_summarize() %>% 
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
          ungroup() %>%
          left_join(CountryRegion, by = "Country")
        
        write_csv(
          corrected_summary,
          file = paste0("results/SCto2030/", str_to_title(group_name), "_", str_to_title(user_group_1), ".csv")
        ) 
      }
    }
  }
  }

##############################################
##############################################
############## REMOVE CARBON #################
##############################################
##############################################
## parse dates

FY23_ImpactIndicators_Other_Sites$CI_Start_Date <- parse_date_time(FY23_ImpactIndicators_Other_Sites$CI_Start_Date, orders = c("dmy", "ymd", "mdy"))
FY23_ImpactIndicators_Other_Sites$Tree_Planting_Date <- parse_date_time(FY23_ImpactIndicators_Other_Sites$Tree_Planting_Date, orders = c("dmy", "ymd", "mdy"))

FY24_ImpactIndicators_Other_Sites$CI_Start_Date <- parse_date_time(FY24_ImpactIndicators_Other_Sites$CI_Start_Date, orders = c("dmy", "ymd", "mdy"))
FY24_ImpactIndicators_Other_Sites$Tree_Planting_Date <- parse_date_time(FY24_ImpactIndicators_Other_Sites$Tree_Planting_Date, orders = c("dmy", "ymd", "mdy"))


######## 
FY21_restore_remove <- FY21_ImpactIndicators_Other_Sites %>%
  filter(new_or_c_1 == "New this year" | new_or_c_1 == "Continued this year") %>%
  filter(!Restoration_Type == "Not Applicable") %>%
  mutate(remove = case_when(
    tree_plant > as.Date('2018-01-01') & 
      tree_plant < as.Date('2021-06-30') ~ "keep",
    Restoration_Type %in% c("Rangeland Restoration - Planned Grazing", "Natural Regeneration") &
      ci_start_d > as.Date('2017-12-31') & 
      ci_start_d < as.Date('2021-06-30') ~ "keep",
    T ~ "remove"
  ))


FY21_restore_remove_overlaps <- FY21_ImpactIndicators_Other_Overlaps %>%
  filter(map_lgl(Restoration_Type, ~ length(unique(.x)) == 1)) %>%
  filter(map_lgl(ci_id, ~all(.x %in% FY21_restore_remove$ci_id)))

############ 

FY22_restore_remove <- FY22_ImpactIndicators_Other_Sites %>%
  filter(new_or_c_1 == "New this year" | new_or_c_1 == "Continued this year") %>%
  filter(!Restoration_Type == "Not Applicable") %>%
  mutate(remove = case_when(
    tree_plant > as.Date('2018-01-01') & 
      tree_plant < as.Date('2022-06-30') ~ "keep",
    Restoration_Type %in% c("Rangeland Restoration", "Natural Regeneration") &
      ci_start_d > as.Date('2018-01-01') & 
      ci_start_d < as.Date('2022-06-30') ~ "keep",
    T ~ "remove"
  )) %>%
  filter(remove == "keep")


FY22_restore_remove_overlaps <- FY22_ImpactIndicators_Other_Overlaps %>%
  filter(map_lgl(Restoration_Type, ~ length(unique(.x)) == 1)) %>%
  filter(map_lgl(ci_id, ~all(.x %in% FY22_restore_remove$ci_id)))

############

FY23_restore_remove <- FY23_ImpactIndicators_Other_Sites %>%
  filter(New_or_Continued_1 == "New this year" | New_or_Continued_1 == "Continued this year") %>%
  filter(!Restoration_Type == "Not Applicable") %>%
  mutate(remove = case_when(
    Tree_Planting_Date > as.Date('2018-01-01') & 
      Tree_Planting_Date < as.Date('2024-06-30') ~ "keep",
    Restoration_Type %in% c("Rangeland Restoration - Planned Grazing", "Natural Regeneration") &
      CI_Start_Date > as.Date('2017-12-31') & 
      CI_Start_Date < as.Date('2024-06-30') ~ "keep",
    T ~ "remove"
  )) %>% 
  filter(remove == "keep")


FY23_restore_remove_overlaps <- FY23_ImpactIndicators_Other_Overlaps %>%
  filter(map_lgl(Restoration_Type, ~ length(unique(.x)) == 1)) %>%
  filter(map_lgl(CI_ID, ~all(.x %in% FY23_restore_remove$CI_ID)))


FY24_restore_remove <- FY24_ImpactIndicators_Other_Sites %>%
  filter(Site_Status == "New this year" | Site_Status == "Continued this year") %>%
  filter(grepl("Restoration", Intervention_Type)) %>%
  mutate(remove = case_when(
    Tree_Planting_Date > as.Date('2018-01-01') & 
      Tree_Planting_Date < as.Date('2024-06-30') ~ "keep",
    Intervention_Type %in% c("Rangeland Restoration - Planned Grazing", "Restoration - Natural Regeneration") &
      CI_Start_Date > as.Date('2018-01-01') & 
      CI_Start_Date < as.Date('2024-06-30') ~ "keep",
    T ~ "remove"
  )) %>%
  filter(remove == "keep")


FY24_restore_remove_overlaps <- FY24_ImpactIndicators_Other_Overlaps %>%
  filter(map_lgl(Intervention_Type, ~ length(unique(.x)) == 1)) %>%
  filter(map_lgl(CI_ID, ~all(.x %in% FY24_restore_remove$CI_ID)))

##########
## Weird ZAF Rangeland Proxies
##########

FY21_ZAF_proxy_restore_remove <- FY22_SouthAfrica_Rangeland_FY21Proxy

FY21_ZAF_proxy_restore_remove_overlaps <-  FY21_ImpactIndicators_Other_Overlaps %>%
  filter(map_lgl(Restoration_Type, ~ length(unique(.x)) == 1)) %>%
  filter(map_lgl(ci_id, ~all(.x %in% FY21_ZAF_proxy_restore_remove$ci_id)))


FY23_ZAF_proxy_restore_remove <- FY22_SouthAfrica_Rangeland_FY23Proxy

FY23_ZAF_proxy_restore_remove_overlaps <- FY23_ImpactIndicators_Other_Overlaps %>%
  filter(map_lgl(Restoration_Type, ~ length(unique(.x)) == 1)) %>%
  filter(map_lgl(CI_ID, ~all(.x %in% FY23_ZAF_proxy_restore_remove$CI_ID)))


years <- c("FY21","FY23")


groups <- c("ZAF_proxy_restore_remove")

vars_1 <- c(
  "Country")

vars_2 <- c(
  "Restoration_Type")

years <- c("FY24")

vars_2 <- c(
  "Intervention_Type")

for (y in seq_along(years)) { 
  
  for (n in seq_along(groups)) {
    group_name <- paste0(years[y], "_", groups[n])
    overlaps_name <- paste0(years[y], "_", groups[n], "_overlaps")
    
    sites_df <- get(group_name)
    overlaps_df <- get(overlaps_name)
    
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
        overlaps_summary <- overlaps_df %>% 
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
          primary_summarize() %>% 
          ungroup() 
        
        # combine the two data frames and summarize to subtract the overlaps
        corrected_summary <- sites_summary %>% 
          bind_rows(overlaps_summary) %>% 
          group_by_at(user_groups) %>% 
          primary_summarize() %>% 
          ungroup() %>%
          left_join(CountryRegion, by = "Country")
        
        write_csv(
          corrected_summary,
          file = paste0("results/SCto2030/", str_to_title(group_name), "_", str_to_title(user_group_1), ".csv")
        ) 
      }
    }
  }
}































# 
# ## filter for Sequestration Dates <-- apply Other analysis; By Intervention Type, Country
# 
# ## FY21 & FY22
# year <- "FY22"
# user_group_1 <- "country"
# user_group_2 <- "restoratio"
# 
# user_groups <- c(user_group_1, user_group_2)
# 
# sites_summary <- FY22_ImpactIndicators_Other_Sites %>% 
#   # filter(restoratio != "Not Applicable") %>% 
#   # # remove sites as requested
#   # mutate(remove = case_when(
#   #   tree_plant > as.Date('2018-01-01') & 
#   #     tree_plant < as.Date('2021-06-30') ~ "keep",
#   #   restoratio %in% c("Rangeland Restoration - Planned Grazing", "Natural Regeneration") &
#   #     ci_start_d > as.Date('2018-01-01') & 
#   #     ci_start_d < as.Date('2021-06-30') ~ "keep",
#   #   T ~ "remove"
#   # )) %>% 
#   # filter(remove == "keep") %>% 
#   # dplyr::select(!remove) %>% 
#   group_by_at(user_groups) %>% 
#   other_summarize() %>% 
#   ungroup() 
# 
# remove_ids <- FY22_ImpactIndicators_Other_Sites %>% 
#   # remove sites as requested
#   mutate(remove = case_when(
#     tree_plant > as.Date('2018-01-01') & 
#       tree_plant < as.Date('2022-06-30') ~ "keep",
#     restoratio %in% c("Rangeland Restoration - Planned Grazing", "Natural Regeneration") &
#       ci_start_d > as.Date('2018-01-01') & 
#       ci_start_d < as.Date('2022-06-30') ~ "keep",
#     T ~ "remove"
#   )) %>% 
#   filter(remove == "remove") %>% 
#   pluck("ci_id") %>% 
#   unique()
# 
# # summarize overlaps
# overlaps_summary <- FY22_ImpactIndicators_Other_Overlaps %>% 
#   rowwise() %>% 
#   mutate(remove = any(unlist(ci_id) %in% remove_ids)) %>%
#   ungroup() %>% 
#   filter(!remove == TRUE) %>% 
#   dplyr::select((!!as.symbol(user_group_1)), (!!as.symbol(user_group_2)), 
#                 area_ha, rest_area, population, tstor_woody, 
#                 tstor_soil, tstor_total, carbon_seq_potl) %>% 
#   rowwise() %>% 
#   # get a list of duplicated values
#   mutate(duplicates = list((!!as.symbol(user_group_2))[duplicated((!!as.symbol(user_group_2)))])) %>% 
#   ungroup() %>% 
#   # remove those with no duplicates
#   filter(!duplicates == "character(0)") %>% 
#   # multiply values of interest by the number of times its duplicated & make negative 
#   rowwise() %>%
#   mutate({{ user_group_2 }} := list(unique(duplicates))) %>% 
#   unnest({{ user_group_2 }}) %>% 
#   ungroup() %>% 
#   rowwise() %>% 
#   mutate(duplicate_n = sum((!!as.symbol(user_group_2)) == duplicates, na.rm = TRUE)) %>% 
#   mutate(across(
#     .cols = c(area_ha, rest_area, population, tstor_woody, 
#               tstor_soil, tstor_total, carbon_seq_potl),
#     ~ .x * -duplicate_n))  %>% 
#   
#   # do again for second grouping if needed
#   mutate(duplicates = list((!!as.symbol(user_group_1))[duplicated((!!as.symbol(user_group_1)))])) %>% 
#   ungroup() %>% 
#   filter(!duplicates == "character(0)") %>% 
#   rowwise() %>%
#   mutate({{ user_group_1 }} := unique(duplicates)) %>% 
#   #      unnest({{ user_group_1 }}) %>% 
#   ungroup() %>% 
#   # summarize as you did with sites
#   group_by_at(user_groups) %>% 
#   other_summarize() %>% 
#   ungroup() %>% 
#   filter(!restoratio == "Not Applicable")
# 
# # combine the two data frames and summarize to subtract the overlaps
# corrected_summary <- sites_summary %>% 
#   bind_rows(overlaps_summary) %>% 
#   group_by_at(user_groups) %>% 
#   other_summarize() %>% 
#   ungroup() %>% 
#   dplyr::select(c(1,2,9)) %>%
#   left_join(CountryRegion, by = c("country" = "Country"))
# 
# write_csv(
# corrected_summary,
# file = paste0("results/SCto2030/", year, "ImpactIndicators", "_", str_to_title(user_group_1),
#               "_", str_to_title(user_group_2), "_carbonSeqPost2018.csv")
# )
# 
#   
# ## FY23
# year <- "FY23"
# user_group_1 <- "Country"
# user_group_2 <- "Restoration_Type"
# 
# FY23_restore_remove <- FY23_ImpactIndicators_Other_Sites %>%
#   filter(grepl("Restoration", Intervention_Type))
# 
# user_groups <- c(user_group_1, user_group_2)
# 
# sites_summary <- FY23_ImpactIndicators_Other_Sites %>% 
#   # filter(restoratio != "Not Applicable") %>% 
#   # # remove sites as requested
#   # mutate(remove = case_when(
#   #   tree_plant > as.Date('2018-01-01') & 
#   #     tree_plant < as.Date('2021-06-30') ~ "keep",
#   #   restoratio %in% c("Rangeland Restoration - Planned Grazing", "Natural Regeneration") &
#   #     ci_start_d > as.Date('2018-01-01') & 
#   #     ci_start_d < as.Date('2021-06-30') ~ "keep",
#   #   T ~ "remove"
#   # )) %>% 
#   # filter(remove == "keep") %>% 
#   # dplyr::select(!remove) %>% 
#   group_by_at(user_groups) %>% 
#   other_summarize() %>% 
#   ungroup() 
# 
# remove_ids <- FY23_ImpactIndicators_Other_Sites %>% 
#   # remove sites as requested
#   mutate(remove = case_when(
#     Tree_Planting_Date > as.Date('2018-01-01') & 
#       Tree_Planting_Date < as.Date('2023-06-30') ~ "keep",
#     Restoration_Type %in% c("Rangeland Restoration - Planned Grazing", "Natural Regeneration") &
#       CI_Start_Date > as.Date('2018-01-01') & 
#       CI_Start_Date < as.Date('2023-06-30') ~ "keep",
#     T ~ "remove"
#   )) %>% 
#   filter(remove == "remove") %>% 
#   pluck("CI_ID") %>% 
#   unique()
# 
# 
# # summarize overlaps
# overlaps_summary <- FY23_ImpactIndicators_Other_Overlaps %>% 
#   rowwise() %>% 
#   mutate(remove = any(unlist(CI_ID) %in% remove_ids)) %>%
#   ungroup() %>% 
#   filter(!remove == TRUE) %>% 
#   dplyr::select((!!as.symbol(user_group_1)), (!!as.symbol(user_group_2)), 
#                 Area_ha, rest_area, population, tstor_woody, 
#                 tstor_soil, tstor_total, carbon_seq_potl) %>% 
#   rowwise() %>% 
#   # get a list of duplicated values
#   mutate(duplicates = list((!!as.symbol(user_group_2))[duplicated((!!as.symbol(user_group_2)))])) %>% 
#   ungroup() %>% 
#   # remove those with no duplicates
#   filter(!duplicates == "character(0)") %>% 
#   # multiply values of interest by the number of times its duplicated & make negative 
#   rowwise() %>%
#   mutate({{ user_group_2 }} := list(unique(duplicates))) %>% 
#   unnest({{ user_group_2 }}) %>% 
#   ungroup() %>% 
#   rowwise() %>% 
#   mutate(duplicate_n = sum((!!as.symbol(user_group_2)) == duplicates, na.rm = TRUE)) %>% 
#   mutate(across(
#     .cols = c(Area_ha, rest_area, population, tstor_woody, 
#               tstor_soil, tstor_total, carbon_seq_potl),
#     ~ .x * -duplicate_n))  %>% 
#   
#   # do again for second grouping if needed
#   mutate(duplicates = list((!!as.symbol(user_group_1))[duplicated((!!as.symbol(user_group_1)))])) %>% 
#   ungroup() %>% 
#   filter(!duplicates == "character(0)") %>% 
#   rowwise() %>%
#   mutate({{ user_group_1 }} := unique(duplicates)) %>% 
#   #      unnest({{ user_group_1 }}) %>% 
#   ungroup() %>% 
#   # summarize as you did with sites
#   group_by_at(user_groups) %>% 
#   other_summarize() %>% 
#   ungroup() %>% 
#   filter(!Restoration_Type == "Not Applicable")
# 
# # combine the two data frames and summarize to subtract the overlaps
# corrected_summary <- sites_summary %>% 
#   bind_rows(overlaps_summary) %>% 
#   group_by_at(user_groups) %>% 
#   other_summarize() %>% 
#   ungroup() %>% 
#   dplyr::select(c(1,2,9)) %>%
#   left_join(CountryRegion, by = c("Country"))
# 
# write_csv(
#   corrected_summary,
#   file = paste0("results/SCto2030/", year, "ImpactIndicators", "_", str_to_title(user_group_1),
#                 "_", str_to_title(user_group_2), "_carbonSeqPost2018.csv")
# )
# 
# 
# ## FY24
# year <- "FY24"
# user_group_1 <- "Country"
# user_group_2 <- "Intervention_Type"
# 
# 
# FY24_ImpactIndicators_Other_Sites$CI_End_Date <- parse_date_time(FY24_ImpactIndicators_Other_Sites$CI_End_Date, orders = c("dmy", "ymd", "mdy"))
# FY24_ImpactIndicators_Other_Sites$CI_Start_Date <- parse_date_time(FY24_ImpactIndicators_Other_Sites$CI_Start_Date, orders = c("dmy", "ymd", "mdy"))
# FY24_ImpactIndicators_Other_Sites$Gazettement_Date <- parse_date_time(FY24_ImpactIndicators_Other_Sites$Gazettement_Date, orders = c("dmy", "ymd", "mdy"))
# FY24_ImpactIndicators_Other_Sites$Tree_Planting_Date <- parse_date_time(FY24_ImpactIndicators_Other_Sites$Tree_Planting_Date, orders = c("dmy", "ymd", "mdy"))
# 
# FY24_restore_remove <- FY24_ImpactIndicators_Other_Sites %>%
#   filter(grepl("Restoration", Intervention_Type))
#   
# 
# user_groups <- c(user_group_1, user_group_2)
# 
# sites_summary <- FY24_ImpactIndicators_Other_Sites_restor %>% 
#   # filter(restoratio != "Not Applicable") %>% 
#   # # remove sites as requested
#   # mutate(remove = case_when(
#   #   tree_plant > as.Date('2018-01-01') & 
#   #     tree_plant < as.Date('2021-06-30') ~ "keep",
#   #   restoratio %in% c("Rangeland Restoration - Planned Grazing", "Natural Regeneration") &
#   #     ci_start_d > as.Date('2018-01-01') & 
#   #     ci_start_d < as.Date('2021-06-30') ~ "keep",
#   #   T ~ "remove"
#   # )) %>% 
#   # filter(remove == "keep") %>% 
#   # dplyr::select(!remove) %>% 
#   group_by_at(user_groups) %>% 
#   other_summarize() %>% 
#   ungroup() 
# 
# remove_ids <- FY24_ImpactIndicators_Other_Sites_restor %>% 
#   # remove sites as requested
#   mutate(remove = case_when(
#     Tree_Planting_Date > as.Date('2018-01-01') & 
#       Tree_Planting_Date < as.Date('2024-06-30') ~ "keep",
#     Intervention_Type %in% c("Rangeland Restoration - Planned Grazing", "Restoration - Natural Regeneration") &
#       CI_Start_Date > as.Date('2018-01-01') & 
#       CI_Start_Date < as.Date('2024-06-30') ~ "keep",
#     T ~ "remove"
#   )) %>% 
#   filter(remove == "remove") %>% 
#   pluck("CI_ID") %>% 
#   unique()
# 
# FY24_ImpactIndicators_Other_Overlaps$Area_ha <- FY24_ImpactIndicators_Other_Overlaps$Area_ha_R
# 
# # summarize overlaps
# overlaps_summary <- FY24_ImpactIndicators_Other_Overlaps %>% 
#   rowwise() %>% 
#   mutate(remove = any(unlist(CI_ID) %in% remove_ids)) %>%
#   ungroup() %>% 
#   filter(!remove == TRUE) %>% 
#   dplyr::select((!!as.symbol(user_group_1)), (!!as.symbol(user_group_2)), 
#                 Area_ha, rest_area, population, tstor_woody, 
#                 tstor_soil, tstor_total, carbon_seq_potl) %>% 
#   rowwise() %>% 
#   # get a list of duplicated values
#   mutate(duplicates = list((!!as.symbol(user_group_2))[duplicated((!!as.symbol(user_group_2)))])) %>% 
#   ungroup() %>% 
#   # remove those with no duplicates
#   filter(!duplicates == "character(0)") %>% 
#   # multiply values of interest by the number of times its duplicated & make negative 
#   rowwise() %>%
#   mutate({{ user_group_2 }} := list(unique(duplicates))) %>% 
#   unnest({{ user_group_2 }}) %>% 
#   ungroup() %>% 
#   rowwise() %>% 
#   mutate(duplicate_n = sum((!!as.symbol(user_group_2)) == duplicates, na.rm = TRUE)) %>% 
#   mutate(across(
#     .cols = c(Area_ha, rest_area, population, tstor_woody, 
#               tstor_soil, tstor_total, carbon_seq_potl),
#     ~ .x * -duplicate_n))  %>% 
#   
#   # do again for second grouping if needed
#   mutate(duplicates = list((!!as.symbol(user_group_1))[duplicated((!!as.symbol(user_group_1)))])) %>% 
#   ungroup() %>% 
#   filter(!duplicates == "character(0)") %>% 
#   rowwise() %>%
#   mutate({{ user_group_1 }} := unique(duplicates)) %>% 
#   #      unnest({{ user_group_1 }}) %>% 
#   ungroup() %>% 
#   # summarize as you did with sites
#   group_by_at(user_groups) %>% 
#   other_summarize() %>% 
#   ungroup()
# 
# # combine the two data frames and summarize to subtract the overlaps
# corrected_summary <- sites_summary %>% 
#   bind_rows(overlaps_summary) %>% 
#   group_by_at(user_groups) %>% 
#   other_summarize() %>% 
#   ungroup() %>% 
#   dplyr::select(c(1,2,9)) %>%
#   left_join(CountryRegion, by = c("Country"))
# 
# write_csv(
#   corrected_summary,
#   file = paste0("results/SCto2030/", year, "ImpactIndicators", "_", str_to_title(user_group_1),
#                 "_", str_to_title(user_group_2), "_carbonSeqPost2018.csv")
# )
