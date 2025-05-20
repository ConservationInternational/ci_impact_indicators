# Summarize results from indicators analysis
# Contact: Anna Ballasiotes
# Last updated: 04/07/2025

library(tidyverse)

year <- "2024"

#####
##### Read in data
#####

#read in results tables as defined in script 2
#irrecoverable carbon broken up by ecosystem

sites_ic <- read_csv(
  paste0("results/FY",
         str_sub(year, start = 3),
         "_ImpactIndicators_IrrecoverableCarbon_Sites.csv")) %>%
  dplyr::select(!tonnes_ha_ic) # can't sum



overlaps_ic <- readRDS(
  paste0("results/FY",
         str_sub(year, start = 3),
         "_ImpactIndicators_IrrecoverableCarbon_Overlaps.rds")) %>%
  dplyr::select(!tonnes_ha_ic) # can't sum

#####
##### CORRECTIONS
#####

overlaps_ic <- overlaps_ic %>%
  rename(Area_ha_list = Area_ha)


# note fields are equal within overlapping polygons for country and ci_divsio, but not always for ci_divis_1 or sls s

####### FUNCTIONS

## IC SITES: "Area_ha"
## IC OVERLAPS: "Area_ha"


## Summarizing Function when adding to country area. 


## This is how we add the IC stuff otgetehr
ic_summarize <- function(df){
  df %>% 
    summarize(across(
      .cols = contains("_ic"),
      sum, na.rm = TRUE),
      .groups = "keep") 
}



#####
##### Tidy and define functions
#####

# Irrecoverable carbon
# create csv for aggregating all sites
# multiply value by (# of overlaps - 1) to account for triple/quadruple/etc counting 
# and make values negative to easily sum
# note we only use this for variables we know don't vary within overlapping sites - for example,
# all overlapping sites are within the same country, but they are not all within the same sls

## Make these negative

## So to correct for the overlaps,
## We need to make them "negative" ^^ see above. 
overlaps_ic_corrected <- overlaps_ic %>% 
  mutate(across(
    .cols = contains("_ic"),
    ~ as.numeric(.x) * -(n.overlaps - 1))) %>% 
  dplyr::select(!n.overlaps)


###########################
### RESTORATION SEPARATION
###########################

sites_ic_restor <- sites_ic %>%
  filter(Intervention_Category == "Restoration Areas")

write_csv(sites_ic_restor, 
          paste0("results/summaries_", year, "/ImpactIndicators_Sites_Restoration_IC.csv"))

### OVERLAPS:
overlaps_ic_corrected_restor <- overlaps_ic_corrected %>%
  filter(sapply(Intervention_Category, function(x) all(x == "Restoration Areas"))) %>%
  filter(sapply(Intervention_Type, function(x) any(duplicated(x))))


## BY RESTORATION TYPE:


# irrecoverable carbon

rest_type_sites_ic <- sites_ic_restor %>% 
  group_by(Intervention_Type) %>% 
  ic_summarize() %>% 
  ungroup()

## OVERLAPS
rest_type_overlaps_ic <- overlaps_ic_corrected_restor %>%
  rowwise() %>% 
  mutate(Intervention_Type = unique(Intervention_Type)) %>% # all countries within polygons are the same
  ungroup() %>% 
  group_by(Intervention_Type) %>% 
  ic_summarize() %>% 
  ungroup()

# Rest Type Summary
rest_type_corrected_ic <- rest_type_sites_ic %>% 
  bind_rows(rest_type_overlaps_ic) %>% 
  group_by(Intervention_Type) %>% 
  ic_summarize() %>%  # can sum bc overlap values are negative
  ungroup()

write_csv(rest_type_corrected_ic, 
          paste0("results/summaries_", year, "/ImpactIndicators_RestTypeSummary_IC.csv"))


#####
##### Country-level
#####

# irrecoverable carbon

country_sites_ic <- sites_ic %>% 
  group_by(Country) %>% 
  ic_summarize() %>% 
  ungroup()

## OVERLAPS
country_overlaps_ic <- overlaps_ic_corrected %>%
  rowwise() %>% 
  mutate(Country = unique(Country)) %>% # all countries within polygons are the same
  ungroup() %>% 
  group_by(Country) %>% 
  ic_summarize() %>% 
  ungroup()


# Country and Overlap Combo
country_corrected_ic <- country_sites_ic %>% 
  bind_rows(country_overlaps_ic) %>% 
  group_by(Country) %>% 
  ic_summarize() %>%  # can sum bc overlap values are negative
  ungroup()

write_csv(country_corrected_ic, 
          paste0("results/summaries_", year, "/ImpactIndicators_CountrySummary_IC.csv"))

#### irrecoverable carbon with ecosystem

country_sites_eco_ic <- sites_ic %>% 
  group_by(Country, ecosystem) %>% 
  ic_summarize() %>% 
  ungroup()

## OVERLAPS
country_overlaps_eco_ic <- overlaps_ic_corrected %>%
  rowwise() %>% 
  mutate(Country = unique(Country)) %>% # all countries within polygons are the same
  ungroup() %>% 
  group_by(Country, ecosystem) %>% 
  ic_summarize() %>% 
  ungroup()


# Country and Overlap Combo
country_corrected_eco_ic <- country_sites_eco_ic %>% 
  bind_rows(country_overlaps_eco_ic) %>% 
  group_by(Country, ecosystem) %>% 
  ic_summarize() %>%  # can sum bc overlap values are negative
  ungroup()

country_corrected_eco_ic_across <- country_corrected_eco_ic %>%
  group_by(Country, ecosystem) %>%
  pivot_wider(names_from = ecosystem, values_from = c(tstor_ic, ha_high_ic, ha_ic, tstor_blue_ic), names_sep = "_") 

write_csv(country_corrected_eco_ic, 
          paste0("results/summaries_", year, "/ImpactIndicators_CountrySummary_ecosystem_IC.csv"))


#####
##### New-Or-Old Level
#####

# IRRECOVERABLE CARBON

## By New/Old

new_sites_ic <- sites_ic %>% dplyr::filter(Site_Status == "New this year")

## I really just want to get the new sites, and... 
## Subtract any overlapping areas from the new sites only.
## So I should just extract out the sites that are NEW.
## And then perform a normal overlap subtraction, because
## I don't want any overlapping in general.
## And if the site still exists to have an overlap, then it will be 
## subtracted. If it doesn't, then no harm no foul

# 
overlaps_ic_corrected_new <- overlaps_ic_corrected %>%
  filter(map_lgl(Site_Status, ~ !any(grepl("Continued this year", .))))

## OVERLAPS
new_overlaps_ic <- overlaps_ic_corrected_new %>%
  rowwise() %>%
  mutate(Country = unique(Country)) %>% # all countries within polygons are the same
  ungroup() %>%
  group_by(Site_Status) %>%
  ic_summarize() %>%
  ungroup()

# Country and Overlap Combo
country_corrected_ic_new <- new_sites_ic %>% 
  bind_rows(new_overlaps_ic) %>% 
  group_by(Country) %>% 
  ic_summarize() %>%  # can sum bc overlap values are negative
  ungroup()

write_csv(country_corrected_ic_new, 
          paste0("results/summaries_", year, "/ImpactIndicators_CountrySummary_New_IC.csv"))


#####
##### Division-level
#####

# irrecoverable carbon

division_sites_ic <- sites_ic %>% 
  group_by(Division_Primary_CI, ecosystem) %>% 
  ic_summarize() %>% 
  ungroup()

division_overlaps_ic <- overlaps_ic_corrected %>%
  rowwise() %>% 
  mutate(Division_Primary_CI = unique(Division_Primary_CI)) %>% 
  ungroup() %>% 
  group_by(Division_Primary_CI, ecosystem) %>% 
  ic_summarize() %>% 
  ungroup()

division_corrected_ic <- division_sites_ic %>% 
  bind_rows(division_overlaps_ic) %>% 
  group_by(Division_Primary_CI, ecosystem) %>% 
  ic_summarize() %>% 
  ungroup()

write_csv(division_corrected_ic, 
          paste0("results/summaries_", year, "/ImpactIndicators_DivisionSummary_IC.csv"))

# other indicators

division_sites_other <- sites_other %>% 
  group_by(Division_Primary_CI) %>% 
  other_summarize() %>% 
  ungroup()

all_identical <- function(x) {
  length(unique(x)) == 1
}

#####
##### New-Old-level
#####

# New-Continued For IC By Ecosystem
# irrecoverable carbon

new_old_sites_ic <- sites_ic %>% 
  group_by(Site_Status, ecosystem) %>% 
  ic_summarize() %>% 
  ungroup()

## Ok so for a site that is... Continued This Year and overlapping with
## A site that is New This Year... how many? 
# only 3; i exlcuded those sites bc we don't want to subtract that to count new things


new_old_overlaps_ic <- overlaps_ic %>% 
  rowwise() %>% 
  # get a list of duplicated values
  mutate(duplicate_new_old = list(New_or_Continued_1[duplicated(New_or_Continued_1)])) %>% 
  ungroup() %>% 
  # remove those without any duplicates 
  # (i.e. birds head & west papua have some overlaps that aren't double counting when we group by sls)
  filter(!duplicate_new_old == "character(0)") %>% 
  # multiply values of interest by the number of times its duplicated & make negative 
  rowwise() %>% 
  mutate(across(
    .cols = contains("_ic"),
    ~ .x * -length(duplicate_new_old))) %>%
  anti_join(continuednew) %>%
  mutate(New_or_Continued_1 = unique(duplicate_new_old)[1]) %>% # is adding the "1" here correct
  # using duplicated[] already accounts for the -1
  ungroup() %>%
  group_by(New_or_Continued_1, ecosystem) %>%
  ic_summarize() %>%
  ungroup()

print(unique(new_old_overlaps_ic$duplicate_new_old))

continuednew <- new_old_overlaps_ic %>%
  filter("Continued this year" %in% duplicate_new_old & "New this year" %in% duplicate_new_old)


new_old_corrected_ic <- new_old_sites_ic %>% 
  bind_rows(new_old_overlaps_ic) %>% 
  group_by(New_or_Continued_1, ecosystem) %>% 
  ic_summarize() %>% 
  ungroup()

write_csv(new_old_corrected_ic, 
          paste0("results_v2/v2_summaries_", year, "/ImpactIndicators_NewContinuedSummary_IC.csv"))

#####

##### Site-level
#####

# overlaps don't occur at site level
# due to that, note can't sum for total (or will have double counting)

write_csv(sites_ic, 
          paste0("results/summaries_", year, "/ImpactIndicators_SiteSummary_IC.csv"))

