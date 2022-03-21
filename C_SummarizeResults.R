# Summarize results from indicators analysis
# Contact: Cameryn Brock
# Last updated: 3/20/2022

library(tidyverse)

#####

##### Read in data
#####

# read in results tables as defined in script B
# irrecoverable carbon broken up by ecosystem
sites_ic <- read_csv("results/FY21_ImpactIndicators_IrrecoverableCarbon_Sites.csv")
overlaps_ic <- read_csv("results/FY21_ImpactIndicators_IrrecoverableCarbon_Overlaps.csv")
overlaps_ic_fields <- readRDS("results/FY21_ImpactIndicators_IrrecoverableCarbon_Overlaps.rds")

# other indicators = population, woody & soil carbon, carbon sequestration potential
sites_other <- read_csv("results/FY21_ImpactIndicators_Other_Sites.csv")
overlaps_other <- read_csv("results/FY21_ImpactIndicators_Other_Overlaps.csv")
overlaps_other_fields <- readRDS("results/FY21_ImpactIndicators_Other_Overlaps.rds")

# note fields are equal within overlapping polygons for country and ci_divsio, but not always for ci_divis_1 or sls s


#####

##### Tidy and define functions
#####

# Irrecoverable carbon
# create csv for aggregating all sites
# multiply value by (# of overlaps - 1) to account for triple/quadruple/etc counting 
# and make values negative to easily sum
# note we only use this for variables we know don't vary within overlapping sites - for example, all overlapping sites are within the same country, but they are not all within the same sls
overlaps_ic_corrected <- overlaps_ic_fields %>% 
  mutate(across(
    .cols = contains("_ic"),
    ~ .x * -(n_overlaps - 1))) %>% 
  dplyr::select(!n_overlaps)

# define function for summarizing
ic_summarize <- function(df){
  df %>% 
    summarize(across(
      .cols = contains("_ic"),
      sum, na.rm = TRUE),
      .groups = "keep") 
}

overlaps_other_corrected <- overlaps_other_fields %>% 
  mutate(across(
    .cols = c(area_ha, rest_area, population, 
              tstor_woody, tstor_soil, tstor_total, 
              carbon_seq_potl),
    ~ .x * -(n_overlaps - 1))) %>% 
  dplyr::select(!n_overlaps)

other_summarize <- function(df){
  df %>% 
    summarize(across(
      .cols = c(area_ha, rest_area, population, 
                tstor_woody, tstor_soil, tstor_total, 
                carbon_seq_potl),
      sum, na.rm = TRUE),
      .groups = "keep") 
}



#####

##### Country-level
#####

# irrecoverable carbon

country_sites_ic <- sites_ic %>% 
  group_by(country) %>% 
  ic_summarize() %>% 
  ungroup()

country_overlaps_ic <- overlaps_ic_corrected %>%
  rowwise() %>% 
  mutate(country = unique(country)) %>% # all countries within polygons are the same
  ungroup() %>% 
  group_by(country) %>% 
  ic_summarize() %>% 
  ungroup()

country_corrected_ic <- country_sites_ic %>% 
  bind_rows(country_overlaps_ic) %>% 
  group_by(country) %>% 
  ic_summarize() %>%  # can sum bc overlap values are negative
  ungroup()

write_csv(country_corrected_ic, "results/summaries/ImpactIndicators_CountrySummary_IC.csv")

# other indicators

country_sites_other <- sites_other %>% 
  group_by(country) %>% 
  other_summarize() %>% 
  ungroup() 
  
country_overlaps_other <- overlaps_other_corrected %>%
  rowwise() %>% 
  mutate(country = unique(country)) %>% # all countries within polygons are the same
  ungroup() %>% 
  group_by(country) %>% 
  other_summarize() %>% 
  ungroup() 

country_corrected_other <- country_sites_other %>% 
  bind_rows(country_overlaps_other) %>% 
  group_by(country) %>% 
  other_summarize() %>%  # can sum bc overlap values are negative
  ungroup() 

write_csv(country_corrected_other, "results/summaries/ImpactIndicators_CountrySummary_OtherIndicators.csv")



#####

##### Division-level
#####

# irrecoverable carbon

division_sites_ic <- sites_ic %>% 
  group_by(ci_divisio, ecosystem) %>% 
  ic_summarize() %>% 
  ungroup()

division_overlaps_ic <- overlaps_ic_corrected %>%
  rowwise() %>% 
  mutate(ci_divisio = unique(ci_divisio)) %>% 
  ungroup() %>% 
  group_by(ci_divisio, ecosystem) %>% 
  ic_summarize() %>% 
  ungroup()

division_corrected_ic <- division_sites_ic %>% 
  bind_rows(division_overlaps_ic) %>% 
  group_by(ci_divisio, ecosystem) %>% 
  ic_summarize() %>% 
  ungroup()

write_csv(division_corrected_ic, "results/summaries/ImpactIndicators_DivisionSummary_IC.csv")

# other indicators

division_sites_other <- sites_other %>% 
  group_by(ci_divisio) %>% 
  other_summarize() %>% 
  ungroup()

division_overlaps_other <- overlaps_other_corrected %>%
  rowwise() %>% 
  mutate(ci_divisio = unique(ci_divisio)) %>% 
  ungroup() %>% 
  group_by(ci_divisio) %>% 
  other_summarize() %>% 
  ungroup()

division_corrected_other <- division_sites_other %>% 
  bind_rows(division_overlaps_other) %>% 
  group_by(ci_divisio) %>% 
  other_summarize() %>% 
  ungroup()

write_csv(division_corrected_other, "results/summaries/ImpactIndicators_DivisionSummary_OtherIndicators.csv")


#####

##### SlS-level
#####

# sls is more complicated because there are overlapping sites with different sls values

# irrecoverable carbon

sls_sites_ic <- sites_ic %>% 
  group_by(ci_sls_1, ecosystem) %>% 
  ic_summarize() %>% 
  ungroup()

sls_overlaps_ic <- overlaps_ic_fields %>% 
  rowwise() %>% 
  # get a list of duplicated values
  mutate(duplicate_sls = list(ci_sls_1[duplicated(ci_sls_1)])) %>% 
  ungroup() %>% 
  # remove those without any duplicates 
  # (i.e. birds head & west papua have some overlaps that aren't double counting when we group by sls)
  filter(!duplicate_sls == "character(0)") %>% 
  # multiply values of interest by the number of times its duplicated & make negative 
  rowwise() %>% 
  mutate(across(
    .cols = contains("_ic"),
    ~ .x * -length(duplicate_sls))) %>% 
  mutate(ci_sls_1 = unique(duplicate_sls)) %>% 
  # using duplicated[] already accounts for the -1
  ungroup() %>% 
  group_by(ci_sls_1, ecosystem) %>% 
  ic_summarize() %>% 
  ungroup()

sls_corrected_ic <- sls_sites_ic %>% 
  bind_rows(sls_overlaps_ic) %>% 
  group_by(ci_sls_1, ecosystem) %>% 
  ic_summarize() %>% 
  ungroup()

write_csv(sls_corrected_ic, "results/summaries/ImpactIndicators_SLSSummary_IC.csv")

# other indicators

sls_sites_other <- sites_other %>% 
  group_by(ci_sls_1) %>% 
  other_summarize() %>% 
  ungroup()

sls_overlaps_other <- overlaps_other_fields %>% 
  rowwise() %>% 
  mutate(duplicate_sls = list(ci_sls_1[duplicated(ci_sls_1)])) %>% 
  ungroup() %>% 
  filter(!duplicate_sls == "character(0)") %>% 
  rowwise() %>% 
  mutate(across(
    .cols = c(area_ha, rest_area, population, 
              tstor_woody, tstor_soil, tstor_total, 
              carbon_seq_potl),
    ~ .x * -length(duplicate_sls))) %>% 
  mutate(ci_sls_1 = unique(duplicate_sls)) %>% 
  ungroup() %>% 
  group_by(ci_sls_1) %>% 
  other_summarize() %>% 
  ungroup()

sls_corrected_other <- sls_sites_other %>% 
  bind_rows(sls_overlaps_other) %>% 
  group_by(ci_sls_1) %>% 
  other_summarize() %>% 
  ungroup()

write_csv(sls_corrected_other, "results/summaries/ImpactIndicators_SLSSummary_OtherIndicators.csv")

#####

##### Site-level
#####

# overlaps don't occur at site level
# due to that, note can't sum for total (or will have double counting)

write_csv(sites_ic, "results/summaries/ImpactIndicators_SiteSummary_IC.csv")
write_csv(sites_other, "results/summaries/ImpactIndicators_SiteSummary_OtherIndicators.csv")
