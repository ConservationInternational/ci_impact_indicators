# Summarize results from indicators analysis
# Contact: Anna Ballasiotes
# Last updated: 04/07/2025

library(tidyverse)

year <- "2024"


#####
##### Read in data
#####

#read in results tables as defined in script 2

sites_primary <- read_csv(
  paste0("results/FY", str_sub(year, start = 3), "_ImpactIndicators_Primary_Sites_ALL.csv"))

overlaps_primary <- readRDS(
  paste0("results/FY", str_sub(year, start = 3), "_ImpactIndicators_Primary_Overlaps.rds"))


### avoided emissions - by site, already prepped
ae <- readRDS("data/avoided_emissions/ae_by_site_by_year_2024.rds") %>% 
  dplyr::select(CI_ID, emissions_avoided_MgCO2e_2023)


#####
##### CORRECTIONS
#####

# NOTE: I commented this out, but these are some manual corrections
# that I have made in years past; you can use some of this code if you
# also need to make manual corrections. 

# sites_primary <- sites_primary %>%
#   filter(!CI_ID %in% c("PHL1105", "PHL1106", "PHL1107", "PHL1108", "PHL1109", "PHL1110"))
#   
# overlaps_primary <- overlaps_primary %>%
#   filter(!map_chr(CI_ID, ~ .x[1]) %in% c("BNA1022", "BNA1007", "BNA1057", "BNA1073", "BNA1024", "BNA1052")) %>%
#   filter(!map_chr(CI_ID, ~ .x[2]) %in% c("BNA1018", "PER1044")) %>%
#   filter(!map_lgl(CI_ID, ~ any(.x %in% c("PHL1105", "PHL1106", "PHL1107", "PHL1108", "PHL1109", "PHL1110"))))

overlaps_primary <- overlaps_primary %>%
  rename(Area_ha_list = Area_ha) %>%
  rename(Area_ha = Area_ha_R)

#####
####### FUNCTIONS
#####
## This is how we add the primary parameters together 

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


## This is how the OVERLAPS are summarized 

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

#####
##### Tidy and define functions
#####

# create csv for aggregating all sites
# multiply value by (# of overlaps - 1) to account for triple/quadruple/etc counting 
# and make values negative to easily sum
# note we only use this for variables we know don't vary within overlapping sites - for example,
# all overlapping sites are within the same country, but they are not all within the same sls

## To correct for the overlaps,
## We need to make them "negative". 

## Correct to make negative
overlaps_primary_corrected <- overlaps_primary %>% 
  mutate(across(
    .cols = c(Area_ha, rest_area, population, 
              tstor_woody, tstor_soil, tstor_total, 
              carbon_seq_potl),
    ~ as.numeric(.x) * - (n.overlaps - 1))) %>% 
  dplyr::select(!n.overlaps)

#####
##### Country-level
#####

## ALL SITES
country_sites_primary <- sites_primary %>% 
  group_by(Country) %>% 
  primary_summarize() %>% 
  ungroup() 

## OVERLAPS
country_overlaps_primary_wo_ae <- overlaps_primary_corrected %>%
  rowwise() %>% 
  mutate(Country = (unique(Country))) %>% # all countries within polygons are the same
  ungroup() %>% 
  group_by(Country) %>% 
  primary_overlap_summarize() %>% 
  ungroup() 

# bring in emissions avoided to overlap
# break up ae by % area

# this is using the AEs to ... calculate 
# how much avoided emission happens per country
total_aes <- country_sites_primary %>% 
  dplyr::select(Country, Area_ha, emissions_avoided_MgCO2e) %>% 
  rename(total_area = Area_ha)

# Then to calculate the AEs in the overlaps, it's 
# still just grouping essentially the overlaps within a country
# and the percentage of... how much overlap happens
# multiplied by how much avoided emissions happens
# so if we grouped by Country AND whether or not emissions were
# added or subtracted.. and the area in which it was added or subtracted
# We.... well, you can't do that because you'd need to know
# by ID -- how much overlap is there per ID?
# so basically i'd have to calculate... the overlap. I don't think it's actually
# really possible. 

ae_ref_country <- country_overlaps_primary_wo_ae %>%  ## Take the 
  dplyr::select(Country, Area_ha) %>% 
  left_join(total_aes, by = "Country") %>% 
  mutate(pct_area = Area_ha/total_area) %>% 
  mutate(emissions_avoided_MgCO2e = emissions_avoided_MgCO2e * pct_area) %>% 
  dplyr::select(Country, emissions_avoided_MgCO2e)

country_overlaps_primary <- country_overlaps_primary_wo_ae %>% 
  left_join(ae_ref_country, by = 'Country')

country_corrected_primary <- country_sites_primary %>% 
  bind_rows(country_overlaps_primary) %>% 
  group_by(Country) %>% 
  primary_summarize() %>%  # can sum bc overlap values are negative
  ungroup() 

write_csv(country_corrected_primary, 
          paste0("results/summaries_", year, "/ImpactIndicators_CountrySummary_PrimaryIndicators.csv"))


## COUNTRY + AVOIDED EMISSIONS BY COUNTRY

# primary indicators

country_sites_primary_avoid <- sites_primary %>%
  filter(!Intervention_Type %in% c("Extractives - Mining", "REDD+")) %>%
  filter(!Biome %in% c("Terrestrial + Marine", "Marine")) %>%
  group_by(Country, AvoidedOrAdded) %>% 
  primary_summarize() %>% 
  ungroup() 

country_sites_primary_avoid_land <- sites_primary %>%
  filter(!Intervention_Type %in% c("Extractives - Mining", "REDD+")) %>%
  filter(!Biome %in% c("Terrestrial + Marine", "Marine")) %>%
  group_by(Country, AvoidedOrAdded) %>% 
  primary_summarize() %>% 
  ungroup() 

overlaps_primary_corrected <- overlaps_primary_corrected %>%
  filter(!map_lgl(Intervention_Type, ~any(.x %in% c("Extractives - Mining", "REDD+"))))


overlaps_primary_corrected_land <- overlaps_primary_corrected %>%
  filter(!map_lgl(Intervention_Type, ~any(.x %in% c("Extractives - Mining", "REDD+")))) %>%
  filter(!map_lgl(Biome, ~any(.x %in% c("Terrestrial + Marine", "Marine")))) 

country_overlaps_primary_wo_ae_avoid_land <- overlaps_primary_corrected_land %>%
  rowwise() %>% 
  mutate(Country = (unique(Country))) %>% # all countries within polygons are the same
  ungroup() %>% 
  group_by(Country) %>% 
  primary_overlap_summarize() %>% 
  ungroup() 

# bring in emissions avoided to overlap
# break up ae by % area
total_aes <- country_sites_primary_avoid_land %>% 
  dplyr::select(Country, Area_ha, emissions_avoided_MgCO2e, AvoidedOrAdded) %>% 
  rename(total_area = Area_ha)

ae_ref <- country_overlaps_primary_wo_ae_avoid_land %>% 
  dplyr::select(Country, Area_ha) %>% 
  left_join(total_aes, by = "Country") %>% 
  mutate(pct_area = Area_ha/total_area) %>% 
  mutate(emissions_avoided_MgCO2e = emissions_avoided_MgCO2e * pct_area) %>% 
  dplyr::select(Country, emissions_avoided_MgCO2e)

ae_ref$AvoidedOrAdded <- ifelse(ae_ref$emissions_avoided_MgCO2e >= 0, "Added", "Avoided")


country_overlaps_primary_avoid_land <- country_overlaps_primary_wo_ae_avoid_land %>% 
  left_join(ae_ref, by = 'Country')

# Each row in `x` is expected to match at most 1 row in `y`.
# Row 1 of `x` matches multiple rows.
# If multiple matches are expected, set `multiple = "all"` to silence this warnin

country_corrected_primary_avoid_land <- country_sites_primary_avoid_land %>% 
  bind_rows(country_overlaps_primary_avoid_land) %>% 
  group_by(Country, AvoidedOrAdded) %>% 
  primary_summarize() %>%  # can sum bc overlap values are negative
  ungroup() 

country_corrected_primary_bycountry_land <- country_sites_primary_avoid_land %>% 
  bind_rows(country_overlaps_primary_avoid_land) %>% 
  group_by(Country) %>% 
  primary_summarize() %>%  # can sum bc overlap values are negative
  ungroup() 

write_csv(country_corrected_primary_avoid_land, 
          paste0("results/summaries_", year, "/reports/FY24_CountryEmissions_PrimaryIndicators_AvoidAdd_TerrestrialOnly.csv"))

write_csv(country_corrected_primary_bycountry, 
          paste0("results/summaries_", year, "/ImpactIndicators_CountryEmissions_PrimaryIndicators.csv"))



#####
##### Site Status + Country
#####


## New

new_sites_primary <- sites_primary %>% dplyr::filter(Site_Status == "New this year")

overlaps_corrected_new <- overlaps_primary_corrected %>%
  filter(map_lgl(Site_Status, ~ any(grepl("New this year", .))))

## OVERLAPS
new_overlaps_primary <- overlaps_corrected_new %>%
  rowwise() %>% 
  mutate(Country = unique(Country)) %>% # all countries within polygons are the same
  ungroup() %>% 
  group_by(Country) %>% 
  primary_summarize() %>% 
  ungroup()

country_corrected_primary_new <- new_sites_primary %>% 
  bind_rows(new_overlaps_primary) %>% 
  group_by(Country) %>% 
  primary_summarize() %>% 
  ungroup()

write_csv(country_corrected_primary_new, 
          paste0("results/summaries_", year, "/ImpactIndicators_CountrySummary_New_Primary.csv"))

### Continued
contd_sites_primary <- sites_primary %>% dplyr::filter(Site_Status == "Continued this year")

overlaps_corrected_contd <- overlaps_primary_corrected %>%
  filter(map_lgl(Site_Status, ~ any(grepl("Continued this year", .))))

## OVERLAPS
contd_overlaps_primary <- overlaps_corrected_contd %>%
  rowwise() %>% 
  mutate(Country = unique(Country)) %>% # all countries within polygons are the same
  ungroup() %>% 
  group_by(Country) %>% 
  primary_overlap_summarize() %>% 
  ungroup()

country_corrected_primary_contd <- contd_sites_primary %>% 
  bind_rows(contd_overlaps_primary) %>% 
  group_by(Country) %>%
  primary_summarize() %>%
  ungroup()


write_csv(country_corrected_primary_contd, 
          paste0("results/summaries_", year, "/ImpactIndicators_CountrySummary_Contd_Primary.csv"))


#####
##### Division-level
#####

# primary indicators

division_sites_primary <- sites_primary %>% 
  group_by(Division_Primary_CI) %>% 
  primary_summarize() %>% 
  ungroup()

all_identical <- function(x) {
  length(unique(x)) == 1
}

overlaps_primary_corrected_division <- overlaps_primary_corrected %>%
  filter(map_lgl(Division_Primary_CI, all_identical))


division_overlaps_primary_wo_ae <- overlaps_primary_corrected_division %>%
  rowwise() %>% 
  mutate(Division_Primary_CI = unique(Division_Primary_CI)) %>% 
  ungroup() %>% 
  group_by(Division_Primary_CI) %>% 
  primary_summarize() %>% 
  ungroup()

# bring in emissions avoided to overlap
# break up ae by % area
total_aes <- division_sites_primary %>% 
  dplyr::select(Division_Primary_CI, Area_ha, emissions_avoided_MgCO2e) %>% 
  rename(total_area = Area_ha)

ae_ref <- division_overlaps_primary_wo_ae %>% 
  dplyr::select(Division_Primary_CI, Area_ha) %>% 
  left_join(total_aes, by = "Division_Primary_CI") %>% 
  mutate(pct_area = Area_ha/total_area) %>% 
  mutate(emissions_avoided_MgCO2e = emissions_avoided_MgCO2e * pct_area) %>% 
  dplyr::select(Division_Primary_CI, emissions_avoided_MgCO2e)

division_overlaps_primary <- division_overlaps_primary_wo_ae %>% 
  left_join(ae_ref, by = 'Division_Primary_CI')

division_corrected_primary <- division_sites_primary %>% 
  bind_rows(division_overlaps_primary_wo_ae) %>% 
  group_by(Division_Primary_CI) %>% 
  primary_summarize() %>% 
  ungroup()

write_csv(division_corrected_primary, 
          paste0("results/summaries_", year, "/ImpactIndicators_DivisionSummary_PrimaryIndicators.csv"))

## BIOME

# New-Old primarys By Biome
new_old_sites_primary <- sites_primary %>% 
  group_by(Site_Status, Biome) %>% 
  primary_summarize() %>% 
  ungroup()

# There are unfortunately a LOT of sites overlapping with different Biomes...
# print(unique(new_old_overlaps_primary_wo_ae$duplicate_biome))

## DO THIS PART TO FIND THE "CONTINUED + NEW" overlap
new_old_overlaps_primary_wo_ae <- overlaps_primary %>% 
  rowwise() %>% 
  # get a list of duplicated values
  mutate(duplicate_new_cont = list(Site_Status[duplicated(Site_Status)])) %>% 
  ungroup() %>% 
  # remove those without any duplicates 
  # (i.e. birds head & west papua have some overlaps that aren't double counting when we group by sls)
  filter(!duplicate_new_cont == "character(0)") %>%
  rowwise() %>%
  mutate(across(
    .cols = c(Area_ha, rest_area, population,
              tstor_woody, tstor_soil, tstor_total,
              carbon_seq_potl),
    ~ .x * -length(duplicate_new_cont)))

print(unique(new_old_overlaps_primary_wo_ae$duplicate_new_cont))

## create df to filter out 
continuednew_primary <- new_old_overlaps_primary_wo_ae %>%
  filter("Continued this year" %in% duplicate_new_cont & "New this year" %in% duplicate_new_cont)


new_old_overlaps_primary_wo_ae_biome <- overlaps_primary %>% 
  rowwise() %>% 
  # get a list of duplicated values
  mutate(duplicate_biome = list(Biome[duplicated(Biome)])) %>%
  ungroup() %>%
  # remove those without any duplicates
  # (i.e. birds head & west papua have some overlaps that aren't double counting when we group by sls)
  filter(!duplicate_biome == "character(0)") %>%
  mutate(across(
    .cols = c(Area_ha, rest_area, population,
              tstor_woody, tstor_soil, tstor_total,
              carbon_seq_potl),
    ~ .x * -length(duplicate_biome)))

print(unique(new_old_overlaps_primary_wo_ae_biome$duplicate_biome))




## Calculate actual overlaps

new_old_overlaps_primary_wo_ae <- overlaps_primary %>% 
  rowwise() %>% 
  # get a list of duplicated values
  mutate(duplicate_new_cont = list(Site_Status[duplicated(Site_Status)])) %>% 
  ungroup() %>% 
  # remove those without any duplicates 
  # (i.e. birds head & west papua have some overlaps that aren't double counting when we group by sls)
  filter(!duplicate_new_cont == "character(0)") %>% 
  mutate(duplicate_biome = list(Biome[duplicated(Biome)])) %>%
  ungroup() %>%
  # remove those without any duplicates
  # (i.e. birds head & west papua have some overlaps that aren't double counting when we group by sls)
  filter(!duplicate_biome == "character(0)") %>%
  ungroup() %>%
  rowwise() %>%
  mutate(across(
    .cols = c(Area_ha, rest_area, population,
              tstor_woody, tstor_soil, tstor_total,
              carbon_seq_potl),
    ~ .x * -length(duplicate_new_cont))) %>%
  anti_join(continuednew_primary) %>%
  mutate(Site_Status = unique(duplicate_new_cont[])) %>%
  ungroup() %>%
  #mutate(Biome = unique(Biome[1])) %>%
  #ungroup() %>%
  group_by(Site_Status, Biome) %>%
  primary_summarize() %>%
  ungroup()

# bring in emissions avoided to overlap
# break up ae by % area
total_aes <- new_old_sites_primary %>% 
  dplyr::select(Site_Status, Area_ha, emissions_avoided_MgCO2e, Biome) %>% 
  rename(total_area = Area_ha)

ae_ref <- new_old_overlaps_primary_wo_ae %>% 
  dplyr::select(Biome, Area_ha) %>% 
  left_join(total_aes, by = "Biome") %>% 
  mutate(pct_area = Area_ha/total_area) %>% 
  mutate(emissions_avoided_MgCO2e = emissions_avoided_MgCO2e * pct_area) %>% 
  dplyr::select(Biome, emissions_avoided_MgCO2e)

biome_overlaps_primary <- sls_overlaps_primary_wo_ae %>% 
  left_join(ae_ref, by = 'Biome')


#####
##### Restoration + Country
#####


## Site Level - Restoration
sites_primary_restor <- sites_primary %>%
  filter(Intervention_Category == "Restoration Areas")

write_csv(sites_primary_restor, 
          paste0("results/summaries_", year, "/FY", str_sub(year, start = 3), "_ImpactIndicators_Sites_Restoration.csv"))

overlaps_primary_corrected_restor <- overlaps_primary_corrected %>%
  filter(sapply(Intervention_Category, function(x) all(x == "Restoration Areas"))) %>%
  filter(sapply(Intervention_Type, function(x) any(duplicated(x))))

## BY RESTORATION TYPE:

rest_type_sites_primary <- sites_primary_restor %>% 
  group_by(Intervention_Type) %>% 
  primary_summarize() %>% 
  ungroup() 

rest_type_overlaps_primary_wo_ae <- overlaps_primary_corrected_restor %>%
  rowwise() %>% 
  mutate(Intervention_Type = (unique(Intervention_Type))) %>% # all countries within polygons are the same
  ungroup() %>% 
  group_by(Intervention_Type) %>% 
  primary_overlap_summarize() %>% 
  ungroup() 

rest_type_corrected_primary <- rest_type_sites_primary %>% 
  bind_rows(rest_type_overlaps_primary_wo_ae) %>% 
  group_by(Intervention_Type) %>% 
  primary_summarize() %>%  # can sum bc overlap values are negative
  ungroup() 

write_csv(rest_type_corrected_primary, 
          paste0("results/summaries_", year, "/ImpactIndicators_RestTypeSummary_Primary.csv"))

## BY COUNTRY

rest_sites_country_primary <- sites_primary_restor %>% 
  group_by(Country) %>% 
  primary_summarize() %>% 
  ungroup() 

rest_overlaps_country_primary_wo_ae <- overlaps_primary_corrected_restor %>%
  rowwise() %>% 
  mutate(Country = (unique(Country))) %>% # all countries within polygons are the same
  ungroup() %>% 
  group_by(Country) %>% 
  primary_overlap_summarize() %>% 
  ungroup() 

rest_by_country_corrected_primary <- rest_sites_country_primary %>% 
  bind_rows(rest_overlaps_country_primary_wo_ae) %>% 
  group_by(Country) %>% 
  primary_summarize() %>%  # can sum bc overlap values are negative
  ungroup() 

write_csv(rest_by_country_corrected_primary, 
          paste0("results/summaries_", year, "/ImpactIndicators_Rest_ByCountry.csv"))


## BY COUNTRY & REST TYPE

rest_type_country_primary <- sites_primary_restor %>% 
  group_by(Country, Intervention_Type) %>% 
  primary_summarize() %>% 
  ungroup() 

rest_type_country_overlaps_wo_ae <- overlaps_primary_corrected_restor %>%
  rowwise() %>% 
  mutate(Country = (unique(Country))) %>% # all countries within polygons are the same
  mutate(Intervention_Type = (unique(Intervention_Type))) %>%
  ungroup() %>% 
  group_by(Country,Intervention_Type) %>% 
  primary_overlap_summarize() %>% 
  ungroup() 

rest_type_country_corrected_primary <- rest_type_country_primary %>% 
  bind_rows(rest_type_country_overlaps_wo_ae) %>% 
  group_by(Country, Intervention_Type) %>% 
  primary_summarize() %>%  # can sum bc overlap values are negative
  ungroup() 

write_csv(rest_type_country_corrected_primary, 
          paste0("results/summaries_", year, "/ImpactIndicators_Rest_Type_ByCountry.csv"))



#####
##### Mastercard PPC
#####

## Site Level - Restoration
sites_PPC_restor <- sites_primary %>%
  filter(Intervention_Category == "Restoration Areas") %>%
  filter(CI_Portfolio == "Mastercard PPC")

overlaps_corrected_PPC_restor <- overlaps_primary_corrected %>%
  filter(sapply(Intervention_Category, function(x) all(x == "Restoration Areas"))) %>%
  filter(sapply(Intervention_Type, function(x) any(duplicated(x)))) %>%
  filter(sapply(CI_Portfolio, function(x) all(x == "Mastercard PPC")))


## BY COUNTRY & REST TYPE

rest_type_country_primary_PPC <- sites_PPC_restor %>% 
  group_by(Country, Intervention_Type) %>% 
  primary_summarize() %>% 
  ungroup() 

rest_type_country_overlaps_PPC <- overlaps_corrected_PPC_restor %>%
  rowwise() %>% 
  mutate(Country = (unique(Country))) %>% # all countries within polygons are the same
  mutate(Intervention_Type = (unique(Intervention_Type))) %>%
  ungroup() %>% 
  group_by(Country,Intervention_Type) %>% 
  primary_overlap_summarize() %>% 
  ungroup() 

rest_type_country_corrected_primary_PPC <- rest_type_country_primary_PPC %>% 
  bind_rows(rest_type_country_overlaps_PPC) %>% 
  group_by(Country, Intervention_Type) %>% 
  primary_summarize() %>%  # can sum bc overlap values are negative
  ungroup() 

write_csv(rest_type_country_corrected_primary_PPC, 
          paste0("results/summaries_", year, "/ImpactIndicators_Rest_Type_ByCountry_PPC.csv"))


#####
##### Site-level
#####

# overlaps don't occur at site level
# due to that, note can't sum for total (or will have double counting)

write_csv(sites_primary, 
          paste0("results/summaries_", year, "/ImpactIndicators_SiteSummary_PrimaryIndicators.csv"))


