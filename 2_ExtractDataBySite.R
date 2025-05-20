# Impact Indicators Analysis
# Contact: Anna Ballasiotes
# Last Updated: 10/25/2024

library(sf)
library(terra)
library(exactextractr)
library(janitor)
library(googledrive)
library(tidyverse)
library(purrr)

## There are two streams of processes here:
## One for "Other" (which I renamed to Primary late in FY24) variables
## One for "IC" variables (Irrecoverable Carbon)

## Main Outputs:
### Overlaps must be stored as RDS because some columns contain lists

# FY24_ImpactIndicators_Primary_Sites_ALL.csv
# FY24_ImpactIndicators_Primary_Overlaps.rds
# FY24_ImpactIndicators_IrrecoverableCarbon_Sites.csv
# FY24_ImpactIndicators_IrrecoverableCarbon_Overlaps.rds


year <- "2024"

#####
##### Read in data
#####

### irrecoverable carbon with layers masked to each ecosystem
ic_total <- rast("data/irrecoverable_carbon/tstor_ic_ecosystem_stack.tif")
ic_high <- rast("data/irrecoverable_carbon/ha_high_ic_ecosystem_stack.tif") # ha w/ ic over 25
ic_any <- rast("data/irrecoverable_carbon/ha_ic_ecosystem_stack.tif") # ha w/ ic over 0.01
ic_stack <- c(ic_total, ic_high, ic_any)

ic_stack <- rast("data/irrecoverable_carbon/ic_stack_prepped.tif")

### define crs 
my_crs <- crs(ic_stack)


# woody and soil carbon
carbon_woody <- rast(paste0("data/carbon_stored/biomass_prepped_2024", ".tif"))
# agb + bgb w/ hansen forest loss to 2022 masked out aggregated to ~300m
carbon_soil <- rast("data/carbon_stored/OCSTHA_30cm_1km.tif")


### population
pop <- rast("data/population/ppp_2020_1km_Aggregated.tif") 


### carbon sequestration potential
seq <- rast("data/carbon_sequestration/carbon_sequestration_potl_stack.tif")


### avoided emissions - by site, using Alex's ae code
ae <- readRDS("data/avoided_emissions/ae_by_site_by_year_2024.rds") %>% 
  dplyr::select(CI_ID, emissions_avoided_MgCO2e_2023)

ae$AvoidedOrAdded <- ifelse(ae$emissions_avoided_MgCO2e_2023 >= 0, "Added", "Avoided")


#####

##### Prepare shapefiles and intersections
#####

# ignore errors around spherical geometries
sf::sf_use_s2(FALSE)

# calculated project length as beginning > end, and if no end then used today [[need to verify]]
adjusted_end_date <- as.Date(paste0(as.numeric(year) + 1, "-12-31"))

######

FY2024_Sites <- readRDS(paste0("data/ci_sites/FY", year, "_Sites_Clean.rds"))%>% 
  rowwise()

# Correcting date type, if needed

FY2024_Sites$CI_End_Date <- as.Date(FY2024_Sites$CI_End_Date , format = "%Y/%m/%d")
FY2024_Sites$CI_Start_Date <- as.Date(FY2024_Sites$CI_Start_Date , format = "%Y/%m/%d")
FY2024_Sites$Gazettement_Date <- as.Date(FY2024_Sites$Gazettement_Date, format = "%Y/%m/%d")
FY2024_Sites$Tree_Planting_Date <- as.Date(FY2024_Sites$Tree_Planting_Date, format = "%Y/%m/%d")


FY24_sites <- FY2024_Sites %>%
  mutate(project_length = lubridate::time_length(
    difftime(if_else(is.na(CI_End_Date), 
                     adjusted_end_date, 
                     CI_End_Date), 
             CI_Start_Date), "years"))

# remove geometries for faster df 
FY24_sites_df <- FY24_sites %>% st_drop_geometry()

# df of intersections among site polygons
# same changes as for shp
int <- readRDS(paste0("data/ci_sites/FY", year, "_Overlaps_Clean.rds")) %>% 
 dplyr::select("n.overlaps", "origins", "Country", "Area_ha_R", "Division_Primary_CI", "Division_Secondary_CI",
               "CI_SLS", "CI_Role", "Intervention_Category", "Geographic_Priority", "Initial_Reporting_Year",  "CI_ID_1",  "CI_ID_2",
                "CI_ID_3",  "CI_ID_4",  "CI_ID_5", "CI_ID_6")

int_df <- int %>% st_drop_geometry()  

write_csv(int_df, "data/ci_sites/int_df_2024.csv")


#####

##### Extract data to polygons
#####

### exact extract

# population
pop_extract_shp <- exact_extract(
  x = pop,
  y = FY24_sites,
  fun = "sum",
  append_cols = colnames(FY24_sites)) %>% 
  rename("population" = sum)


pop_extract_int <- exact_extract(
  x = pop,
  y = int,
  fun = "sum",
  append_cols = colnames(int)) %>% 
  rename("population" = sum)

# carbon stored - soil
soil_extract_shp <- exact_extract(
  x = carbon_soil,
  y = FY24_sites,
  fun = "sum",
  append_cols = "CI_ID") %>% 
  rename("tstor_soil" = sum)

soil_extract_int <- exact_extract(
  x = carbon_soil,
  y = int,
  fun = "sum",
  append_cols = c("origins")) %>% 
  rename("tstor_soil" = sum)

# carbon stored - woody
woody_extract_shp <- exact_extract(
  x = carbon_woody,
  y = FY24_sites,
  fun = "sum",
  append_cols = "CI_ID") %>% 
  rename("tstor_woody" = sum)

woody_extract_int <- exact_extract(
  x = carbon_woody,
  y = int,
  fun = "sum",
  append_cols = c("origins")) %>% 
  rename("tstor_woody" = sum) 

# carbon sequestration potential
# carbon seq potential depends on restoration type
## assuming mixed can be 'other broadleaf' - no metric for mixed
## for future years may need to add more in here - there are more layers (see names(seq))
## but this is all the restoration activities present in the df this script was made with

seq_extract_shp <- exact_extract(
  x = seq,
  y = FY24_sites,
  fun = "mean",
  append_cols = c("CI_ID", "Area_ha", "Intervention_Type", "CI_End_Date", "CI_Start_Date", "project_length")) %>%  
  rowwise() %>%
  mutate(carbon_seq_potl = case_when(
    Intervention_Type == "Restoration - Agroforestry" & project_length < 20 ~ mean.agfor0020,
    Intervention_Type == "Restoration - Agroforestry" & project_length >= 20 ~ mean.agfor2060,
    Intervention_Type == "Restoration - Mangrove Shrub" & project_length < 20 ~ mean.mshrr0020,
    Intervention_Type == "Restoration - Mangrove Shrub" & project_length >= 20 ~ mean.mshrr2060,
    Intervention_Type == "Restoration - Mangrove Tree" & project_length < 20 ~ mean.mtrer0020,
    Intervention_Type == "Restoration - Mangrove Tree" & project_length >= 20 ~ mean.mtrer2060,
    Intervention_Type == "Restoration - Natural Regeneration" & project_length < 20 ~ mean.natre0020,
    Intervention_Type == "Restoration - Natural Regeneration" & project_length >= 20 ~ mean.natre2060,
    Intervention_Type == "Restoration - Plantations & Woodlots (Eucalyptus)" ~ mean.pweuc0020,
    Intervention_Type == "Restoration - Plantations & Woodlots (Mixed 50/50)" ~ mean.pwobr0020, 
    Intervention_Type == "Restoration - Enrichment Planting/Assisted Natural Regeneration" ~ mean.natre0020,
    Intervention_Type == "Rangeland Restoration - Planned Grazing" ~ 3.67,
    Intervention_Type == "Restoration - Seed Dispersal" ~ mean.natre0020 * 0.6,
    Intervention_Type == "Restoration - Silvopasture" & project_length < 20 ~ mean.agfor0020 * 0.15,
    Intervention_Type == "Restoration - Silvopasture" & project_length >= 20 ~ mean.agfor2060 * 0.15,
    Intervention_Type == "Restoration - Plantations & Woodlots  (Other Broadleaf)" & project_length < 20 ~ mean.pwobr0020,
    Intervention_Type == "Restoration - Plantations & Woodlots (Pine)" & project_length < 20 ~ mean.pwpin0020,
    Intervention_Type == "Restoration - Plantations & Woodlots (Other Conifer)" & project_length < 20 ~ mean.pwoco0020
      )) %>%
  mutate(carbon_seq_potl = carbon_seq_potl * Area_ha) %>% 
  ungroup() %>% 
  dplyr::select(CI_ID, carbon_seq_potl)

saveRDS(seq_extract_shp, 
        paste0("data/carbon_sequestration/FY", str_sub(year, start = 3),
               "_SequestrationExtract_ALL.rds"))


# assumption here - for those that overlap with different restoration activities, will
# calculate carbon sequestration as min (w/ na.rm = FALSE) between the two to avoid any 
# negative numbers when accounting for double counting (this was an issue)

# join restoration activity + length info to ints
rest_info <- FY24_sites_df %>% 
  dplyr::select(CI_ID, Intervention_Type, Intervention_Category, project_length, Area_ha) %>% 
  mutate(under_restoration = (Intervention_Category == 'Restoration Areas'))

# checked and rest activities within 3rd and 4th overlaps are all captured within first two
# if not the case in future analyses, may have to change to account for 3,4,5,etc

#######################

# combine all except ic - ic will be separate bc divided by ecosystem
shp_all_wo_ae <- pop_extract_shp %>% 
  left_join(woody_extract_shp, by = "CI_ID") %>% 
  left_join(soil_extract_shp, by = "CI_ID") %>% 
  left_join(seq_extract_shp, by = "CI_ID") %>% 
  rowwise() %>% 
  mutate(tstor_total = sum(tstor_woody, tstor_soil)) %>% 
  ungroup() %>% 
  relocate(tstor_total, .before = carbon_seq_potl) 

write_csv(shp_all_wo_ae, 
          paste0("results/FY", str_sub(year, start = 3), 
                 "_ImpactIndicators_Primary_Sites_wo_ae.csv"))

# break up ae by % area
site_areas <- shp_all_wo_ae %>% 
  group_by(CI_ID) %>% 
  summarize("total_area" = sum(Area_ha, na.rm = TRUE)) %>% 
  ungroup()

ae_ref <- shp_all_wo_ae %>% 
  left_join(site_areas, by = 'CI_ID') %>% 
  rowwise() %>% 
  mutate(pct_area = Area_ha/total_area) %>% 
  ungroup() %>% 
  left_join(ae, by = 'CI_ID') %>% ## here we lose all BNA/WWF sites - large portion of emissions (almost all)
  mutate(emissions_avoided_MgCO2e = emissions_avoided_MgCO2e_2023 * pct_area, na.rm = TRUE) %>% 
  mutate(emissions_avoided_MgCO2e_2023 = replace_na(emissions_avoided_MgCO2e, 0)) %>% 
  dplyr::select(CI_ID, emissions_avoided_MgCO2e, AvoidedOrAdded)


shp_all<- shp_all_wo_ae %>% 
  cbind(ae_ref$emissions_avoided_MgCO2e) %>%
  cbind(ae_ref$AvoidedOrAdded) %>%
  rename(emissions_avoided_MgCO2e = 'ae_ref$emissions_avoided_MgCO2e') %>%
  rename(AvoidedOrAdded = 'ae_ref$AvoidedOrAdded')

write_csv(shp_all, 
          paste0("results/FY", str_sub(year, start = 3), 
                 "_ImpactIndicators_Primary_Sites_ALL.csv"))

write_csv(shp_all_wo_ae, 
          paste0("results/FY", str_sub(year, start = 3), 
                 "_ImpactIndicators_Primary_Sites_ALL_wo_ae.csv"))

#######################
#### OVERLAPS #########
#######################

# assumption here - for those that overlap with different restoration activities, will
# calculate carbon sequestration as min (w/ na.rm = FALSE) between the two to avoid any negative numbers when accounting for double counting (this was an issue)

int_rest <- int %>% 
  left_join(rest_info, by = c("CI_ID_1" = "CI_ID")) %>% 
  rename(rest_1 = Intervention_Type,
         length_1 = project_length,
         under_rest_1 = under_restoration) %>%
  mutate("rest_area_1" = case_when(
    under_rest_1 == TRUE ~ Area_ha_R,
    T ~ 0
  )) %>% 
  left_join(rest_info, by = c("CI_ID_2" = "CI_ID")) %>% 
  rename(rest_2 = Intervention_Type,
         length_2 = project_length,
         under_rest_2 = under_restoration) %>%
  mutate("rest_area_2" = case_when(
    under_rest_2 == TRUE ~ Area_ha_R,
    T ~ 0
  ))


# assumption here - for those that overlap with different restoration activities, will
# calculate carbon sequestration as min (w/ na.rm = FALSE) between the two to avoid any negative numbers when accounting for double counting (this was an issue)
# note to self: anna -- this is why you use min!!! changed back to min 04/19/2024

seq_extract_int <- exact_extract(
  x = seq,
  y = int_rest,
  fun = "mean",
  append_cols = colnames(int_rest)
)  %>%
  rowwise() %>%
  mutate(seq_potl_1 = case_when(
    rest_1 == "Restoration - Agroforestry" & length_1 < 20 ~ mean.agfor0020,
    rest_1 == "Restoration - Agroforestry" & length_1 >= 20 ~ mean.agfor2060,
    rest_1 == "Mangrove Shrub Restoration" & length_1 < 20 ~ mean.mshrr0020,
    rest_1 == "Mangrove Shrub Restoration" & length_1 >= 20 ~ mean.mshrr2060,
    rest_1 == "Restoration - Mangrove Tree" & length_1 < 20 ~ mean.mtrer0020,
    rest_1 == "Restoration - Mangrove Tree" & length_1 >= 20 ~ mean.mtrer2060,
    rest_1 == "Restoration - Natural Regeneration" & length_1 < 20 ~ mean.natre0020,
    rest_1 == "Restoration - Natural Regeneration" & length_1 >= 20 ~ mean.natre2060,
    rest_1 == "Restoration - Plantations & Woodlots (Eucalyptus)" ~ mean.pweuc0020,
    rest_1 == "Restoration - Plantations & Woodlots (Mixed 50/50)" ~ mean.pwobr0020, 
    rest_1 == "Restoration - Enrichment Planting/Assisted Natural Regeneration" ~ mean.natre0020,
    rest_1 == "Rangeland Restoration - Planned Grazing" ~ 3.67,
    rest_1 == "Restoration - Seed Dispersal" ~ mean.natre0020 * 0.6,
    rest_1 == "Restoration - Silvopasture" & length_1 < 20 ~ mean.agfor0020 * 0.15,
    rest_1 == "Restoration - Silvopasture" & length_1 >= 20 ~ mean.agfor2060 * 0.15,
    rest_1 == "Restoration - Plantations & Woodlots  (Other Broadleaf)" & length_1 < 20 ~ mean.pwobr0020,
    rest_1 == "Restoration - Plantations & Woodlots (Pine)" & length_1 < 20 ~ mean.pwpin0020,
    rest_1 == "Restoration - Plantations & Woodlots (Other Conifer)" & length_1 < 20 ~ mean.pwoco0020
  )) %>%
  mutate(seq_potl_1 = seq_potl_1) %>% 
  mutate(seq_potl_2 = case_when(
    rest_2 == "Restoration - Agroforestry" & length_2 < 20 ~ mean.agfor0020,
    rest_2 == "Restoration - Agroforestry" & length_2 >= 20 ~ mean.agfor2060,
    rest_2 == "Mangrove Shrub Restoration" & length_2 < 20 ~ mean.mshrr0020,
    rest_2 == "Mangrove Shrub Restoration" & length_2 >= 20 ~ mean.mshrr2060,
    rest_2 == "Restoration - Mangrove Tree" & length_2 < 20 ~ mean.mtrer0020,
    rest_2 == "Restoration - Mangrove Tree" & length_2 >= 20 ~ mean.mtrer2060,
    rest_2 == "Restoration - Natural Regeneration" & length_2 < 20 ~ mean.natre0020,
    rest_2 == "Restoration - Natural Regeneration" & length_2 >= 20 ~ mean.natre2060,
    rest_2 == "Restoration - Plantations & Woodlots (Eucalyptus)" ~ mean.pweuc0020,
    rest_2 == "Restoration - Plantations & Woodlots (Mixed 50/50)" ~ mean.pwobr0020, 
    rest_2 == "Restoration - Enrichment Planting/Assisted Natural Regeneration" ~ mean.natre0020,
    rest_2 == "Rangeland Restoration - Planned Grazing" ~ 3.67,
    rest_2 == "Restoration - Seed Dispersal" ~ mean.natre0020 * 0.6,
    rest_2 == "Restoration - Silvopasture" & length_2 < 20 ~ mean.agfor0020 * 0.15,
    rest_2 == "Restoration - Silvopasture" & length_2 >= 20 ~ mean.agfor2060 * 0.15,
    rest_2 == "Restoration - Plantations & Woodlots  (Other Broadleaf)" & length_2 < 20 ~ mean.pwobr0020,
    rest_2 == "Restoration - Plantations & Woodlots (Pine)" & length_2 < 20 ~ mean.pwpin0020,
    rest_2 == "Restoration - Plantations & Woodlots (Other Conifer)" & length_2 < 20 ~ mean.pwoco0020
  )) %>%
  mutate(seq_potl_2 = seq_potl_2) %>% 
  mutate(carbon_seq_potl = min(seq_potl_1, seq_potl_2, na.rm = FALSE), ## in the original code, this took the min... but why? that doesn't make sense.... 
         rest_area = min(rest_area_1, rest_area_2, na.rm = FALSE)) %>% ## in the original code, this took the min... but why?? that doesn't make sense??? if you're not removing NAs
  mutate(carbon_seq_potl = carbon_seq_potl * Area_ha_R) %>%
  # CSA defined annual rate as 1 ton C/year, scaling to CO2 eq
  ungroup() %>%
  filter(!carbon_seq_potl == -Inf) %>% # dropping bc means rest has NA
  dplyr::select(origins, rest_area, carbon_seq_potl)

int_all <- pop_extract_int %>% 
  left_join(woody_extract_int, by = "origins") %>% 
  left_join(soil_extract_int, by = c("origins")) %>% 
  left_join(seq_extract_int, by = c("origins")) %>% 
  rowwise() %>% 
  mutate(tstor_total = sum(tstor_woody, tstor_soil)) %>% 
  ungroup() %>% 
  relocate(Area_ha_R, .before = population) %>% 
  relocate(rest_area, .before = population) %>% 
  relocate(tstor_total, .before = carbon_seq_potl)


saveRDS(int_all, 
        paste0("data/ci_sites/FY", str_sub(year, start = 3),
               "_Overlaps_WithData.rds"))

fields <- FY24_sites_df %>% 
  dplyr::select(!c(undr_rest, rest_area))

int_all_minus_fields <-  int_all %>% dplyr::select(!c(Country, Division_Primary_CI, Division_Secondary_CI, 
                                                      CI_SLS, CI_Role, Intervention_Category, Geographic_Priority, Initial_Reporting_Year))

# If you need to change some of the field names here
# As not great as it is, I may recommend copy/pasting
# and asking ChatGPT to do it. 

int_w_fields <- int_all_minus_fields %>%
  left_join(fields, by = c("CI_ID_1" = "CI_ID")) %>% 
  rename_with(function(x){paste0(x, "_1")}, .cols = Area_Name:project_length) %>% 
  left_join(fields, by = c("CI_ID_2" = "CI_ID")) %>% 
  rename_with(function(x){paste0(x, "_2")}, .cols = Area_Name:project_length) %>% 
  left_join(fields, by = c("CI_ID_3" = "CI_ID")) %>% 
  rename_with(function(x){paste0(x, "_3")}, .cols = Area_Name:project_length) %>% 
  left_join(fields, by = c("CI_ID_4" = "CI_ID")) %>% 
  rename_with(function(x){paste0(x, "_4")}, .cols = Area_Name:project_length) %>% 
  left_join(fields, by = c("CI_ID_5" = "CI_ID")) %>% 
  rename_with(function(x){paste0(x, "_5")}, .cols = Area_Name:project_length)  %>%
  left_join(fields, by = c("CI_ID_6" = "CI_ID")) %>% 
  rename_with(function(x){paste0(x, "_6")}, .cols = Area_Name:project_length)  %>% 
  rowwise() %>% 
  mutate(
    CI_ID = list(na.omit(c(CI_ID_1, CI_ID_2, CI_ID_3, CI_ID_4, CI_ID_5, CI_ID_6))),
    Biome = list(na.omit(c(Biome_1, Biome_2, Biome_3, Biome_4, Biome_5, Biome_6))),
    Division_Primary_CI = list(na.omit(c(Division_Primary_CI_1, Division_Primary_CI_2, Division_Primary_CI_3, Division_Primary_CI_4, Division_Primary_CI_5, Division_Primary_CI_6))),
    Division_Secondary_CI = list(na.omit(c(Division_Secondary_CI_1, Division_Secondary_CI_2, Division_Secondary_CI_3, Division_Secondary_CI_4, Division_Secondary_CI_5, Division_Secondary_CI_6))),
    Country = list(na.omit(c(Country_1, Country_2, Country_3, Country_4, Country_5, Country_6))),
    CI_Portfolio = list(na.omit(c(CI_Portfolio_1, CI_Portfolio_2, CI_Portfolio_3, CI_Portfolio_4, CI_Portfolio_5, CI_Portfolio_6))),
    Governance_Type = list(na.omit(c(Governance_Type_1, Governance_Type_2, Governance_Type_3, Governance_Type_4, Governance_Type_5, Governance_Type_6))),
    CI_Role = list(na.omit(c(CI_Role_1, CI_Role_2, CI_Role_3, CI_Role_4, CI_Role_5, CI_Role_6))),
    Internal_External = list(na.omit(c(Internal_External_1, Internal_External_2, Internal_External_3, Internal_External_4, Internal_External_5, Internal_External_6))),
    Star_Tag = list(na.omit(c(Star_Tag_1, Star_Tag_2, Star_Tag_3, Star_Tag_4, Star_Tag_5, Star_Tag_6))),
    Geographic_Priority = list(na.omit(c(Geographic_Priority_1, Geographic_Priority_2, Geographic_Priority_3, Geographic_Priority_4, Geographic_Priority_5, Geographic_Priority_6))),
    Area_Name = list(na.omit(c(Area_Name_1, Area_Name_2, Area_Name_3, Area_Name_4, Area_Name_5, Area_Name_6))),
    Intervention_Type = list(na.omit(c(Intervention_Type_1, Intervention_Type_2, Intervention_Type_3, Intervention_Type_4, Intervention_Type_5, Intervention_Type_6))),
    Intervention_Category = list(na.omit(c(Intervention_Category_1, Intervention_Category_2, Intervention_Category_3, Intervention_Category_4, Intervention_Category_5, Intervention_Category_6))),
    Gazettement_Date = list(na.omit(c(Gazettement_Date_1, Gazettement_Date_2, Gazettement_Date_3, Gazettement_Date_4, Gazettement_Date_5, Gazettement_Date_6))),
    IUCN_Category = list(na.omit(c(IUCN_Category_1, IUCN_Category_2, IUCN_Category_3, IUCN_Category_4, IUCN_Category_5, IUCN_Category_6))),
    Tree_Planting_Date = list(na.omit(c(Tree_Planting_Date_1, Tree_Planting_Date_2, Tree_Planting_Date_3, Tree_Planting_Date_4, Tree_Planting_Date_5, Tree_Planting_Date_6))),
    Forest_Age = list(na.omit(c(Forest_Age_1, Forest_Age_2, Forest_Age_3, Forest_Age_4, Forest_Age_5, Forest_Age_6))),
    Comments = list(na.omit(c(Comments_1, Comments_2, Comments_3, Comments_4, Comments_5, Comments_6))),
    Local_Land_Designation = list(na.omit(c(Local_Land_Designation_1, Local_Land_Designation_2, Local_Land_Designation_3, Local_Land_Designation_4, Local_Land_Designation_5, Local_Land_Designation_6))),
    Site_Status_Sequestration = list(na.omit(c(Site_Status_Sequestration_1, Site_Status_Sequestration_2, Site_Status_Sequestration_3, Site_Status_Sequestration_4, Site_Status_Sequestration_5, Site_Status_Sequestration_6))),
    Site_Status = list(na.omit(c(Site_Status_1, Site_Status_2, Site_Status_3, Site_Status_4, Site_Status_5, Site_Status_6))),
    CI_SLS = list(na.omit(c(CI_SLS_1, CI_SLS_2, CI_SLS_3, CI_SLS_4, CI_SLS_5, CI_SLS_6))),
    Initial_Reporting_Year = list(na.omit(c(Initial_Reporting_Year_1, Initial_Reporting_Year_2, Initial_Reporting_Year_3, Initial_Reporting_Year_4, Initial_Reporting_Year_5, Initial_Reporting_Year_6))),
    Improved_Management = list(na.omit(c(Improved_Management_1, Improved_Management_2, Improved_Management_3, Improved_Management_4, Improved_Management_5, Improved_Management_6))),
    CI_Start_Date = list(na.omit(c(CI_Start_Date_1, CI_Start_Date_2, CI_Start_Date_3, CI_Start_Date_4, CI_Start_Date_5, CI_Start_Date_6))),
    CI_End_Date = list(na.omit(c(CI_End_Date_1, CI_End_Date_2, CI_End_Date_3, CI_End_Date_4, CI_End_Date_5, CI_End_Date_6))),
    CI_NPE = list(na.omit(c(CI_NPE_1, CI_NPE_2, CI_NPE_3, CI_NPE_4, CI_NPE_5, CI_NPE_6))),
    #rest_area = list(na.omit(c(rest_area_1, rest_area_2, rest_area_3, rest_area_4, rest_area_5, rest_area_6))),
    project_length = list(na.omit(c(project_length_1, project_length_2, project_length_3, project_length_4, project_length_5, project_length_6))),
    origin = list(na.omit(c(origin_1, origin_2, origin_3, origin_4, origin_5, origin_6))),
    Area_ha = list(na.omit(c(Area_ha_1, Area_ha_2, Area_ha_3, Area_ha_4, Area_ha_5, Area_ha_6)))
  ) %>%
  ungroup() %>%
  dplyr::select(n.overlaps, colnames(fields), Area_ha_R, rest_area, population, tstor_woody, tstor_soil, tstor_total, carbon_seq_potl)


# needs to be rds to maintain lists
saveRDS(int_w_fields, 
        paste0("results/FY", str_sub(year, start = 3),
               "_ImpactIndicators_Other_Overlaps.rds"))

#######################
## Irr carbon ----
#######################

# irrecoverable carbon
ic_extract_shp <- exact_extract(
  x = ic_stack, 
  y = FY24_sites,
  fun = "sum", 
  append_cols = colnames(FY24_sites)
) 

# tidy up and keep relevant, summable values
# note only returns sites with ic

ic_shp_tidy <- ic_extract_shp %>% 
  pivot_longer(cols = starts_with("sum."),
               names_to = "colname",
               values_to = "value") %>% 
  mutate(colname = str_sub(colname, start = 5)) %>%
  mutate(ecosystem = case_when(
    str_detect(colname, "primary_forest") ~ "Primary forest",
    str_detect(colname, "secondary_forest") ~ "Secondary forest",
    str_detect(colname, "grassland") ~ "Grassland",
    str_detect(colname, "wetlands") ~ "Wetlands",
    str_detect(colname, "mangroves") ~ "Mangroves",
    str_detect(colname, "salt_marsh") ~ "Salt marsh",
    str_detect(colname, "seagrass") ~ "Seagrass",
    str_detect(colname, "peatland") ~ "Peatland"
  )) %>%
  drop_na(ecosystem) %>% # because others are totals which can be obtained via summing
  mutate(colname = case_when(
    str_detect(colname, "ha_high_ic") ~ "ha_high_ic",
    str_detect(colname, "ha_ic") ~ "ha_ic",
    T ~ "tstor_ic")) %>%
  pivot_wider(names_from = colname, values_from = value, values_fn = sum) %>%
  rowwise() %>% 
  mutate(tonnes_ha_ic = tstor_ic/Area_ha,
         tstor_blue_ic = case_when(
           ecosystem %in% c("Mangroves", "Salt marsh", "Seagrass") ~ sum(tstor_ic),
           T ~ 0
         )) %>% 
  filter(!tstor_ic == 0)


saveRDS(ic_shp_tidy, paste0("results/FY", str_sub(year, start = 3),
                            "_ImpactIndicators_IrrecoverableCarbon_Sites.rds"))

write_csv(ic_shp_tidy, 
          paste0("results/FY", str_sub(year, start = 3), 
                "_ImpactIndicators_IrrecoverableCarbon_Sites.csv"))

##########
### OVERLAPS
##########

ic_extract_int <- exact_extract(
  x = ic_stack,
  y = int,
  fun = "sum",
  append_cols = colnames(int)
) 

ic_int_tidy <- ic_extract_int %>% 
  dplyr::select(c(n.overlaps, starts_with("CI_ID_"), Area_ha_R, starts_with("sum."))) %>% 
  pivot_longer(cols = starts_with("sum."),
               names_to = "colname",
               values_to = "value") %>% 
  mutate(colname = str_sub(colname, start = 5)) %>%
  mutate(ecosystem = case_when(
    str_detect(colname, "primary_forest") ~ "Primary forest",
    str_detect(colname, "secondary_forest") ~ "Secondary forest",
    str_detect(colname, "grassland") ~ "Grassland",
    str_detect(colname, "wetlands") ~ "Wetlands",
    str_detect(colname, "mangroves") ~ "Mangroves",
    str_detect(colname, "salt_marsh") ~ "Salt marsh",
    str_detect(colname, "seagrass") ~ "Seagrass",
    str_detect(colname, "peatland") ~ "Peatland"
  )) %>%
  drop_na(ecosystem) %>% # because others are totals which can be obtained via summing
  mutate(colname = case_when(
    str_detect(colname, "ha_high_ic") ~ "ha_high_ic",
    str_detect(colname, "ha_ic") ~ "ha_ic",
    T ~ "tstor_ic")) %>% 
  pivot_wider(names_from = colname, values_from = value) %>% 
  rowwise() %>% 
  mutate(tonnes_ha_ic = sum(tstor_ic)/ sum(Area_ha_R),
         tstor_blue_ic = case_when(
           ecosystem %in% c("Mangroves", "Salt marsh", "Seagrass") ~ sum(tstor_ic),
           T ~ 0
         )) %>% 
  filter(!sum(tstor_ic) == 0) %>% 
  dplyr::select(!c(Area_ha_R)) 

fields2 <- shp_all_wo_ae %>% 
  dplyr::select(c(CI_ID:project_length)) %>%
  dplyr::select(!c(undr_rest, rest_area))


# add fields of interest
ic_int_w_fields <- ic_int_tidy %>%
  left_join(fields2, by = c("CI_ID_1" = "CI_ID")) %>% 
  rename_with(function(x){paste0(x, "_1")}, .cols = Area_Name:project_length) %>% 
  left_join(fields2, by = c("CI_ID_2" = "CI_ID")) %>% 
  rename_with(function(x){paste0(x, "_2")}, .cols = Area_Name:project_length) %>% 
  left_join(fields2, by = c("CI_ID_3" = "CI_ID")) %>% 
  rename_with(function(x){paste0(x, "_3")}, .cols = Area_Name:project_length) %>% 
  left_join(fields2, by = c("CI_ID_4" = "CI_ID")) %>% 
  rename_with(function(x){paste0(x, "_4")}, .cols = Area_Name:project_length) %>% 
  left_join(fields2, by = c("CI_ID_5" = "CI_ID")) %>% 
  rename_with(function(x){paste0(x, "_5")}, .cols = Area_Name:project_length)  %>%
  left_join(fields2, by = c("CI_ID_6" = "CI_ID")) %>% 
  rename_with(function(x){paste0(x, "_6")}, .cols = Area_Name:project_length)  %>% 
  rowwise() %>% 
  mutate(
    CI_ID = list(na.omit(c(CI_ID_1, CI_ID_2, CI_ID_3, CI_ID_4, CI_ID_5, CI_ID_6))),
    Biome = list(na.omit(c(Biome_1, Biome_2, Biome_3, Biome_4, Biome_5, Biome_6))),
    Division_Primary_CI = list(na.omit(c(Division_Primary_CI_1, Division_Primary_CI_2, Division_Primary_CI_3, Division_Primary_CI_4, Division_Primary_CI_5, Division_Primary_CI_6))),
    Division_Secondary_CI = list(na.omit(c(Division_Secondary_CI_1, Division_Secondary_CI_2, Division_Secondary_CI_3, Division_Secondary_CI_4, Division_Secondary_CI_5, Division_Secondary_CI_6))),
    Country = list(na.omit(c(Country_1, Country_2, Country_3, Country_4, Country_5, Country_6))),
    CI_Portfolio = list(na.omit(c(CI_Portfolio_1, CI_Portfolio_2, CI_Portfolio_3, CI_Portfolio_4, CI_Portfolio_5, CI_Portfolio_6))),
    Governance_Type = list(na.omit(c(Governance_Type_1, Governance_Type_2, Governance_Type_3, Governance_Type_4, Governance_Type_5, Governance_Type_6))),
    CI_Role = list(na.omit(c(CI_Role_1, CI_Role_2, CI_Role_3, CI_Role_4, CI_Role_5, CI_Role_6))),
    Internal_External = list(na.omit(c(Internal_External_1, Internal_External_2, Internal_External_3, Internal_External_4, Internal_External_5, Internal_External_6))),
    Star_Tag = list(na.omit(c(Star_Tag_1, Star_Tag_2, Star_Tag_3, Star_Tag_4, Star_Tag_5, Star_Tag_6))),
    Geographic_Priority = list(na.omit(c(Geographic_Priority_1, Geographic_Priority_2, Geographic_Priority_3, Geographic_Priority_4, Geographic_Priority_5, Geographic_Priority_6))),
    Area_Name = list(na.omit(c(Area_Name_1, Area_Name_2, Area_Name_3, Area_Name_4, Area_Name_5, Area_Name_6))),
    Intervention_Type = list(na.omit(c(Intervention_Type_1, Intervention_Type_2, Intervention_Type_3, Intervention_Type_4, Intervention_Type_5, Intervention_Type_6))),
    Intervention_Category = list(na.omit(c(Intervention_Category_1, Intervention_Category_2, Intervention_Category_3, Intervention_Category_4, Intervention_Category_5, Intervention_Category_6))),
    Gazettement_Date = list(na.omit(c(Gazettement_Date_1, Gazettement_Date_2, Gazettement_Date_3, Gazettement_Date_4, Gazettement_Date_5, Gazettement_Date_6))),
    IUCN_Category = list(na.omit(c(IUCN_Category_1, IUCN_Category_2, IUCN_Category_3, IUCN_Category_4, IUCN_Category_5, IUCN_Category_6))),
    Tree_Planting_Date = list(na.omit(c(Tree_Planting_Date_1, Tree_Planting_Date_2, Tree_Planting_Date_3, Tree_Planting_Date_4, Tree_Planting_Date_5, Tree_Planting_Date_6))),
    Forest_Age = list(na.omit(c(Forest_Age_1, Forest_Age_2, Forest_Age_3, Forest_Age_4, Forest_Age_5, Forest_Age_6))),
    Comments = list(na.omit(c(Comments_1, Comments_2, Comments_3, Comments_4, Comments_5, Comments_6))),
    Local_Land_Designation = list(na.omit(c(Local_Land_Designation_1, Local_Land_Designation_2, Local_Land_Designation_3, Local_Land_Designation_4, Local_Land_Designation_5, Local_Land_Designation_6))),
    Site_Status_Sequestration = list(na.omit(c(Site_Status_Sequestration_1, Site_Status_Sequestration_2, Site_Status_Sequestration_3, Site_Status_Sequestration_4, Site_Status_Sequestration_5, Site_Status_Sequestration_6))),
    Site_Status = list(na.omit(c(Site_Status_1, Site_Status_2, Site_Status_3, Site_Status_4, Site_Status_5, Site_Status_6))),
    CI_SLS = list(na.omit(c(CI_SLS_1, CI_SLS_2, CI_SLS_3, CI_SLS_4, CI_SLS_5, CI_SLS_6))),
    Initial_Reporting_Year = list(na.omit(c(Initial_Reporting_Year_1, Initial_Reporting_Year_2, Initial_Reporting_Year_3, Initial_Reporting_Year_4, Initial_Reporting_Year_5, Initial_Reporting_Year_6))),
    Improved_Management = list(na.omit(c(Improved_Management_1, Improved_Management_2, Improved_Management_3, Improved_Management_4, Improved_Management_5, Improved_Management_6))),
    CI_Start_Date = list(na.omit(c(CI_Start_Date_1, CI_Start_Date_2, CI_Start_Date_3, CI_Start_Date_4, CI_Start_Date_5, CI_Start_Date_6))),
    CI_End_Date = list(na.omit(c(CI_End_Date_1, CI_End_Date_2, CI_End_Date_3, CI_End_Date_4, CI_End_Date_5, CI_End_Date_6))),
    CI_NPE = list(na.omit(c(CI_NPE_1, CI_NPE_2, CI_NPE_3, CI_NPE_4, CI_NPE_5, CI_NPE_6))),
    #rest_area = list(na.omit(c(rest_area_1, rest_area_2, rest_area_3, rest_area_4, rest_area_5, rest_area_6))),
    project_length = list(na.omit(c(project_length_1, project_length_2, project_length_3, project_length_4, project_length_5, project_length_6))),
    origin = list(na.omit(c(origin_1, origin_2, origin_3, origin_4, origin_5, origin_6))),
    Area_ha = list(na.omit(c(Area_ha_1, Area_ha_2, Area_ha_3, Area_ha_4, Area_ha_5, Area_ha_6)))
  ) %>% 
  ungroup() %>% 
  dplyr::select(n.overlaps, 
                colnames(fields2),
                ecosystem, tstor_ic, ha_high_ic, ha_ic, tonnes_ha_ic, tstor_blue_ic)

saveRDS(ic_int_w_fields, 
        paste0("results/FY", str_sub(year, start = 3),
               "_ImpactIndicators_IrrecoverableCarbon_Overlaps.rds"))
