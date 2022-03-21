# Impact Indicators Analysis
# Contact: Cameryn Brock
# Last Updated: 3/19/2022

library(sf)
library(terra)
library(exactextractr)
library(janitor)
library(tidyverse)



#####

##### Read in data
#####

### irrecoverable carbon with layers masked to each ecosystem
ic_total <- rast("data/irrecoverable_carbon/tstor_ic_ecosystem_stack.tif")
ic_high <- rast("data/irrecoverable_carbon/ha_high_ic_ecosystem_stack.tif") # ha w/ ic over 25
ic_any <- rast("data/irrecoverable_carbon/ha_ic_ecosystem_stack.tif") # ha w/ ic over 0.01

ic_stack <- c(ic_total, ic_high, ic_any)

### define crs 
my_crs <- crs(ic_stack)



# woody and soil carbon
carbon_woody <- rast("data/carbon_stored/biomass_prepped.tif")
# agb + bgb w/ hansen forest loss to 2020 masked out aggregated to ~300m
carbon_soil <- rast("data/carbon_stored/OCSTHA_30cm_1km.tif")



### population
pop <- rast("data/population/ppp_2010_1km_Aggregated.tif") 



### carbon sequestration potential
seq <- rast("data/carbon_sequestration/carbon_sequestration_potl_stack.tif")





#####

##### Prepare shapefiles and intersections
#####

# ignore errors around spherical geometries
sf::sf_use_s2(FALSE)

# site polygons with site-level info
# removed wwf, bna, & slw sites by request
# removed proposed PAs by request
# calculated project length as beggining > end, and if no end then used today [[need to verify]]
adjusted_end_date <- as.Date("2022-03-18")

shp <- read_sf(
  dsn = "data/ci_sites",
  layer = "FY2021_Sites")  %>% 
  rowwise() %>% 
  mutate(project_length = lubridate::time_length(
    difftime(if_else(is.na(ci_end_dat), 
                     adjusted_end_date, 
                     ci_end_dat), 
             ci_start_d), "years")) %>% 
  ungroup()

# remove geometries for faster df 
shp_df <- shp %>% st_drop_geometry()

# df of intersections among site polygons
# same changes as for shp
int <- readRDS("data/ci_sites/FY2021_Overlaps_Clean.rds") 

int_df <- int %>% st_drop_geometry()  


#####

##### Extract data to polygons
#####

### exact extract

# population
pop_extract_shp <- exact_extract(
  x = pop,
  y = shp,
  fun = "sum",
  append_cols = colnames(shp)) %>% 
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
  y = shp,
  fun = "sum",
  append_cols = "ci_id") %>% 
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
  y = shp,
  fun = "sum",
  append_cols = "ci_id") %>% 
  rename("tstor_woody" = sum)

woody_extract_int <- exact_extract(
  x = carbon_woody,
  y = int,
  fun = "sum",
  append_cols = c("origins")) %>% 
  rename("tstor_woody" = sum)



# carbon sequestration potential
# carbon seq potential depends on restoration type and proj length
## assuming mixed can be 'other broadleaf' - no metric for mixed
## for future years may need to add more in here - there are more layers (see names(seq)) but this is all the restoration activities present in the df this script was made with
seq_extract_shp <- exact_extract(
  x = seq,
  y = shp,
  fun = "sum",
  append_cols = c("ci_id", "restoratio", "ci_end_dat", "ci_start_d", "project_length")) %>%
  rowwise() %>%
  mutate(carbon_seq_potl = case_when(
    restoratio == "Agroforestry" & project_length < 20 ~ sum.agfor0020,
    restoratio == "Agroforestry" & project_length >= 20 ~ sum.agfor2060,
    restoratio == "Mangrove Shrub Restoration" & project_length < 20 ~ sum.mshrr0020,
    restoratio == "Mangrove Shrub Restoration" & project_length >= 20 ~ sum.mshrr2060,
    restoratio == "Mangrove Tree Restoration" & project_length < 20 ~ sum.mtrer0020,
    restoratio == "Mangrove Tree Restoration" & project_length >= 20 ~ sum.mtrer2060,
    restoratio == "Natural Regeneration" & project_length < 20 ~ sum.natre0020,
    restoratio == "Natural Regeneration" & project_length >= 20 ~ sum.natre2060,
    restoratio == "Plantations & Woodlots - Eucalyptus" ~ sum.pweuc0020,
    restoratio == "Plantations & Woodlots - Mixed 50/50" ~ sum.pwobr0020, 
    restoratio == "Enrichment Planting/Assisted Natural Regeneration" ~ sum.natre0020,
    restoratio == "Rangeland Restoration - Planned Grazing" ~ sum.agfor0020 * 0.15,
    restoratio == "Seed Dispersal" ~ sum.natre0020 * 0.6
  )) %>%
  mutate(carbon_seq_potl = carbon_seq_potl * project_length * 3.67) %>%  
  # CSA defined annual rate as 1 ton C/year, scaling to CO2 eq
  ungroup() %>% 
  dplyr::select(ci_id, carbon_seq_potl)

# assumption here - for those that overlap with different restoration activities, will
# calculate carbon sequestration as min (w/ na.rm = FALSE) between the two to avoid any negative numbers when accounting for double counting (this was an issue)

# join restoration activity + length info to ints
rest_info <- shp_df %>% 
  dplyr::select(ci_id, restoratio, project_length) %>% 
  mutate(under_restoration = (restoratio != ' ' & 
                                restoratio != 'Not Applicable' & 
                                !is.na(restoratio)))

# checked and rest activities within 3rd and 4th overlaps are all captured within first two
# if not the case in future analyses, may have to change to account for 3,4,5,etc

int_rest <- int %>% 
  left_join(rest_info, by = c("ci_id_1" = "ci_id")) %>% 
  rename(rest_1 = restoratio,
         length_1 = project_length,
         under_rest_1 = under_restoration) %>%
  mutate("rest_area_1" = case_when(
    under_rest_1 == TRUE ~ area_ha,
    T ~ 0
  )) %>% 
  left_join(rest_info, by = c("ci_id_2" = "ci_id")) %>% 
  rename(rest_2 = restoratio,
         length_2 = project_length,
         under_rest_2 = under_restoration) %>%
  mutate("rest_area_2" = case_when(
    under_rest_2 == TRUE ~ area_ha,
    T ~ 0
  ))

seq_extract_int <- exact_extract(
  x = seq,
  y = int_rest,
  fun = "sum",
  append_cols = colnames(int_rest)
)  %>%
  rowwise() %>%
  mutate(seq_potl_1 = case_when(
    rest_1 == "Agroforestry" & length_1 < 20 ~ sum.agfor0020,
    rest_1 == "Agroforestry" & length_1 >= 20 ~ sum.agfor2060,
    rest_1 == "Mangrove Shrub Restoration" & length_1 < 20 ~ sum.mshrr0020,
    rest_1 == "Mangrove Shrub Restoration" & length_1 >= 20 ~ sum.mshrr2060,
    rest_1 == "Mangrove Tree Restoration" & length_1 < 20 ~ sum.mtrer0020,
    rest_1 == "Mangrove Tree Restoration" & length_1 >= 20 ~ sum.mtrer2060,
    rest_1 == "Natural Regeneration" & length_1 < 20 ~ sum.natre0020,
    rest_1 == "Natural Regeneration" & length_1 >= 20 ~ sum.natre2060,
    rest_1 == "Plantations & Woodlots - Eucalyptus" ~ sum.pweuc0020,
    rest_1 == "Plantations & Woodlots - Mixed 50/50" ~ sum.pwobr0020, 
    rest_1 == "Enrichment Planting/Assisted Natural Regeneration" ~ sum.natre0020,
    rest_1 == "Rangeland Restoration - Planned Grazing" ~ sum.agfor0020 * 0.15,
    rest_1 == "Seed Dispersal" ~ sum.natre0020 * 0.6
  )) %>%
  mutate(seq_potl_1 = seq_potl_1 * length_1 * 3.67) %>% 
  mutate(seq_potl_2 = case_when(
    rest_2 == "Agroforestry" & length_2 < 20 ~ sum.agfor0020,
    rest_2 == "Agroforestry" & length_2 >= 20 ~ sum.agfor2060,
    rest_2 == "Mangrove Shrub Restoration" & length_2 < 20 ~ sum.mshrr0020,
    rest_2 == "Mangrove Shrub Restoration" & length_2 >= 20 ~ sum.mshrr2060,
    rest_2 == "Mangrove Tree Restoration" & length_2 < 20 ~ sum.mtrer0020,
    rest_2 == "Mangrove Tree Restoration" & length_2 >= 20 ~ sum.mtrer2060,
    rest_2 == "Natural Regeneration" & length_2 < 20 ~ sum.natre0020,
    rest_2 == "Natural Regeneration" & length_2 >= 20 ~ sum.natre2060,
    rest_2 == "Plantations & Woodlots - Eucalyptus" ~ sum.pweuc0020,
    rest_2 == "Plantations & Woodlots - Mixed 50/50" ~ sum.pwobr0020, 
    rest_2 == "Enrichment Planting/Assisted Natural Regeneration" ~ sum.natre0020,
    rest_2 == "Rangeland Restoration - Planned Grazing" ~ sum.agfor0020 * 0.15,
    rest_2 == "Seed Dispersal" ~ sum.natre0020 * 0.6
  )) %>%
  mutate(seq_potl_2 = seq_potl_2 * length_2 * 3.67) %>% 
  mutate(carbon_seq_potl = min(seq_potl_1, seq_potl_2),
         rest_area = min(rest_area_1, rest_area_2)) %>% 
  # CSA defined annual rate as 1 ton C/year, scaling to CO2 eq
  ungroup() %>% 
  filter(!carbon_seq_potl == Inf) %>% # dropping bc means rest has NA
  dplyr::select(origins, rest_area, carbon_seq_potl)


# combine all except ic - ic will be separate bc divied by ecosystem
shp_all <- pop_extract_shp %>% 
  left_join(woody_extract_shp, by = "ci_id") %>% 
  left_join(soil_extract_shp, by = "ci_id") %>% 
  left_join(seq_extract_shp, by = "ci_id") %>% 
  dplyr::select(c(ci_id, biome, ci_divisio,
                  ci_divis_1, country,
                  new_or_con, new_or_c_1,
                  ci_sls_1, ci_sls_2, 
                  area_ha, rest_area,
                  population, tstor_woody, 
                  tstor_soil, carbon_seq_potl)) %>% 
  rowwise() %>% 
  mutate(tstor_total = sum(tstor_woody, tstor_soil)) %>% 
  ungroup() %>% 
  relocate(tstor_total, .before = carbon_seq_potl)

write_csv(shp_all, "results/FY21_ImpactIndicators_Other_Sites.csv")

int_all <- pop_extract_int %>% 
  left_join(woody_extract_int, by = "origins") %>% 
  left_join(soil_extract_int, by = c("origins")) %>% 
  left_join(seq_extract_int, by = c("origins")) %>% 
  rowwise() %>% 
  mutate(tstor_total = sum(tstor_woody, tstor_soil)) %>% 
  ungroup() %>% 
  relocate(area_ha, .before = population) %>% 
  relocate(rest_area, .before = population) %>% 
  relocate(tstor_total, .before = carbon_seq_potl)

write_csv(int_all, "results/FY21_ImpactIndicators_Other_Overlaps.csv")

# join with field data
# interested in country, division, sls, site
fields <- shp_df %>% 
  dplyr::select(ci_id, country, ci_divisio, ci_divis_1, ci_sls_1, ci_sls_2)

int_w_fields <- int_all %>%
  left_join(fields, by = c("ci_id_1" = "ci_id")) %>%  
  rename(country_1 = country,
         ci_div_1 = ci_divisio,
         ci_div1_1 = ci_divis_1,
         ci_sls1_1 = ci_sls_1,
         ci_sls2_1 = ci_sls_2) %>%
  left_join(fields, by = c("ci_id_2" = "ci_id")) %>% 
  rename(country_2 = country,
         ci_div_2 = ci_divisio,
         ci_div1_2 = ci_divis_1,
         ci_sls1_2 = ci_sls_1,
         ci_sls2_2 = ci_sls_2) %>% 
  left_join(fields, by = c("ci_id_3" = "ci_id")) %>% 
  rename(country_3 = country,
         ci_div_3 = ci_divisio,
         ci_div1_3 = ci_divis_1,
         ci_sls1_3 = ci_sls_1,
         ci_sls2_3 = ci_sls_2) %>% 
  left_join(fields, by = c("ci_id_4" = "ci_id")) %>% 
  rename(country_4 = country,
         ci_div_4 = ci_divisio,
         ci_div1_4 = ci_divis_1,
         ci_sls1_4 = ci_sls_1,
         ci_sls2_4 = ci_sls_2) %>% 
  left_join(fields, by = c("ci_id_5" = "ci_id")) %>% 
  rename(country_5 = country,
         ci_div_5 = ci_divisio,
         ci_div1_5 = ci_divis_1,
         ci_sls1_5 = ci_sls_1,
         ci_sls2_5 = ci_sls_2) %>% 
  rowwise() %>% 
  mutate(
    ci_id = list(na.omit(c(ci_id_1, ci_id_2, ci_id_3, ci_id_4, ci_id_5))),
    country = list(na.omit(c(country_1, country_2, country_3, country_4, country_5))),
    ci_divisio = list(na.omit(c(ci_div_1, ci_div_2, ci_div_3, ci_div_4, ci_div_5))),
    ci_divis_1 = list(na.omit(c(ci_div1_1, ci_div1_2, ci_div1_3, ci_div1_4, ci_div1_5))),
    ci_sls_1 = list(na.omit(c(ci_sls1_1, ci_sls1_2, ci_sls1_3, ci_sls1_4, ci_sls1_5))),
    ci_sls_2 = list(na.omit(c(ci_sls2_1, ci_sls2_2, ci_sls2_3, ci_sls2_4, ci_sls2_5)))) %>% 
  ungroup() %>% 
  dplyr::select(n_overlaps, ci_id, country, ci_divisio, ci_divis_1, ci_sls_1, ci_sls_2,
                area_ha, rest_area,
                population, tstor_woody, tstor_soil, tstor_total, carbon_seq_potl)

# needs to be rds to maintain lists
saveRDS(int_w_fields, "results/FY21_ImpactIndicators_Other_Overlaps.rds")


## Irr carbon

# irrecoverable carbon
ic_extract_shp <- exact_extract(
  x = ic_stack, 
  y = shp,
  fun = "sum", 
  append_cols = colnames(shp)
) 

ic_extract_int <- exact_extract(
  x = ic_stack,
  y = int,
  fun = "sum",
  append_cols = colnames(int)
) 

# tidy up and keep relevant, summable values
# note only returns sites with ic
ic_shp_tidy <- ic_extract_shp %>% 
  dplyr::select(c(ci_id, biome, ci_divisio,
                  ci_divis_1, country,
                  new_or_con, new_or_c_1,
                  ci_sls_1, ci_sls_2, 
                  area_ha, starts_with("sum."))) %>% 
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
  mutate(tonnes_ha_ic = tstor_ic/area_ha,
         tstor_blue_ic = case_when(
           ecosystem %in% c("Mangroves", "Salt marsh", "Seagrass") ~ tstor_ic,
           T ~ 0
         )) %>% 
  filter(!tstor_ic == 0) %>% 
  dplyr::select(!c(area_ha))

write_csv(ic_shp_tidy, "results/FY21_ImpactIndicators_IrrecoverableCarbon_Sites.csv")



ic_int_tidy <- ic_extract_int %>% 
  dplyr::select(c(n_overlaps, starts_with("ci_id_"), area_ha, starts_with("sum."))) %>% 
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
  mutate(tonnes_ha_ic = tstor_ic/area_ha,
         tstor_blue_ic = case_when(
           ecosystem %in% c("Mangroves", "Salt marsh", "Seagrass") ~ tstor_ic,
           T ~ 0
         )) %>% 
  filter(!tstor_ic == 0) %>% 
  dplyr::select(!c(area_ha)) 

write_csv(ic_int_tidy, "results/FY21_ImpactIndicators_IrrecoverableCarbon_Overlaps.csv")

# add fields of interest
ic_int_w_fields <- ic_int_tidy %>%
  left_join(fields, by = c("ci_id_1" = "ci_id")) %>%  
  rename(country_1 = country,
         ci_div_1 = ci_divisio,
         ci_div1_1 = ci_divis_1,
         ci_sls1_1 = ci_sls_1,
         ci_sls2_1 = ci_sls_2) %>%
  left_join(fields, by = c("ci_id_2" = "ci_id")) %>% 
  rename(country_2 = country,
         ci_div_2 = ci_divisio,
         ci_div1_2 = ci_divis_1,
         ci_sls1_2 = ci_sls_1,
         ci_sls2_2 = ci_sls_2) %>% 
  left_join(fields, by = c("ci_id_3" = "ci_id")) %>% 
  rename(country_3 = country,
         ci_div_3 = ci_divisio,
         ci_div1_3 = ci_divis_1,
         ci_sls1_3 = ci_sls_1,
         ci_sls2_3 = ci_sls_2) %>% 
  left_join(fields, by = c("ci_id_4" = "ci_id")) %>% 
  rename(country_4 = country,
         ci_div_4 = ci_divisio,
         ci_div1_4 = ci_divis_1,
         ci_sls1_4 = ci_sls_1,
         ci_sls2_4 = ci_sls_2) %>% 
  left_join(fields, by = c("ci_id_5" = "ci_id")) %>% 
  rename(country_5 = country,
         ci_div_5 = ci_divisio,
         ci_div1_5 = ci_divis_1,
         ci_sls1_5 = ci_sls_1,
         ci_sls2_5 = ci_sls_2) %>% 
  rowwise() %>% 
  mutate(
    ci_id = list(na.omit(c(ci_id_1, ci_id_2, ci_id_3, ci_id_4, ci_id_5))),
    country = list(na.omit(c(country_1, country_2, country_3, country_4, country_5))),
    ci_divisio = list(na.omit(c(ci_div_1, ci_div_2, ci_div_3, ci_div_4, ci_div_5))),
    ci_divis_1 = list(na.omit(c(ci_div1_1, ci_div1_2, ci_div1_3, ci_div1_4, ci_div1_5))),
    ci_sls_1 = list(na.omit(c(ci_sls1_1, ci_sls1_2, ci_sls1_3, ci_sls1_4, ci_sls1_5))),
    ci_sls_2 = list(na.omit(c(ci_sls2_1, ci_sls2_2, ci_sls2_3, ci_sls2_4, ci_sls2_5)))) %>% 
  ungroup() %>% 
  dplyr::select(n_overlaps, ci_id, country, ci_divisio, ci_divis_1, ci_sls_1, ci_sls_2,
                ecosystem, tstor_ic, ha_high_ic, ha_ic, tonnes_ha_ic, tstor_blue_ic)


saveRDS(ic_int_w_fields, "results/FY21_ImpactIndicators_IrrecoverableCarbon_Overlaps.rds")
