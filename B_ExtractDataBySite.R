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
pop <- rast("data/population/ppp_2020_1km_Aggregated.tif") 


### carbon sequestration potential
seq <- rast("data/carbon_sequestration/carbon_sequestration_potl_stack.tif")





#####

##### Prepare shapefiles and intersections
#####

# ignore errors around spherical geometries
sf::sf_use_s2(FALSE)


# join to Howard's intervention categories
inter_codes <- read_csv("misc/intervention_lookup.csv") %>% 
  dplyr::select(!id) %>% 
  rename(interventi = intervention_type,
         inter_cat = intervention_category)

# site polygons with site-level info
# removed wwf, bna, & slw sites by request
# removed proposed PAs by request
# calculated project length as beggining > end, and if no end then used today [[need to verify]]
adjusted_end_date <- as.Date("2022-12-31")

shp <- read_sf(
  dsn = "data/ci_sites",
  layer = "FY2021_Sites")  %>% 
  rowwise() %>% 
  mutate(project_length = lubridate::time_length(
    difftime(if_else(is.na(ci_end_dat), 
                     adjusted_end_date, 
                     ci_end_dat), 
             ci_start_d), "years")) %>% 
  ungroup() %>% 
  left_join(inter_codes, by = "interventi") %>% 
  relocate(inter_cat, .before = interventi)

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
  rename("tstor_woody" = sum) %>% 
  rowwise() %>% 
  mutate(tstor_woody = tstor_woody/2) %>% # adjustment for missing a step in data prep
  ungroup()

woody_extract_int <- exact_extract(
  x = carbon_woody,
  y = int,
  fun = "sum",
  append_cols = c("origins")) %>% 
  rename("tstor_woody" = sum) %>% 
  rowwise() %>% 
  mutate(tstor_woody = tstor_woody/2) %>% 
  ungroup()



# carbon sequestration potential
# carbon seq potential depends on restoration type and proj length
## assuming mixed can be 'other broadleaf' - no metric for mixed
## for future years may need to add more in here - there are more layers (see names(seq)) but this is all the restoration activities present in the df this script was made with
seq_extract_shp <- exact_extract(
  x = seq,
  y = shp,
  fun = "mean",
  append_cols = c("ci_id", "area_ha", "restoratio", "ci_end_dat", "ci_start_d", "project_length")) %>%  
  rowwise() %>%
  mutate(carbon_seq_potl = case_when(
    restoratio == "Agroforestry" & project_length < 20 ~ mean.agfor0020,
    restoratio == "Agroforestry" & project_length >= 20 ~ mean.agfor2060,
    restoratio == "Mangrove Shrub Restoration" & project_length < 20 ~ mean.mshrr0020,
    restoratio == "Mangrove Shrub Restoration" & project_length >= 20 ~ mean.mshrr2060,
    restoratio == "Mangrove Tree Restoration" & project_length < 20 ~ mean.mtrer0020,
    restoratio == "Mangrove Tree Restoration" & project_length >= 20 ~ mean.mtrer2060,
    restoratio == "Natural Regeneration" & project_length < 20 ~ mean.natre0020,
    restoratio == "Natural Regeneration" & project_length >= 20 ~ mean.natre2060,
    restoratio == "Plantations & Woodlots - Eucalyptus" ~ mean.pweuc0020,
    restoratio == "Plantations & Woodlots - Mixed 50/50" ~ mean.pwobr0020, 
    restoratio == "Enrichment Planting/Assisted Natural Regeneration" ~ mean.natre0020,
    restoratio == "Rangeland Restoration - Planned Grazing" ~ 3.67,
    restoratio == "Seed Dispersal" ~ mean.natre0020 * 0.6
  )) %>%
  mutate(carbon_seq_potl = carbon_seq_potl * area_ha) %>% 
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
  fun = "mean",
  append_cols = colnames(int_rest)
)  %>%
  rowwise() %>%
  mutate(seq_potl_1 = case_when(
    rest_1 == "Agroforestry" & length_1 < 20 ~ mean.agfor0020,
    rest_1 == "Agroforestry" & length_1 >= 20 ~ mean.agfor2060,
    rest_1 == "Mangrove Shrub Restoration" & length_1 < 20 ~ mean.mshrr0020,
    rest_1 == "Mangrove Shrub Restoration" & length_1 >= 20 ~ mean.mshrr2060,
    rest_1 == "Mangrove Tree Restoration" & length_1 < 20 ~ mean.mtrer0020,
    rest_1 == "Mangrove Tree Restoration" & length_1 >= 20 ~ mean.mtrer2060,
    rest_1 == "Natural Regeneration" & length_1 < 20 ~ mean.natre0020,
    rest_1 == "Natural Regeneration" & length_1 >= 20 ~ mean.natre2060,
    rest_1 == "Plantations & Woodlots - Eucalyptus" ~ mean.pweuc0020,
    rest_1 == "Plantations & Woodlots - Mixed 50/50" ~ mean.pwobr0020, 
    rest_1 == "Enrichment Planting/Assisted Natural Regeneration" ~ mean.natre0020,
    rest_1 == "Rangeland Restoration - Planned Grazing" ~ 3.67,
    rest_1 == "Seed Dispersal" ~ mean.natre0020 * 0.6
  )) %>%
  mutate(seq_potl_1 = seq_potl_1) %>% 
  mutate(seq_potl_2 = case_when(
    rest_2 == "Agroforestry" & length_2 < 20 ~ mean.agfor0020,
    rest_2 == "Agroforestry" & length_2 >= 20 ~ mean.agfor2060,
    rest_2 == "Mangrove Shrub Restoration" & length_2 < 20 ~ mean.mshrr0020,
    rest_2 == "Mangrove Shrub Restoration" & length_2 >= 20 ~ mean.mshrr2060,
    rest_2 == "Mangrove Tree Restoration" & length_2 < 20 ~ mean.mtrer0020,
    rest_2 == "Mangrove Tree Restoration" & length_2 >= 20 ~ mean.mtrer2060,
    rest_2 == "Natural Regeneration" & length_2 < 20 ~ mean.natre0020,
    rest_2 == "Natural Regeneration" & length_2 >= 20 ~ mean.natre2060,
    rest_2 == "Plantations & Woodlots - Eucalyptus" ~ mean.pweuc0020,
    rest_2 == "Plantations & Woodlots - Mixed 50/50" ~ mean.pwobr0020, 
    rest_2 == "Enrichment Planting/Assisted Natural Regeneration" ~ mean.natre0020,
    rest_2 == "Rangeland Restoration - Planned Grazing" ~ 3.67,
    rest_2 == "Seed Dispersal" ~ mean.natre0020 * 0.6
  )) %>%
  mutate(seq_potl_2 = seq_potl_2) %>% 
  mutate(carbon_seq_potl = min(seq_potl_1, seq_potl_2),
         rest_area = min(rest_area_1, rest_area_2)) %>% 
  mutate(carbon_seq_potl = carbon_seq_potl * area_ha) %>% 
  # CSA defined annual rate as 1 ton C/year, scaling to CO2 eq
  ungroup() %>% 
  filter(!carbon_seq_potl == Inf) %>% # dropping bc means rest has NA
  dplyr::select(origins, rest_area, carbon_seq_potl)


# combine all except ic - ic will be separate bc divied by ecosystem
shp_all <- pop_extract_shp %>% 
  left_join(woody_extract_shp, by = "ci_id") %>% 
  left_join(soil_extract_shp, by = "ci_id") %>% 
  left_join(seq_extract_shp, by = "ci_id") %>% 
  dplyr::select(!c(country_is, star_tag_t, area_name, comment, star_tag_p, star_tag_s, origin, project_length)) %>% 
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


# join with field data
# interested in country, division, sls, site
fields <- shp_df %>% 
  dplyr::select(!c(area_ha, country_is, star_tag_t, area_name, comment, star_tag_p, star_tag_s, rest_area, origin, project_length))

int_w_fields <- int_all %>%
  left_join(fields, by = c("ci_id_1" = "ci_id")) %>% 
  rename_with(function(x){paste0(x, "_1")}, .cols = biome:ci_end_dat) %>% 
  left_join(fields, by = c("ci_id_2" = "ci_id")) %>% 
  rename_with(function(x){paste0(x, "_2")}, .cols = biome:ci_end_dat) %>% 
  left_join(fields, by = c("ci_id_3" = "ci_id")) %>% 
  rename_with(function(x){paste0(x, "_3")}, .cols = biome:ci_end_dat) %>% 
  left_join(fields, by = c("ci_id_4" = "ci_id")) %>% 
  rename_with(function(x){paste0(x, "_4")}, .cols = biome:ci_end_dat) %>% 
  left_join(fields, by = c("ci_id_5" = "ci_id")) %>% 
  rename_with(function(x){paste0(x, "_5")}, .cols = biome:ci_end_dat)  %>% 
  rowwise() %>% 
  mutate(
    ci_id = list(na.omit(c(ci_id_1, ci_id_2, ci_id_3, ci_id_4, ci_id_5))),
    biome = list(na.omit(c(biome_1, biome_2, biome_3, biome_4, biome_5))),
    ci_divisio = list(na.omit(c(ci_divisio_1, ci_divisio_2, ci_divisio_3, ci_divisio_4, ci_divisio_5))),
    ci_divis_1 = list(na.omit(c(ci_divis_1_1, ci_divis_1_2, ci_divis_1_3, ci_divis_1_4, ci_divis_1_5))),
    country = list(na.omit(c(country_1, country_2, country_3, country_4, country_5))),
    geographic = list(na.omit(c(geographic_1, geographic_2, geographic_3, geographic_4, geographic_5))),
    inter_cat = list(na.omit(c(inter_cat_1, inter_cat_2, inter_cat_3, inter_cat_4, inter_cat_5))),
    interventi = list(na.omit(c(interventi_1, interventi_2, interventi_3, interventi_4, interventi_5))),
    interven_1 = list(na.omit(c(interven_1_1, interven_1_2, interven_1_3, interven_1_4, interven_1_5))),
    interven_2 = list(na.omit(c(interven_2_1, interven_2_2, interven_2_3, interven_2_4, interven_2_5))),
    gazettemen = list(na.omit(c(gazettemen_1, gazettemen_2, gazettemen_3, gazettemen_4, gazettemen_5))),
    iucn_categ = list(na.omit(c(iucn_categ_1, iucn_categ_2, iucn_categ_3, iucn_categ_4, iucn_categ_5))),
    tree_plant = list(na.omit(c(tree_plant_1, tree_plant_2, tree_plant_3, tree_plant_4, tree_plant_5))),
    restoratio = list(na.omit(c(restoratio_1, restoratio_2, restoratio_3, restoratio_4, restoratio_5))),
    forest_age = list(na.omit(c(forest_age_1, forest_age_2, forest_age_3, forest_age_4, forest_age_5))),
    restor_fin = list(na.omit(c(restor_fin_1, restor_fin_2, restor_fin_3, restor_fin_4, restor_fin_5))),
    local_land = list(na.omit(c(local_land_1, local_land_2, local_land_3, local_land_4, local_land_5))),
    new_or_con = list(na.omit(c(new_or_con_1, new_or_con_2, new_or_con_3, new_or_con_4, new_or_con_5))),
    new_or_c_1 = list(na.omit(c(new_or_c_1_1, new_or_c_1_2, new_or_c_1_3, new_or_c_1_4, new_or_c_1_5))),  
    ci_sls_1 = list(na.omit(c(ci_sls_1_1, ci_sls_1_2, ci_sls_1_3, ci_sls_1_4, ci_sls_1_5))),
    ci_sls_2 = list(na.omit(c(ci_sls_2_1, ci_sls_2_2, ci_sls_2_3, ci_sls_2_4, ci_sls_2_5))),
    reporting = list(na.omit(c(reporting_1, reporting_2, reporting_3, reporting_4, reporting_5))),
    improved_m = list(na.omit(c(improved_m_1, improved_m_2, improved_m_3, improved_m_4, improved_m_5))),
    ci_start_d = list(na.omit(c(ci_start_d_1, ci_start_d_2, ci_start_d_3, ci_start_d_4, ci_start_d_5))),
    ci_end_dat = list(na.omit(c(ci_end_dat_1, ci_end_dat_2, ci_end_dat_3, ci_end_dat_4, ci_end_dat_5)))
    ) %>% 
  ungroup() %>% 
  dplyr::select(n_overlaps,
                colnames(fields),
                area_ha, rest_area,
                population, tstor_woody, tstor_soil, tstor_total, carbon_seq_potl)

# needs to be rds to maintain lists
saveRDS(int_w_fields, "results/FY21_ImpactIndicators_Other_Overlaps.rds")


## Irr carbon ----

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
  dplyr::select(!c(country_is, star_tag_t, area_name, rest_area, comment,
                   star_tag_p, star_tag_s, origin, project_length)) %>% 
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


# add fields of interest
ic_int_w_fields <- ic_int_tidy %>%
  left_join(fields, by = c("ci_id_1" = "ci_id")) %>% 
  rename_with(function(x){paste0(x, "_1")}, .cols = biome:ci_end_dat) %>% 
  left_join(fields, by = c("ci_id_2" = "ci_id")) %>% 
  rename_with(function(x){paste0(x, "_2")}, .cols = biome:ci_end_dat) %>% 
  left_join(fields, by = c("ci_id_3" = "ci_id")) %>% 
  rename_with(function(x){paste0(x, "_3")}, .cols = biome:ci_end_dat) %>% 
  left_join(fields, by = c("ci_id_4" = "ci_id")) %>% 
  rename_with(function(x){paste0(x, "_4")}, .cols = biome:ci_end_dat) %>% 
  left_join(fields, by = c("ci_id_5" = "ci_id")) %>% 
  rename_with(function(x){paste0(x, "_5")}, .cols = biome:ci_end_dat)  %>% 
  rowwise() %>% 
  mutate(
    ci_id = list(na.omit(c(ci_id_1, ci_id_2, ci_id_3, ci_id_4, ci_id_5))),
    biome = list(na.omit(c(biome_1, biome_2, biome_3, biome_4, biome_5))),
    ci_divisio = list(na.omit(c(ci_divisio_1, ci_divisio_2, ci_divisio_3, ci_divisio_4, ci_divisio_5))),
    ci_divis_1 = list(na.omit(c(ci_divis_1_1, ci_divis_1_2, ci_divis_1_3, ci_divis_1_4, ci_divis_1_5))),
    country = list(na.omit(c(country_1, country_2, country_3, country_4, country_5))),
    geographic = list(na.omit(c(geographic_1, geographic_2, geographic_3, geographic_4, geographic_5))),
    inter_cat = list(na.omit(c(inter_cat_1, inter_cat_2, inter_cat_3, inter_cat_4, inter_cat_5))),
    interventi = list(na.omit(c(interventi_1, interventi_2, interventi_3, interventi_4, interventi_5))),
    interven_1 = list(na.omit(c(interven_1_1, interven_1_2, interven_1_3, interven_1_4, interven_1_5))),
    interven_2 = list(na.omit(c(interven_2_1, interven_2_2, interven_2_3, interven_2_4, interven_2_5))),
    gazettemen = list(na.omit(c(gazettemen_1, gazettemen_2, gazettemen_3, gazettemen_4, gazettemen_5))),
    iucn_categ = list(na.omit(c(iucn_categ_1, iucn_categ_2, iucn_categ_3, iucn_categ_4, iucn_categ_5))),
    tree_plant = list(na.omit(c(tree_plant_1, tree_plant_2, tree_plant_3, tree_plant_4, tree_plant_5))),
    restoratio = list(na.omit(c(restoratio_1, restoratio_2, restoratio_3, restoratio_4, restoratio_5))),
    forest_age = list(na.omit(c(forest_age_1, forest_age_2, forest_age_3, forest_age_4, forest_age_5))),
    restor_fin = list(na.omit(c(restor_fin_1, restor_fin_2, restor_fin_3, restor_fin_4, restor_fin_5))),
    local_land = list(na.omit(c(local_land_1, local_land_2, local_land_3, local_land_4, local_land_5))),
    new_or_con = list(na.omit(c(new_or_con_1, new_or_con_2, new_or_con_3, new_or_con_4, new_or_con_5))),
    new_or_c_1 = list(na.omit(c(new_or_c_1_1, new_or_c_1_2, new_or_c_1_3, new_or_c_1_4, new_or_c_1_5))),  
    ci_sls_1 = list(na.omit(c(ci_sls_1_1, ci_sls_1_2, ci_sls_1_3, ci_sls_1_4, ci_sls_1_5))),
    ci_sls_2 = list(na.omit(c(ci_sls_2_1, ci_sls_2_2, ci_sls_2_3, ci_sls_2_4, ci_sls_2_5))),
    reporting = list(na.omit(c(reporting_1, reporting_2, reporting_3, reporting_4, reporting_5))),
    improved_m = list(na.omit(c(improved_m_1, improved_m_2, improved_m_3, improved_m_4, improved_m_5))),
    ci_start_d = list(na.omit(c(ci_start_d_1, ci_start_d_2, ci_start_d_3, ci_start_d_4, ci_start_d_5))),
    ci_end_dat = list(na.omit(c(ci_end_dat_1, ci_end_dat_2, ci_end_dat_3, ci_end_dat_4, ci_end_dat_5)))
  ) %>% 
  ungroup() %>% 
  dplyr::select(n_overlaps, 
                colnames(fields),
                ecosystem, tstor_ic, ha_high_ic, ha_ic, tonnes_ha_ic, tstor_blue_ic)

saveRDS(ic_int_w_fields, "results/FY21_ImpactIndicators_IrrecoverableCarbon_Overlaps.rds")
