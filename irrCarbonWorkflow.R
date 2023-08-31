# Prep data -----

library(tidyverse)
library(raster)
library(terra)
library(sf)
library(janitor)
library(exactextractr)

eco <- rast("data/ecosystems/ecosystems.tif") 

# total tonnes of carbon per hectare
total_ic <- rast("data/irrecoverable_carbon/Irrecoverable_C_Total_2018.tif") %>% 
  resample(eco, method = "bilinear")
# checked multiple methods of resampling + combining with area and this was the most accurate
# at achieving 139B tonnes

# convert to tonnes per pixel
area <- raster::area(raster(total_ic)) * 100 #km2 to ha

tonnes_ic <- raster(total_ic) * area

# ic total with vals below 25 tonnes/ha as NA
high_ic <- rast("data/irrecoverable_carbon/ic_above_25.tif") %>% 
  resample(eco, method = "bilinear")

# make raster of 'any' ic (above 0.01)
any_ic <- total_ic
any_ic[any_ic > 0.01] <- 1
any_ic[any_ic <= 0.01] <- NA

# change ic_high to binary
high_ic[high_ic > 0] <- 1

# need 'high' and 'any' to reflect has/pixel
# change ic_high and ic_any to represent area

high_area <- cellSize(high_ic, unit = "ha")

any_area <- cellSize(any_ic, unit = "ha")


# create stack
ic_stack <- c(total_ic, high_area, any_area)  

names(ic_stack) <- c("tstor_ic", "ha_high_ic", "ha_ic")

# check total to confirm it matches with our 140B tonnes total in Noon et al
global(ic_stack$tstor_ic, "sum", na.rm = TRUE)

# divvy up irrecoverable carbon by ecosystem

### ecosystems to mask irr carbon to

eco <- rast("data/ecosystems/ecosystems.tif") 
eco_codes <- data.frame(cats(eco$class)) %>%
  rename('ecosystem_name' = class)  %>% 
  drop_na()

eco_seg <- eco %>%
  terra::segregate()

eco_seg[eco_seg == 0] <- NA

# for each ic layer
for (i in 1:nlyr(ic_stack)){
  
  ic_lyr <- ic_stack[[i]]
  og_name <- names(ic_lyr)
  ic_eco_stack <- ic_lyr
  
  # mask irr carbon to each ecosystem layer then add to stack
  for(e in 1:8){
    eco_mask <- eco_seg[[e]]
    class <- eco_codes$ecosystem_name[e] %>% 
      make_clean_names()
    ic_masked <- ic_lyr[[1]] %>%
      terra::mask(mask = eco_mask)
    names(ic_masked) <- paste0(og_name, "_", class)
    add(ic_eco_stack) <- ic_masked
    remove(eco_mask)
    remove(ic_masked)
    gc()
  }
  
  # then save stack for that ic 
  writeRaster(
    ic_eco_stack, 
    filename = paste0("data/irrecoverable_carbon/", og_name, "_ecosystem_stack.tif"),
    overwrite = TRUE)
}



# Extract data by site ----

### irrecoverable carbon with layers masked to each ecosystem
ic_total <- rast("data/irrecoverable_carbon/tstor_ic_ecosystem_stack.tif")
ic_high <- rast("data/irrecoverable_carbon/ha_high_ic_ecosystem_stack.tif") # ha w/ ic over 25
ic_any <- rast("data/irrecoverable_carbon/ha_ic_ecosystem_stack.tif") # ha w/ ic over 0.01

# note this is a new 'ic_stack' with layers for each ecosystem
ic_stack <- c(ic_total, ic_high, ic_any)

### define crs 
my_crs <- crs(ic_stack)

##### Prepare shapefiles and intersections
#####

# ignore errors around spherical geometries
sf::sf_use_s2(FALSE)

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

