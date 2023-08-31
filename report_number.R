library(tidyverse)

csv <- read_csv("avoided_emissions/output_raw_by_site.csv")

true <- csv %>% 
  filter(treatment == TRUE) %>% 
  rename(true = Emissions_MgCO2e)

false <- csv %>% 
  filter(treatment == FALSE) %>% 
  rename(false = Emissions_MgCO2e) %>% 
  dplyr::select(CI_ID, Data_Year, year, false)

comb <- true %>% 
  left_join(false) %>% 
  mutate(value = true - false) %>% 
  filter(value > 0) %>% 
  filter(year == 2022)

# number that Howard reports
scales::comma(sum(comb$value))
