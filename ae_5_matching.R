library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(tidyverse)
library(foreach)
library(optmatch)
library(lubridate)
library(biglm)
library(tictoc)

# Creates "ae_sites_by_year" file
## and 'ae_output_raw_by_site_2024' df 

options("optmatch_max_problem_size"=Inf)
# This package trys to prevent us from running ridiculous parameters
# But we are trying to circumvent that limit 

MAX_TREATMENT <- 1000
CONTROL_MULTIPLIER <- 50
 

###############################  
######## FUNCTIONS ############
###############################  

# Function to allow rbinding dataframes with foreach even when some dataframes 
# may not have any rows
## adding / tacking on rows
## so maybe empty dataframes, so a foreach, you can't bind an empty dataframe
## it's letting you combine empty dataframes with those that do 

foreach_rbind <- function(d1, d2) {
    if (is.null(d1) & is.null(d2)) {
        return(NULL)
    } else if (!is.null(d1) & is.null(d2)) {
        return(d1)
    } else if (is.null(d1) & !is.null(d2)) {
        return(d2)
    } else  {
        return(bind_rows(d1, d2))
    }
}

# Basic function to extract variable names from a formula object
## which variables to use -- this variables let's me construct a formula name
## from the other variables -- if you knew in advance what covariates you're using,
## you can just write the equation
get_names <- function(f) {
    f <- paste0(as.character(f), collapse=' ')
    v <- strsplit(f, split='[+ ~]')[[1]]
    v <- v[v != '']
    gsub('strata\\(([a-zA-Z_]*)\\)', '\\1', v)
}


get_matches <- function(d, dists) {
    # If the controls are too far from the treatments (due to a caliper) then 
    # the matching may fail. Can test for this by seeing if subdim runs 
    # successfully
  ## Folks like Seb and Anand -- within matching, you can include a calipher -- we will only 
  ## include a match if it's within a certain distance -- our current version doesn't , but
  ## a calipher makes it more likely to NOT have a match -- which then leadas to not having any matches
  ## and having matches fail -- if you run everything in parallel, trying to uderstand if it's possible to not have a match for something
    subdim_works <- tryCatch(is.data.frame(subdim(dists)),
                             error=function(e) return(FALSE))
    if (subdim_works) {
        #m <- pairmatch(dists, data=d, remove.unmatchables = TRUE)
        m <- fullmatch(dists, min.controls=1, max.controls=1, data=d)
        d <- d[matched(m), ]
    } else {
        d <- data.frame()
    }
    return(d)
}

match_ae <- function(d, f) {
  ## Group involves the 3 matches -- region, ecoregion, in/out PA
    m <- foreach(this_group=unique(d$group), .combine=foreach_rbind) %do% {
        this_d <- filter(d, group == this_group)
        # Calculate propensity scores with a GLM, or else use Mahalanobis 
        # distance if there aren't enough points to run a glm
        if (sum(this_d$treatment) > 30) {
            model <- glm(f, data=this_d, family=binomial()) #main line -- used generalized linear model, logistic, binomial using formula F (using getnames function)
            dists <- match_on(model, data=this_d) #main line 2 -- this is where the match happens
        } else {
            dists <- match_on(f, data=this_d) ## second way -- if you have more than 30 points, run a linear regression; if you have 
            ## fewer than 30 points, then less chance of having a linear model, so just calculating within the data, how far it is 
            ## one thing to think about -- those sites are likely to have fewer poly matches -- interesting to look at what those
            ## reuslts look like
        }
        return(get_matches(this_d, dists))
    }
    # Need to handle the possibility that there were no matches for this 
    # treatment, meaning d will be an empty data.frame
    if (nrow(m) == 0) {
        return(NULL)
    } else {
        return(m)
    }
}

#############################################  
######## LOAD SITES & COVARIATES ############
#############################################

#Load treatment key
treatment_key <- readRDS('ae_output/treatment_cell_key.RDS')

#Load sites
readRDS('data/avoided_emissions/sites_FY2024.RDS') %>%
    dplyr::select(-geometry) %>%
    as_tibble() -> sites

sites <- sites %>% mutate(CI_Start_Year=year(CI_Start_Date_clean),
       CI_End_Year=ifelse(is.na(year(CI_End_Date_clean)), 2099, year(CI_End_Date_clean)))

# Filter to include only values of group that appear in the treatment pixels, 
# and to not include values that appear only in the treatment pixels

## BC of the way we pulled the data in the previous script -- we may have
## only want to compare pixels that are iN PAs -- filter out all of the stuff
## that's extra data that we don't need -- saving memory

filter_groups <- function(vals) {
    vals$group <- interaction(vals$region, vals$ecoregion, vals$pa)
    vals <- filter(vals, group %in% unique(filter(vals, treatment)$group))
    treatment_groups <- unique(filter(vals, treatment)$group)
    control_groups <- unique(filter(vals, !treatment)$group)
    vals <- filter(vals, group %in% treatment_groups[treatment_groups %in% control_groups])
    # Filter out values of group that appear ONLY in the treatment pixels
    vals$group <- droplevels(vals$group)
    return(vals)
}



################################## 
######## RUN MATCHING ############
##################################

sites <- sites %>%
  filter(map_lgl(CI_ID, ~ !any(grepl("BNA", .))))

treatment_key <- readRDS('ae_output/treatment_cell_key.RDS')

# Add "Data_Year" column to treatments
treatment_key$Data_Year <- 2024

treatment_key <- treatment_key %>% mutate(CI_End_Year=ifelse(is.na(year(CI_End_Year)), 2099, CI_End_Year))


## this is randomized component --
## this maintains the same matches
## this makes it repeateable
set.seed(31)


## This is calling the "get matches" function
## Looping across years (don't really need this anymore)

ae <- foreach(this_year=unique(treatment_key$Data_Year),
              .combine=foreach_rbind, .inorder=FALSE) %do% {
        foreach(this_CI_ID=unique(treatment_key$CI_ID), ## loop across each CI site
            .combine=foreach_rbind, .inorder=FALSE) %do% {
        tic() #Measures time for length of processing
        ###############
        
        set.seed(31)
              # Load datasets
        
        site <- filter(sites,
                       CI_ID == this_CI_ID,
                       Data_Year == this_year)
        if (file.exists(paste0('ae_output/m_', this_CI_ID, '_', this_year, '.RDS'))) {
            print(paste0('Skipping ', this_CI_ID, ' for year ', this_year, '. Already processed.'))
            return(NULL)
        } else {
            print(paste0('Processing ', this_CI_ID, ' for year ', this_year, '.'))
        }

        treatment_cell_IDs <- filter(treatment_key,
                                     CI_ID == this_CI_ID,
                                     !is.na(region),
                                     Data_Year == this_year)
        n_treatment_cells_total <- nrow(treatment_cell_IDs)
        if (n_treatment_cells_total == 0) {
            print(paste0('Skipping ', this_CI_ID, ' for year ', this_year, '. No treatment cells.'))
            return(NULL)
        }
        
        ## Reading in covariates
        vals <- foreach(this_region = unique(treatment_cell_IDs$region),
                        .combine=rbind) %do% {
            v <- readRDS(paste0('ae_output/treatments_and_controls_', this_region, 
                                '.RDS'))
            filter(v, region == this_region)
        }

        ## joining treatment cells and shows
        ## which cells are treatment (treatment true vs. cells where we didn't)
        vals %>% full_join(
                treatment_cell_IDs %>%
                    dplyr::select(cell, Data_Year) %>%
                    mutate(treatment=TRUE)
                , by='cell') -> vals
        vals$treatment <- as.logical(vals$treatment)
        vals$treatment[is.na(vals$treatment)] <- FALSE
        vals$Data_Year <- this_year

        # Remove areas falling within another CI site from the control sample 
        # (but DON'T remove those areas falling within this site)
        # We don't want to measure sites that we work in a DIFFERENT CI site
        filter(vals,
               !(cell %in% filter(treatment_key,
                                  !(cell %in% filter(treatment_key,
                                                     CI_ID == this_CI_ID)$cell))$cell)) -> vals

        ################
        # Setup grouping
        
        # Eliminate any pixels with NAs in group variables (happens occasionally 
        # where polygons overlap some ocean, leading to undefined ecoregion, for 
        # example)
        n_filtered <- nrow(vals)
        print(paste0("Total rows:",n_filtered))
        
        ## areas where we don't have region, ecoregion, or PA
        vals <- filter(vals,
                       !is.na(region),
                       !is.na(ecoregion),
                       !is.na(pa))
        print(paste0("Filtered out Pixels:",nrow(vals)))
        n_filtered <- n_filtered - nrow(vals)
        
        ## tells how many rows we filtered out
        if (n_filtered > 0) {
            print(paste0(this_CI_ID, ': Filtered ', n_filtered, ' rows due to missing data in grouping variables.'))
        }
        # Eliminate any groups that are only in the control pixels, or only in the 
        # treatment pixels
        vals <- filter_groups(vals)
        sample_sizes <- vals %>%
            count(treatment, group)
 
        # Sample the treatment cells if there are more than MAX_TREATMENT pixels,
        # and the control cells if there are more than CONTROL_MULTIPLIER *
        # MAX_TREATMENT pixels
        
        # control multiplier -- if you have 1000 treatment cells, then pull out
        # 50,000 cells to compare it to 
        # 1:1 match -- if we have 1000 cells, we want 1000 matches; but we draw from
        # a pool of 50,000 -- max that you'd pull 
        # either bigger or smaller -- anything you change would change the result 
        # not a worse result, just different 
        set.seed(31)
        
        treatment <-  filter(vals, treatment)  %>%
          group_by(group) %>%
          sample_n(min(MAX_TREATMENT, n()))
        
        set.seed(31)
        
        controls <-filter(vals, !treatment)  %>%
          group_by(this_group=group) %>%
          sample_n(min(CONTROL_MULTIPLIER * filter(sample_sizes,
                                                   treatment == TRUE,
                                                   group == this_group[1])$n,
                       n()))


        bind_rows(treatment, controls) %>%
            ungroup() %>%
            dplyr::select(-this_group) -> vals
        # Refilter in case any groups were lost due to the sampling
        vals <- filter_groups(vals)
  
      # Project all items to cylindrical equal area
      # d_crop <- projectRaster(d_crop, crs=CRS('+proj=cea'), method='ngb')
  
      ################
      # Add defor data
        
        ## before in the powerpoint -- matching on 5 years of deforestation
        ## going to be noise from deforestation data
        
      # For sites that were established in or after 2005, match on the five years
      # of deforestation data preceding the year of establishment. For sites
      # estab prior to 2005, don't match on defor rate
      
        
      estab_year <- year(site$CI_Start_Date_clean)
      f <- readRDS('ae_output/formula.RDS')
      if (estab_year >= 2005) {
          init <- vals[, grepl(paste0('fc20', substr(estab_year - 5, 3, 4)), names(vals))]
          final <- vals[,grepl(paste0('fc20', substr(estab_year, 3, 4)), names(vals))]
          defor_pre_intervention <- ((final - init) / init) * 100
          # Correct for division by zero in places that had no forest cover in
          # year 0
          defor_pre_intervention[init == 0] <- 0
          names(defor_pre_intervention) <- 'defor_pre_intervention'
          vals <- cbind(vals, defor_pre_intervention)
          f <- update(f, ~ . + defor_pre_intervention)
      }
      #vals <- vals %>% select(-starts_with('fc_'), -starts_with('fcc_'))

      ## this is just calculating to print -- it's already happened above but this just prints it
      sample_sizes <- vals %>%
          count(treatment, group)
      print(paste0(this_CI_ID, ': ', paste(filter(sample_sizes, treatment)$n, collapse=', '), ' treatment pixels'))
      print(paste0(this_CI_ID, ': ', paste(filter(sample_sizes, !treatment)$n, collapse=', '), ' control pixels'))
  
      ##############
      # Run matching
      
      ## actually running the matching
      
      ## if we've gotten this far and don't have anything -- nothing left to do
      if (nrow(filter(vals, treatment)) == 0) {
          print(paste0(this_CI_ID, ': No treatment values remaining after filtering'))
          return(NULL)
      } else {
        # otherwise, pulls out the data
          m <- match_ae(vals, f)
          print(paste0(this_CI_ID, ': Formatting output'))
          if (is.null(m)) {
              print(paste0(this_CI_ID, ': no matches'))
          } else {
              m$CI_ID <- this_CI_ID
              m$Data_Year <- this_year
              m <- m %>% dplyr::select(CI_ID, everything()) # m is what you get back from matching code
              # just list of pixels (adding CI ID and data year, adding back in )
              print(paste0(this_CI_ID, ': saving output'))
              m$sampled_fraction <- sum(vals$treatment) / n_treatment_cells_total ## This is the sampled fraction
              # if it's only 20% of the site, then we need to scale the numbers up 
              # Count how many cells in your treatment you have (1000 over some total number of cells -- say 10,000 cells
              # in CI site -- or 0.1 -- so then when we actually use it, we are saying that we want to multiply
              # forest loss by 10 to account for only looking at 10% of the site)
              saveRDS(m, paste0('ae_output/m_', this_CI_ID, '_', this_year, '.RDS'))
          }
      }
      toc()
      return(m)
     }
}

###############################################################################
# Load all output matches and resave in one file
###############################################################################

file_list <- list.files('ae_output', pattern ='^m_[A-Z]{3}[0-9]*_[0-9]{4}.RDS$')

m <- foreach(f=list.files('ae_output', pattern ='^m_[A-Z]{3}[0-9]*_[0-9]{4}.RDS$'), 
              .combine=foreach_rbind) %do% {
    readRDS(paste0('ae_output/', f))
}
saveRDS(m, 'ae_output/m_ALL.RDS')


###############################################################################
# Summarize results by site
###############################################################################


readRDS('data/avoided_emissions/sites_FY2024.RDS') %>%
  dplyr::select(CI_ID,
                Data_Year,
                CI_Start_Date_clean, 
                CI_End_Date_clean) %>%
  mutate(CI_Start_Year=year(CI_Start_Date_clean),
         CI_End_Year=ifelse(is.na(year(CI_End_Date_clean)), 2099, year(CI_End_Date_clean))) %>%
  dplyr::select(-CI_Start_Date_clean, -CI_End_Date_clean) -> sites


# sites %>%
#   dplyr::select(CI_ID,
#                 Data_Year,
#                 CI_Start_Date_clean, 
#                 CI_End_Date_clean) %>%
#   mutate(CI_Start_Year=year(CI_Start_Date_clean), ## THE END YEAR LINE IS SUPER IMPORTANT
#          CI_End_Year=ifelse(is.na(year(CI_End_Date_clean)), 2099, year(CI_End_Date_clean))) %>%
#   dplyr::select(-CI_Start_Date_clean, -CI_End_Date_clean) -> sites

## Reading the sites in

get_chunk <- function(d, n, n_chunks=10) {
    start_ind <- round(seq(1, length(d), length.out=n_chunks + 1))[1:(n_chunks)]
    end_ind <- c(start_ind[2:n_chunks] - 1, length(d))
    return(d[start_ind[n]:end_ind[n]])
}

# Process in chunks to save memory
n_chunks <- 20

data_files <- list.files('ae_output', pattern ='^m_[A-Z]{3}[0-9]*_[0-9]{4}.RDS$')


m_site <- foreach (i=1:n_chunks, .combine=bind_rows) %do% {
  print(paste0('Progress: ', ((i-1)/n_chunks)*100, '%'))
  tic()
  this_m <- foreach(f=get_chunk(data_files, i, n_chunks),
                    .combine=bind_rows, .inorder=FALSE) %do% {
                      
                      readRDS(paste0('ae_output/', f)) %>%
                        rename_with(~ gsub("^fc2(\\d+)", "fc_2\\1", .), starts_with("fc2")) %>%
                        dplyr::select(cell,
                               CI_ID,
                               Data_Year, 
                               treatment,
                               sampled_fraction,
                               total_biomass,
                               starts_with('fc_')) %>%
                        left_join(sites, by=c('CI_ID', 'Data_Year')) %>%
                        gather(year, forest_at_year_end, starts_with('fc_')) %>%
                        mutate(year= as.numeric(str_replace(year, 'fc_', ''))) %>%
                        group_by(CI_ID, Data_Year, cell, treatment) %>%
                        filter(between(year, CI_Start_Year[1] - 1, CI_End_Year[1])) %>% # include one year prior to project start to get initial forest cover
                        arrange(cell, year) %>%
                        mutate(forest_loss_during_year=c(NA, diff(forest_at_year_end)),
                               forest_frac_remaining = forest_at_year_end / forest_at_year_end[1],
                               biomass_at_year_end = total_biomass * forest_frac_remaining,
                               #  to convert biomass to carbon * .5
                               C_change=c(NA, diff(biomass_at_year_end)) * .5,
                               #  to convert change in C to CO2e * 3.67
                               Emissions_MgCO2e=C_change * -3.67) %>%
                        filter(between(year, CI_Start_Year[1], CI_End_Year[1])) %>% # drop year prior to project start as no longer needed
                        as_tibble() -> x
                    }
  this_m %>%
    group_by(CI_ID, Data_Year, year, treatment) %>%
    summarise(CI_Start_Year=CI_Start_Year[1],
              CI_End_Year=CI_End_Year[1],
              # correct totals for areas where only a partial sample was used 
              # by taking into account the fraction sampled
              forest_loss_ha=sum(forest_loss_during_year, na.rm=TRUE) * (1 / sampled_fraction[1]),
              Emissions_MgCO2e=sum(Emissions_MgCO2e, na.rm=TRUE) * (1 / sampled_fraction[1]),
              n_pixels=n()) %>%
    as_tibble() -> this_m_site
  
  toc()
  gc()
  return(this_m_site)
}


saveRDS(m_site, file='ae_output/output_raw_by_site.RDS')
write_csv(m_site, 'ae_output/output_raw_by_site.csv')

### Make "ae_sites_by_year" file
## Take 'output_raw_by_site_2024' df
## This df has these variables: CI_ID,  Data_Year,  year, treatment, forest_loss_ha, Emissions_MgCO2e
## for each unique year (year column) and each unique CI_ID :
## Subtract the Emissions_MgCO2e values when treatment = "FALSE" treatment from the Emissions_MgCO2e when treatment = "TRUE"
## Put this value in a new column: Emissions_avoided_MgCO2e_2022 or Emissions_avoided_MgCO2e_2023
## For forest_loss_ha, do the same thing
## Final output table should have 5 columns: CI_ID, emissions_avoided_MgCO2e_2022, emissions_avoided_MgCO2e_2023, forest_loss_avoided_ha_2022, forest_loss_avoided_ha_2023 


calculate_avoided_emissions <- function(df) {
  # Filter the data for year 2024
  data_2023 <- df %>%
    filter(year == 2023)
  
  # Calculate the differences for year 2023 and CI_ID
  emissions_avoided_2023 <- data_2023 %>%
    group_by(CI_ID) %>%
    summarize(emissions_avoided_MgCO2e_2023 = sum(Emissions_MgCO2e[treatment == "TRUE"]) - sum(Emissions_MgCO2e[treatment == "FALSE"]),
              forest_loss_avoided_ha_2023 = sum(forest_loss_ha[treatment == "TRUE"]) - sum(forest_loss_ha[treatment == "FALSE"]))
  
  # Filter the data for year 2022
  data_2022 <- df %>%
    filter(year == 2022)
  
  # Calculate the differences for year 2023 and CI_ID
  emissions_avoided_2022 <- data_2022 %>%
    group_by(CI_ID) %>%
    summarize(emissions_avoided_MgCO2e_2022 = sum(Emissions_MgCO2e[treatment == "TRUE"]) - sum(Emissions_MgCO2e[treatment == "FALSE"]),
              forest_loss_avoided_ha_2022 = sum(forest_loss_ha[treatment == "TRUE"]) - sum(forest_loss_ha[treatment == "FALSE"]))
  
  # Merge the results
  result <- left_join(emissions_avoided_2023, emissions_avoided_2022, by = "CI_ID")
  
  # Return the result
  return(result)
}



# Example usage:
output_raw_by_site_2024 <- m_site
output_raw_by_site_2024 <- output_raw_by_site
# Call the function with your dataframe
result <- calculate_avoided_emissions(output_raw_by_site_2024)

result <- ae_by_site_by_year_2024
saveRDS(result, file='data/avoided_emissions/ae_by_site_by_year_2024.RDS')
write_csv(result, 'results/summaries_2024/ae_output_raw_by_site.csv')

