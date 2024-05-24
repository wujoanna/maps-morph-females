#######################
# process daymet data
#
# https://daymet.ornl.gov/
# 
#######################


# input -------------------------------------------------------------------

dir <- '~/Work/Research/Projects/maps-morph-females/'
daymet_process_date <- '2024-05-24'


# load packages -----------------------------------------------------------

library(tidyverse)
library(daymetr)


# download point values for daymet from ORNL ------------------------------

#https://cran.r-project.org/web/packages/daymetr/vignettes/daymetr-vignette.html
#dl data, convert to wide format
dm <- daymetr::download_daymet_batch(file_location = 
                                       paste0(dir, 'Data/L0/daymet_query_locations.csv'),
                                     start = 1980,
                                     end = 2018,
                                     internal = TRUE,
                                     simplify = TRUE) %>% 
  dplyr::filter(measurement %in% 
                  c('prcp..mm.day.',
                    'tmin..deg.c.',
                    'tmax..deg.c.')) %>%
  #dplyr::select(-altitude, -tile, -latitude, -longitude) %>%
  tidyr::pivot_wider(names_from = measurement, values_from = value) %>%
  dplyr::rename(tmax = tmax..deg.c.,
                tmin = tmin..deg.c.,
                prcp = prcp..mm.day.)


# save daymet RDS -----------------------------------------------------------

#daily env for each station
saveRDS(dm, paste0(dir, 'Data/L0/daymet-dl-', daymet_process_date, '.rds'))
# dm <- readRDS(paste0(dir, 'Data/L0/daymet-dl-', daymet_process_date, '.rds'))


# average over specified timeframe ----------------------------------------

#average prcp, tmin, and tmax for each station over specified time frame

#add date
dm$date <- as.Date(paste0(dm$year, '-01-01')) + dm$yday

#non-leap year
Jan_env <- dplyr::group_by(dm, site, year) %>%
  dplyr::filter(yday >= 1, yday <= 31) %>%
  dplyr::summarize(Jan_prcp = mean(prcp),
                   Jan_tmin = mean(tmin),
                   Jan_tmax = mean(tmax)) %>%
  dplyr::ungroup()

Feb_env <- dplyr::group_by(dm, site, year) %>%
  dplyr::filter(yday >= 32, yday <= 59) %>%
  dplyr::summarize(Feb_prcp = mean(prcp),
                   Feb_tmin = mean(tmin),
                   Feb_tmax = mean(tmax)) %>%
  dplyr::ungroup()

Mar_env <- dplyr::group_by(dm, site, year) %>%
  dplyr::filter(yday >= 60, yday <= 90) %>%
  dplyr::summarize(Mar_prcp = mean(prcp),
                   Mar_tmin = mean(tmin),
                   Mar_tmax = mean(tmax)) %>%
  dplyr::ungroup()

Apr_env <- dplyr::group_by(dm, site, year) %>%
  dplyr::filter(yday >= 91, yday <= 120) %>%
  dplyr::summarize(Apr_prcp = mean(prcp),
                   Apr_tmin = mean(tmin),
                   Apr_tmax = mean(tmax)) %>%
  dplyr::ungroup()

May_env <- dplyr::group_by(dm, site, year) %>%
  dplyr::filter(yday >= 121, yday <= 151) %>%
  dplyr::summarize(May_prcp = mean(prcp),
                   May_tmin = mean(tmin),
                   May_tmax = mean(tmax)) %>%
  dplyr::ungroup()

June_env <- dplyr::group_by(dm, site, year) %>%
  dplyr::filter(yday >= 152, yday <= 181) %>%
  dplyr::summarize(June_prcp = mean(prcp),
                   June_tmin = mean(tmin),
                   June_tmax = mean(tmax)) %>%
  dplyr::ungroup()

July_env <- dplyr::group_by(dm, site, year) %>%
  dplyr::filter(yday >= 182, yday <= 212) %>%
  dplyr::summarize(July_prcp = mean(prcp),
                   July_tmin = mean(tmin),
                   July_tmax = mean(tmax)) %>%
  dplyr::ungroup()

Aug_env <- dplyr::group_by(dm, site, year) %>%
  dplyr::filter(yday >= 213, yday <= 243) %>%
  dplyr::summarize(Aug_prcp = mean(prcp),
                   Aug_tmin = mean(tmin),
                   Aug_tmax = mean(tmax)) %>%
  dplyr::ungroup()

Sep_env <- dplyr::group_by(dm, site, year) %>%
  dplyr::filter(yday >= 244, yday <= 273) %>%
  dplyr::summarize(Sep_prcp = mean(prcp),
                   Sep_tmin = mean(tmin),
                   Sep_tmax = mean(tmax)) %>%
  dplyr::ungroup()

Oct_env <- dplyr::group_by(dm, site, year) %>%
  dplyr::filter(yday >= 274, yday <= 304) %>%
  dplyr::summarize(Oct_prcp = mean(prcp),
                   Oct_tmin = mean(tmin),
                   Oct_tmax = mean(tmax)) %>%
  dplyr::ungroup()

Nov_env <- dplyr::group_by(dm, site, year) %>%
  dplyr::filter(yday >= 305, yday <= 334) %>%
  dplyr::summarize(Nov_prcp = mean(prcp),
                   Nov_tmin = mean(tmin),
                   Nov_tmax = mean(tmax)) %>%
  dplyr::ungroup()

Dec_env <- dplyr::group_by(dm, site, year) %>%
  dplyr::filter(yday >= 335, yday <= 365) %>%
  dplyr::summarize(Dec_prcp = mean(prcp),
                   Dec_tmin = mean(tmin),
                   Dec_tmax = mean(tmax)) %>%
  dplyr::ungroup()

env_mrg <- dplyr::left_join(Jan_env, Feb_env) %>%
  dplyr::left_join(Mar_env) %>%
  dplyr::left_join(Apr_env) %>%
  dplyr::left_join(May_env) %>%
  dplyr::left_join(June_env) %>%
  dplyr::left_join(July_env) %>% 
  dplyr::left_join(Aug_env) %>%
  dplyr::left_join(Sep_env) %>% 
  dplyr::left_join(Oct_env) %>% 
  dplyr::left_join(Nov_env) %>% 
  dplyr::left_join(Dec_env)


# save RDS object -----------------------------------------------------

saveRDS(env_mrg, paste0(dir, 'Data/L1/daymet-main-', 
                        daymet_process_date, '.rds'))

