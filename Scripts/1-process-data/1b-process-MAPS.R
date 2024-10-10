####################
# 1b - Process MAPS data
#
####################


# set dirs ----------------------------------------------------------------

RUN_DATE <- "2024-07-15"
dir <- "/Users/joannawu/Library/CloudStorage/Box-Box/Documents/AcademicScience/Projects/MorphFemales/"


# load packages -----------------------------------------------------------

library(tidyverse)
library(foreign)
library(terra)
library(sf)


# read in data -------------------------------------------------------------

#banding data
maps_b <- foreign::read.dbf(paste0(dir, '/Data/L0/MAPS_data/1016-MAPS_band&effort/1016B19.DBF'), as.is = TRUE)

#monitoring data
maps_m <- foreign::read.dbf(paste0(dir, '/Data/L0/MAPS_data/1016-MAPS_band&effort/1016M19.dbf')) %>%
  dplyr::select(station = STATION,
                sp_code = SPEC,
                year = YR,
                yr_br_status = YS)

#station data
maps_s <- read.csv(paste0(dir, '/Data/L0/MAPS_data/CntrlStations/STATIONS.csv'), 
                   skipNul = TRUE, as.is = TRUE) %>%
  dplyr::select(station = STATION,
                name = NAME,
                lat = DECLAT,
                lng = DECLNG,
                stratum = STRATUM,
                bcr = BCR,
                habitat = HABITAT)

#species codes
species_codes <- read.csv(paste0(dir, '/Data/L0/MAPS_data/IBP-AOS-LIST23.csv'),
                          as.is = TRUE) %>%
  dplyr::select(sp_code = SPEC,
                common_name = COMMONNAME,
                sci_name_raw = SCINAME)


# merge and filter --------------------------------------------------------------

maps1 <- dplyr::select(maps_b, 
                       station = STATION,
                       capture_code = C,
                       band_id = BAND,
                       sp_code = SPEC,
                       age_code = AGE,
                       sex = SEX,
                       cloacal_pro = CP,
                       brood_patch = BP,
                       fat_content = F,
                       body_molt = BM,
                       flight_molt = FM,
                       feather_wear = FW,
                       wing_chord = WNG,
                       mass = WEIGHT,
                       date = DATE,
                       anet = ANET,
                       disp = DISP,
                       standard_effort = N) %>%
  dplyr::mutate(ld = lubridate::ymd(date),
                year = lubridate::year(ld),
                #ordinal date
                day = lubridate::yday(ld),
                #add age class - AHY (adult), young of the year (juv), or unknown (NA) 
                age_class = ifelse(age_code %in% c(5,1,6,7,8), 'adult',
                                   ifelse(age_code %in% c(4, 2), 'juv', NA))) %>%
  dplyr::left_join(maps_s, by = 'station') %>%
  dplyr::left_join(maps_m, by = c('station', 'sp_code', 'year')) %>%
  dplyr::left_join(species_codes, by = 'sp_code') %>%
  #only specific capture codes
  dplyr::filter(capture_code %in% c('N', 'R', 'U'),
                #only species that were recorded breeding at station in given year
                yr_br_status == 'B',
                #only during 'breeding season'
                day >= 121, day <= 220,
                #remove clearly erroneous records
                !is.na(wing_chord), wing_chord > 0,
                !is.na(mass), mass > 0,
                lat > 0, lat < 99,
                #only known adult/juvs
                age_class %in% c('adult', 'juv'))

# From MAPS documentation:

#C (capture_code): Only codes N,R,U are used for analyses according to MAPS docs
#N - newly banded bird
#R - recaptured bird
#U - unbanded bird

#N (standard_effort - Indicator to include in productivity and survivorship analyses): Okay to use these for analyses that do not require to control for effort
#O - not caught at MAPS station
#T - time outside normal MAPS operation for that station for that year
#S - caught within MAPS station boundary but not in a MAPS net
#D - date outside of MAPS periods
#- - record examined with current MAPS analytical procedure (taken to be standard capture methods)
#+ - record examined with preliminary MAPS analytical procedure

#Age:
#0 - unknown
#4 - local (young bird incapable of flight)
#2 - hatching-year bird
#1 - after hatching-year bird
#5 - second-year bird
#6 - after second-year bird
#7 - third-year bird
#8 - after third-year bird

#Sex:
#M - male
#F - female
#U - unknown
#X - unattempted

#CP (cloacal_pro - cloacal protuberance):
#0 - none
#1 - small
#2 - medium
#3 - large

#BP (brood_patch - Brood patch:
#0 - none
#1 - smooth
#2 - vascularized
#3 - heavy
#4 - wrinkled
#5 - molting

#F (fat_content):
#0 - none
#1 - trace
#2 - light
#3 - half
#4 - full
#5 - bulging
#6 - greatly bulging
#7 - very excessive

#YS (yr_br_status - year specific breedin status):
#B - breeder (at least one ind was summer resident at station)
#L - likely breeder (at least one ind was suspected summer resident)
#T - transient (station is within breeding rng of species, but no individual of the species was a summer resident at the station)
#A - altitudinal disperser (species breeds only at lower elevations than that of the station and which disperses to higher elevations after breeding)
#H - high altitudinal disperser (species breeds usually designated an altitudinal disperser. However, has resided during the height of the breeding season (not just during the post- breeding period) in a given year above its normal breeding elevation.
#M - migrant (station is not within the breeding range of the species, and the species was not a summer resident)
#E - extralimital breeder (one or more individuals of the species was a summer resident at the station, but the station lies outside of the normal breeding range of the species)
#- - absent (no evidence of species in data; presumably absent from station during year in question)
#? - uncertain species identification or band number (no breeding status assigned)
## - station operated this year, but breeding status determinations were not made for species that were not captured; used only for species without capture records
#D - species was only encountered at the station outside of the MAPS season, but the station lies within breeding range of the species.
#W - species was only encountered at the station outside of the MAPS season, and the station lies outside of the breeding range of species.
#@ - Breeding Status List is missing or incomplete for these species this year.

# DISP: disposition of birds upon release or after capture
# O - old (healed) injury
# M - malformed (deformity such as crossed mandibles)
# W - wing injury
# L - leg injury
# T - tongue injury
# E - eye injury
# B - body injury
# I - illness/infection/disease
# S - stress or shock
# P - predation (death due to predation)
# D - dead (death due to causes other than predation or removed permanently from station)
# R - band removed from bird and then bird released bandless
# " " - blank, bird released alive, uninjured


# true age ----------------------------------------------------------------

#could add true age, but maybe not useful - have code somewhere to do this


# lump sub-species --------------------------------------------------------

#unique species
usp <- sort(unique(maps1$sci_name_raw))
#species with slashes
spsl <- grep('/', usp)
#hybrids
hsp <- grep(' x ', usp)
#remove
usp2 <- usp[-c(spsl, hsp)]

#lump subspecies
l2w <- which(stringr::str_count(usp2, '\\S+') > 2)
ndf1 <- data.frame(sci_name_raw = usp2[l2w], 
                   sci_name = stringr::word(usp2[l2w], 1, 2, sep = ' '))
ndf2 <- rbind(ndf1, data.frame(sci_name_raw = usp2[-l2w], 
                               sci_name = NA))
nm_na <- which(is.na(ndf2$sci_name))
ndf2$sci_name[nm_na] <- ndf2$sci_name_raw[nm_na]

maps2 <- dplyr::left_join(maps1, ndf2, by = 'sci_name_raw') %>%
  dplyr::select(sci_name, 
                common_name,
                station,
                lat,
                lng,
                year, 
                day,
                band_id,
                sex,
                age_class,
                mass,
                wing_chord,
                fat_content,
                brood_patch,
                cloacal_pro,
                habitat,
                bcr) %>%
  dplyr::mutate(sci_name = gsub(' ', '_', sci_name)) %>%
  #only first capture of individual in a season
  dplyr::group_by(band_id, year) %>%
  dplyr::slice_min(day, with_ties = FALSE) %>%
  dplyr::ungroup()


# remove anomalous morphological values for adults and juveniles-----------------------------

thresh <- dplyr::group_by(maps2, sci_name, age_class, sex) %>%
  dplyr::summarize(med_wc = median(wing_chord),
                   med_mass = median(mass),
                   MAD_wc = mad(wing_chord),
                   MAD_mass = mad(mass)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(low_wc = med_wc - (5 * MAD_wc),  ## remove all values more than 3 MAD away (v conservative) from median - Leys et al. 2013
                high_wc = med_wc + (5 * MAD_wc),
                low_mass = med_mass - (5 * MAD_mass),
                high_mass = med_mass + (5 * MAD_mass))

# Filter out outliers
# only first capture of an individual in a season
# only adults
# at least 375 captures for each species
maps3 <- dplyr::left_join(maps2, thresh, 
                          by = c('sci_name', 'age_class', 'sex')) %>%
  dplyr::filter(!(wing_chord < low_wc | 
             wing_chord > high_wc | 
             mass < low_mass | 
             mass > high_mass)) %>%
  dplyr::filter(!is.na(sci_name), 
                !is.na(band_id),
                age_class == 'adult') %>%
  dplyr::group_by(sci_name) %>%
  dplyr::filter(n() >= 375) %>%
  #sort
  dplyr::arrange(sci_name, station, year) %>%
  #add species ID
  dplyr::mutate(sp_id = dplyr::cur_group_id()) %>%
  dplyr::ungroup() %>%
  dplyr::select(-med_wc,
                -med_mass,
                -MAD_wc,
                -MAD_mass,
                -low_wc,
                -high_wc,
                -low_mass,
                -high_mass)


# locations for daymet extraction -----------------------------------------

#write locations to file to extract daymet (temp/precip) data
#arrange according to: https://github.com/ornldaac/daymet-single-pixel-batch
#or use daymetr package

write.table(unique(data.frame(station = maps3$station, 
                              latitude = maps3$lat, 
                              longitude = maps3$lng)),
            paste0(dir, 'Data/L0/daymet_query_locations.csv'), 
            sep = ',', row.names = FALSE, col.names = TRUE)


# add elev ----------------------------------------------------------------

#read in mosaic tif
elev_mosaic <- terra::rast(paste0(dir, 'Data/L1/elev_mosaic.tif'))


# #plot
# plot(elev_mosaic)
# ull <- unique(maps4[, c('lng', 'lat')])
# points(ull$lng, ull$lat, col = rgb(0,0,0,0.5), pch = '.')

#convert lat lng to spatial coords
coords_sp <- sf::st_as_sf(maps3, coords = c('lng', 'lat'), crs = 4326)

#extract from raster and add to maps df
maps3$elev <- terra::extract(elev_mosaic, coords_sp, df = FALSE)$elev_mosaic


# data stats --------------------------------------------------------------

#DATA STATS
# num samples - 555,873
NROW(maps3)
# num stations - 1136
length(unique(maps3$station))
# num species - 142
length(unique(maps3$sci_name))
# num by sex and age class
dplyr::count(maps3, sex, age_class)
# range years - 1989-2018
range(maps3$year)
# range lat - 26.1-69.4
range(maps3$lat)
# range elev - 0m-2996m
range(maps3$elev)


# save rds ----------------------------------------------------------------

#rearrange columns and save out
dplyr::relocate(maps3, sp_id, .before = station) %>%
  dplyr::relocate(elev, .before = year) %>%
  saveRDS(paste0(dir, '/Data/L1/MAPS-main-', RUN_DATE, '.rds'))

