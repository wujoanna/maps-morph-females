#################################
# 1c - Process bird range map data
#################################


# set dirs ----------------------------------------------------------------

RUN_DATE <- '2024-07-17'
dir <- "/Users/joannawu/Library/CloudStorage/Box-Box/Documents/AcademicScience/Projects/MorphFemales/"
maps_date <- '2024-05-24'
run_date <- '2024-05-24'


# Load packages ----------------------------------------------------------

library(sf)
library(tidyverse)


# read in data -------------------------------------------------------------

#read in data
#filter by species
mdata <- readRDS(paste0(dir, 'Data/L1/MAPS-main-', maps_date, '.rds'))
#unique species
MAPS_usp <- unique(gsub('_', ' ', mdata$sci_name))

#range data - takes some time to read in
BL_data <- sf::st_read(dsn = paste0(dir, 'Data/L0/BOTW_2022/BOTW.gdb'))
#unique species
bl_usp <- unique(BL_data$sci_name)


# which don't match -------------------------------------------------------

#where matches are
nin <- which(MAPS_usp %in% bl_usp)

#names df
names_df <- data.frame(MAPS_sci_name = MAPS_usp, BL_sci_name = NA)
names_df$BL_sci_name[nin] <- MAPS_usp[nin]

#no matches
nm <- dplyr::filter(names_df, is.na(BL_sci_name))


# fill in NA ------------------------------------------------

# Dryobates villosus -> Leuconotopicus villosus
bl_usp[grep('villosus', bl_usp)]
idx1 <- which(names_df$MAPS_sci_name == nm$MAPS_sci_name[1])
names_df$BL_sci_name[idx1] <- 'Leuconotopicus villosus'

# Icterus bullockii -> Icterus bullockiorum
bl_usp[grep('Icterus', bl_usp)]
idx2 <- which(names_df$MAPS_sci_name == nm$MAPS_sci_name[2])
names_df$BL_sci_name[idx2] <- 'Icterus bullockiorum'

#write out to csv
write.csv(names_df, file = paste0(dir, 'Data/L1/range_names_key-', 
                                  run_date, '.csv'), row.names = FALSE)


# filter ranges -----------------------------------------------------------

#read in modified csv (names key)
names_df_mod <- read.csv(file = paste0(dir, 
                                       'Data/L1/range_names_key-', run_date, '.csv'))

#species names in format to pass to query
sn <- paste0(dQuote(names_df_mod$BL_sci_name), collapse = ',')

#range data for subset of species
BL_data_filt <- dplyr::filter(BL_data, sci_name %in% names_df_mod$BL_sci_name)

saveRDS(BL_data_filt, paste0(dir, 'Data/L1/BOTW_ranges-', run_date, '.rds'))

