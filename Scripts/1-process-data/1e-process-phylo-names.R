####################
# 1e - process phylo names
#
####################


# set dirs ----------------------------------------------------------------

RUN_DATE <- '2024-05-24'
MAPS_DATE <- '2024-05-24'
dir <- '~/Work/Research/Projects/maps-morph-females/'


# load packages -----------------------------------------------------------

library(tidyverse)


# read in data ------------------------------------------------------------

mdata <- readRDS(paste0(dir, 'Data/L1/MAPS-main-', MAPS_DATE, '.rds'))

MAPS_usp <- unique(mdata$sci_name)
MAPS_ucn <- unique(mdata$common_name)

#load taxonomy info - FROM BIRDNET
bt_phylo <- read.csv(paste0(dir, 'Data/L0/bird_phylo/BLIOCPhyloMasterTax.csv'))


# change names to match birdtree.org --------------------------------------

#where matches are
nin <- which(MAPS_usp %in% bt_phylo$TipLabel)

#names df
names_df <- data.frame(MAPS_sci_name = MAPS_usp, BT_sci_name = NA)
names_df$BT_sci_name[nin] <- MAPS_usp[nin]

#no matches
nm <- dplyr::filter(names_df, is.na(BT_sci_name))

#data.frame of matches
match <- data.frame(MAPS_sci_name = nm$MAPS_sci_name,
                    BT_sci_name = 
                      c('Carduelis_flammea',
                        'Wilsonia_canadensis',
                        'Wilsonia_pusilla',
                        'Regulus_calendula',
                        'Picoides_nuttallii',
                        'Picoides_pubescens',
                        'Picoides_villosus',
                        'Oporornis_formosus',
                        'Oporornis_philadelphia',
                        'Oporornis_tolmiei',
                        'Carpodacus_cassinii',
                        'Carpodacus_mexicanus',
                        'Carpodacus_purpureus',
                        'Zoothera_naevia',
                        'Vermivora_celata',
                        'Vermivora_luciae',
                        'Vermivora_peregrina',
                        'Vermivora_ruficapilla',
                        'Vermivora_virginiae',
                        'Pipilo_crissalis',
                        'Seiurus_motacilla',
                        'Seiurus_noveboracensis',
                        'Parus_atricapillus',
                        'Parus_carolinensis',
                        'Parus_gambeli',
                        'Parus_rufescens',
                        'Parula_americana',
                        'Dendroica_caerulescens',
                        'Wilsonia_citrina',
                        'Dendroica_coronata',
                        'Dendroica_discolor',
                        'Dendroica_magnolia',
                        'Dendroica_nigrescens',
                        'Dendroica_occidentalis',
                        'Dendroica_pensylvanica',
                        'Dendroica_petechia',
                        'Dendroica_townsendi',
                        'Dendroica_virens',
                        'Carduelis_pinus',
                        'Carduelis_psaltria',
                        'Troglodytes_troglodytes',
                        'Vermivora_pinus'))

#fill in gaps
for (i in 1:NROW(match))
{
  #i <- 1
  names_df$BT_sci_name[which(names_df$MAPS_sci_name == match$MAPS_sci_name[i])] <- match$BT_sci_name[i]
}


sci_names_phylo <- dplyr::left_join(names_df, 
                                    bt_phylo, 
                                    by = c('BT_sci_name' = 'TipLabel')) %>%
  dplyr::select(MAPS_sci_name, BT_sci_name, 
                common_name = English, 
                family_latin = BLFamilyLatin, 
                family_english = BLFamilyEnglish, 
                order = IOCOrder)


# write out ---------------------------------------------------------------

write.csv(sci_names_phylo, file = paste0(dir, 'Data/L1/bird_phylo/phylo_names-', 
                                     RUN_DATE, '.csv'), row.names = FALSE)


# tabulate ----------------------------------------------------------------

#number of orders
dplyr::count(sci_names_phylo, order)

#number of families
dplyr::count(sci_names_phylo, family_latin)

