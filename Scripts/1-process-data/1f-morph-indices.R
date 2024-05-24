####################
# 1f - calc morph indices
#
####################


# set dirs ----------------------------------------------------------------

RUN_DATE <- '2024-05-24'
PHYLO_DATE <- '2024-05-24'
MAPS_DATE <- '2024-05-24'
dir <- '~/Work/Research/Projects/maps-morph-females/'


# load packages -----------------------------------------------------------

library(tidyverse)
library(ape)
library(caper)


# read in data ------------------------------------------------------------

mdata <- readRDS(paste0(dir, 'Data/L1/MAPS-main-', MAPS_DATE, '.rds'))


# filter by species -------------------------------------------------------

#only species with >= 375 data points adult males
mdata2 <- dplyr::group_by(mdata, sci_name) %>%
  dplyr::filter(sex == 'M',
                n() >= 375) %>%
  dplyr::ungroup() %>%
  #arrange and add sp_id
  dplyr::arrange(sci_name, station, year) %>%
  dplyr::mutate(sp_id = as.numeric(factor(sci_name)))


# estimate scaling exponent ----------------------------------------------------

# calculate species-level log(mass) and log(wing length)
mn_morph <- dplyr::group_by(mdata2, sci_name) %>%
  dplyr::summarize(mean_l_mass = mean(log(mass)),
                   mean_l_wing_chord = mean(log(wing_chord))) %>%
  dplyr::arrange(desc(mean_l_mass))

#relationship between wing length and body size follows a power law; wing length = mass^s. 
#Log both wing length and mass to linearize the relationship and estimate slope
#data from birdtree.org

#load in full set of trees for BirdTree.org
ptree <- ape::read.tree(paste0(dir, 'Data/L0/bird_phylo/birdtree/data/projects/birdphylo/Tree_sets/Stage2_full_data/CombinedTrees/AllBirdsEricson1.tre'))


#####
#from IBEEM
#species not found in both datasets (species to drop from tree)
nm <- setdiff(bird.phylo$tip.label, bird_df2$species)

#prune specified tips from tree
pr_tree <- ape::drop.tip(bird.phylo, nm)
#####


pnk <- read.csv(paste0(dir, 'Data/L1/bird_phylo/phylo_names-', 
                       PHYLO_DATE, '.csv'))

mn_morph2 <- dplyr::left_join(mn_morph, pnk, 
                              by = c('sci_name' = 'MAPS_sci_name')) %>%
  dplyr::arrange(sci_name) %>%
  as.data.frame()

# run phylogenetic regression for each tree realization
phy_reg_int <- rep(NA, length(ptree))
phy_reg_sl <- rep(NA, length(ptree))
for (i in 1:length(ptree))
{
  #i <- 1
  print(paste0('tree ', i, ' of ', length(ptree)))
  tree_n <- ptree[[i]]
  
  #species not found in both datasets (species to drop from tree)
  nm <- setdiff(tree_n$tip.label, mn_morph2$BT_sci_name)
  
  #prune specified tips from tree
  pr_tree <- ape::drop.tip(tree_n, nm)
  
  td <- caper::comparative.data(pr_tree, mn_morph2, BT_sci_name, 
                                vcv = TRUE, vcv.dim = 3)
  tf <- caper::pgls(mean_l_wing_chord ~ mean_l_mass, td)
  phy_reg_int[i] <- as.numeric(tf$model$coef[1])
  phy_reg_sl[i] <- as.numeric(tf$model$coef[2])
}

#mean intercept of regression from all trees
int <- mean(phy_reg_int)
#mean slope of regressions from all trees
sl <- mean(phy_reg_sl)
sd_sl <- sd(phy_reg_sl)


# reproject data points ---------------------------------------------------

# Rotate so that x-axis is the size of the bird (size index), and y axis is the wingy-ness of the bird (wingy-ness index).

#rotate points by slope of log(wl) log(mass) relationship
#get angle to rotate points (in radians)
angle <- -atan(sl)

#create rotation matrix
#R = \begin{bmatrix} cos(\theta) & -sin(\theta) \\ sin(\theta) & cos(\theta) \end{bmatrix}
#\begin{bmatrix} x' \\ y' \end{bmatrix} = \begin{bmatrix} cos(\theta) & -sin(\theta) \\ sin(\theta) & cos(\theta) \end{bmatrix} \begin{bmatrix} x \\ y \end{bmatrix}
#rotm = [[cos(theta) -sin(theta]
#        [sin(theta)  cos(theta)]]
#90 degrees
# rotm <- matrix(c(0, -1, 
#                  1, 0), ncol = 2)
rotm <- matrix(c(cos(angle), sin(angle), 
                 -sin(angle), cos(angle)), ncol = 2)

#function to log, rotate, then scale data
# std = TRUE -> standardized
# std = FALSE -> percent deviation from mean
# percent deviation suffers from not being so interpretatble bc morph metrics were logged - could use though
rotate_fun <- function(input)
{
  dm <- cbind(log(input$mass), log(input$wing_chord))
  M <- t(rotm %*% t(dm))
  
  mn_x <- mean(M[,1])
  sd_x <- sd(M[,1])
  mn_y <- mean(M[,2])
  sd_y <- sd(M[,2])
  
  M_sc <- apply(M, 2, function(x) scale(x, scale = TRUE))
  
  #mean and sd for each rotated axis
  ol <- list(M_sc = M_sc,
             mn_x = mn_x,
             sd_x = sd_x,
             mn_y = mn_y,
             sd_y = sd_y)
  
  return(ol)
}

#apply rotation matrix to data for each species individually and scale (center and std) to get 'size index' (x-axis of rotated data) and 'wing index' (y-axis of rotated data)
usp <- sort(unique(mdata2$sci_name))
mdata3 <- mdata2
mdata3$size_idx <- NA
mdata3$wingy_idx <- NA
mn_morph2$mn_x <- NA
mn_morph2$sd_x <- NA
mn_morph2$mn_y <- NA
mn_morph2$sd_y <- NA
for (i in 1:length(usp))
{
  #candidates: 37
  #i <- 98
  idx <- which(mdata3$sci_name == usp[i])
  temp <- mdata3[idx, ]
  
  #rotate
  rpd <- rotate_fun(temp)
  
  mdata3$size_idx[idx] <- rpd$M_sc[,1]
  mdata3$wingy_idx[idx] <- rpd$M_sc[,2]
  
  mn_idx <- which(mn_morph2$sci_name == usp[i])
  mn_morph2$mn_x[mn_idx] <- rpd$mn_x
  mn_morph2$sd_x[mn_idx] <- rpd$sd_x
  mn_morph2$mn_y[mn_idx] <- rpd$mn_y
  mn_morph2$sd_y[mn_idx] <- rpd$sd_y
}

# apply rotation matrix to species-level means
M_sp <- t(rotm %*% t(cbind(mn_morph2$mean_l_mass, mn_morph2$mean_l_wing_chord)))
M_sp_sc <- apply(M_sp, 2, function(x) scale(x, scale = TRUE))
mn_morph2$sp_size_idx <- M_sp_sc[,1]
mn_morph2$sp_wingy_idx <- M_sp_sc[,2]


# save rds ----------------------------------------------------------------

saveRDS(mdata3, paste0(dir, 'Data/L2/MAPS-processed-', RUN_DATE, '.rds'))

