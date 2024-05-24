# maps-morph-females

Adapted for ongoing work.

Code originally for Youngflesh et al. 2022 _Nature Ecology and Evolution_

This repository contains code to characterize avian morphological change over time, latitude, elevation, and in response to temperature.

**Associated publications:**

Youngflesh, C, JF Saracco, RB Siegel, MW Tingley. 2022. [Abiotic conditions shape spatial and temporal morphological variation in North American birds.](https://www.nature.com/articles/s41559-022-01893-x) **_Nature Ecology and Evolution_** 6:1860-1870.


**Repository structure:**

* `Scripts/`
  * `1-process-data/`
    * `1a-process-elev.R` - stitch elevation data into single .tif
    * `1b-process-MAPS.R` - process MAPS data to create main MAPS data file
    * `1c-process-range-maps.R` - process range maps for each species
    * `1d-process-daymet.R` - DL and process daymet (temp, precip) data
    * `1e-process-phylo-names.R` - process phylo names (matching)
    * `1f-morph-indices.R` - calculate morph indices
  - `4-si-tle.R` - old
  - `5-wi-tle.R` - old
  - `6-si-temp-space.R` - old
  - `7-si-temp-lag.R` - old
  - `8-temp-cov.R` - old
  - `9-backtransform.R` - old
  - `10-analyze-results.R` - old
  - `11-haldanes.R` - old
  - `model_files/`
    - `morph-temp-space.stan` - old
    - `morph-temp-ss2.stan` - old
    - `morph-tle.stan` - old
    - `temp_sp_cov.stan` - old
* `Data/`
  * `LO/` - Raw data
    * `bird_phylo/` - phylogenetic tree data frmo birdtree.org
    * `BOTW_2022/` - birdlife range maps
    * `daymet_query_locations.csv` - maps station locations for daymet query
    * `daymet-dl-YYYY-MM-DD.rds` - daymet downloaded (pre-processing)
    * `GMTED2010_mean/` - elevation data
    * `MAPS_data/` - raw MAPS data
  * `L1/` - Level 1 processing
    * `bird_phylo` - processed bir phylogenetic names
    * `BOTW_ranges-YYYY-MM-DD.rds` - subset of bird ranges
    * `daymet-main-YYYY-MM-DD.rds` - processed daymet data (monthly)
    * `elev_mosaic.tif` - mosaiced elevation data
    * `MAPS-main-YYYY-MM-DD.rds` - processed MAPS data
    * `range_names_key-YYYY-MM-DD.csv` - range names for processing
  * `L2/` - Level 2 processing
    * `MAPS-processed-YYYY-MM-DD.rds` - MAPS data with morph indices (males only as of 2024-05-24)