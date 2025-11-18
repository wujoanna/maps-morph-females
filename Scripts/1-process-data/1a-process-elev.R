#######################
# 1a - mosaic elevation data
#
# GMTED2010 Global Elevation Product - 30 arc seconds (~ 1km)
#https://www.usgs.gov/centers/eros/science/usgs-eros-archive-digital-elevation-global-multi-resolution-terrain-elevation?qt-science_center_objects=0#qt-science_center_objects
# Region specified on Earth Explorer and downloaded using the EE bulk download tool
# Only mean elevation for each pixel retained to save space (each metric [min, max, std, ...] was downloaded as a different .tiff)
# *** GDAL must be installed on computer ***
#######################


# set dirs -------------------------------------------------------------------

# CY
# dir <- '~/Work/Research/Projects/maps-morph-females/'
# JW
dir <- "/Users/joannawu/Library/CloudStorage/Box-Box/Documents/AcademicScience/Projects/MorphFemales/"
setwd(dir)

# mosaic GMTED tifs -------------------------------------------------------

#mosaic GMTED2010 DEM tiles together with gdal
#call gdal from command line through R
system(paste0('gdalbuildvrt ', 
              dir, 'Data/L0/GMTED2010_mean/elev_mosaic.vrt ', 
              dir, 'Data/L0/GMTED2010_mean/*.tif'))
system(paste0('gdal_translate ', 
              dir, 'Data/L0/GMTED2010_mean/elev_mosaic.vrt ',
              dir, 'Data/L1/elev_mosaic.tif'))
