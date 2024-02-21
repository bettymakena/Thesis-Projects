#Betty Makena
# 11-06-2023- Monday

library(raster)
library(sf)
library(ggplot2)
library(tidyverse)
library(rgdal)
library(dplyr)

setwd("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/CHIRPS/chirps_test")

chirps_mam <- list.files(pattern = ".tif")

# stack the raster files
chirps_mam_stk <- stack(chirps_mam)


#create an empty z-score stack
ch_null_crt <- stack()

for (i in 1:length(chirps_mam_stk)) {
  chirp_lyr <- chirps_mam_stk[[i]] #extract the nth layer
  chirp_lyr[chirp_lyr < 0] <- 0
  ch_null_crt <- stack(ch_null_crt, chirp_lyr)
  
}
#output
output <- ("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/CHIRPS/chirps_test/output")

for (i in 1:nlayers(ch_null_crt )) {
  # Extract the ith layer from the stack
  layer <- ch_null_crt[[i]]
  
  # Define the output file name for this layer (you can customize the naming)
  output_filename <- paste("Chirp_mam_", i, ".tif", sep = "")
  
  # Create the full output path
  output_path <- file.path(output, output_filename)
  
  # Save the layer as a TIFF file
  writeRaster(layer, filename = output_path, format = "GTiff")
}


