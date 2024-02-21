library(raster)
library(sf)
library(ggplot2)
library(tidyverse)
library(rgdal)
library(dplyr)

# MAM analysis ------------------------------------------------------------

setwd("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/INPUTS/CHIRPS/MAM")
output <- ("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/Outputs/CHIRPS/MAM/Kenya_mam_mean_V2")

# load kenya shapefile
kenya <- st_read("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/admin/kenya_edited.shp")

#kenya2 <- kenya$geometry

# List all the raster files in the directory with the specified pattern
MAM_rast <- list.files(pattern = "chirps-v2.0\\.\\d{4}\\.\\d{2}\\.tif") 

# Extract unique years from file names
unique_years <- unique(gsub(".*\\.(\\d{4})\\.\\d{2}\\.tif", "\\1", MAM_rast))

# Create an empty list to store the raster objects for each year
raster_stk <- stack()


# Loop through unique years, read files, calculate mean, and clip for each year
for (year in unique_years) {
  # Filter files for the current year
  files_for_year <- grep(paste0(year, "\\."), MAM_rast, value = TRUE)
  
  # Read raster files for the current year
  rasters_for_year <- lapply(files_for_year, raster)
  
  # Calculate mean for the rasters of the current year
  mean_raster <- mean(stack(rasters_for_year))
  
  # Clip the mean raster 
  # Replace the extent coordinates with your desired values
  
  #mean_raster_clipped <- crop(mean_raster, kenya) # not working
  
  # Add clipped mean raster to the list with the year as the identifier
  raster_stk  <- stack(raster_stk,mean_raster)
  
}

# save the raster layers
#save the 
for (i in 1:length(raster_stk)) {
  # Extract the ith layer from the stack
  layer <- raster_stk[[i]]
  
  # Define the output file name for this layer (you can customize the naming)
  output_filename <- paste("Chirp_MAM_K", i, ".tif", sep = "")
  
  # Create the full output path
  output_path <- file.path(output, output_filename)
  
  # Save the layer as a TIFF file
  writeRaster(layer, filename = output_path, format = "GTiff")
} 
