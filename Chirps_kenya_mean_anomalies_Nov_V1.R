library(raster)
library(sf)
library(ggplot2)
library(tidyverse)
library(rgdal)
library(dplyr)

# MAM analysis ------------------------------------------------------------
#setwd("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/INPUTS/CHIRPS/OND")
setwd("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/data_2024/Chapter 2/Long  term means/OND/Standardized_Aug_Nov/Augu_sept/RAW")

#output <- ("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/data_2024/chapter1/chirps/OND/Mean_chirps")
output <- ("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/data_2024/Chapter 2/Long  term means/OND/Standardized_Aug_Nov/Augu_sept/mean")


# load kenya shapefile
kenya <- st_read("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/admin/kenya_edited.shp")

#kenya2 <- kenya$geometry

# List all the raster files in the directory with the specified pattern
MAM_rast <- ras_sep <- list.files(pattern = ".tif$")
  
#list.files(pattern = "chirps-v2.0\\.\\d{4}\\.\\d{2}\\.tif") 

# Extract unique years from file names
unique_years <- unique(gsub(".*\\.(\\d{4})\\.\\d{2}\\.tif", "\\1", MAM_rast))

# Create an empty list to store the raster objects for each year
raster_list <- list()


# Loop through unique years, read files, calculate mean, and clip for each year
for (year in unique_years) {
  # Filter files for the current year
  files_for_year <- grep(paste0(year, "\\."), MAM_rast, value = TRUE)
  
  # Read raster files for the current year
  rasters_for_year <- lapply(files_for_year, raster)
  
  # Calculate mean for the rasters of the current year
  mean_raster <- mean(stack(rasters_for_year))
  
  # Clip the mean raster to a specific extent or region
  # Replace the extent coordinates with your desired values
  
  #mean_raster_clipped <- crop(mean_raster, kenya)
  
  # Add clipped mean raster to the list with the year as the identifier
  raster_list[[year]] <- mean_raster
}
#save the mean
writeRaster(kch_mean, file.path(output, "mean_chirps_K_mam2.tif"), format = "GTiff")

#calculate the standard deviation
kch_sd <- calc(kmam_stack, fun = sd, na.rm = TRUE)

#save the SD
writeRaster(kch_sd, file.path(output, "SD_chirps_K_mam2.tif"), format = "GTiff")

#write to folder
Czcore_mam <- stack()

#calculate the zscore
for (i in 1:length(kmam_stack)) {
  chirp_layer <- kmam_stack[[i]] #extract the nth layer
  z_score <- (chirp_layer - kch_mean)/kch_sd
  Czcore_mam <- stack(Czcore_mam, z_score)
  
}

#save the layers
for (i in 1:length(Czcore_mam)) {
  # Extract the ith layer from the stack
  layer <- Czcore_mam[[i]]
  
  # Define the output file name for this layer (you can customize the naming)
  output_filename <- paste("Crp_Zscore_K_", i, ".tif", sep = "")
  
  # Create the full output path
  output_path <- file.path(output, output_filename)
  
  # Save the layer as a TIFF file
  writeRaster(layer, filename = output_path, format = "GTiff")
} 




# OND Analysis ------------------------------------------------------------

setwd("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/INPUTS/CHIRPS/OND")
output <- ("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/Outputs/CHIRPS/OND/kenya_OND_mean")
# load kenya shapefile
kenya <- st_read("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/admin/kenya_edited.shp")

#kenya2 <- kenya$geometry

# List all the raster files in the directory with the specified pattern
ond_rast <- list.files(pattern = "chirps-v2.0\\.\\d{4}\\.\\d{2}\\.tif") 

# Extract unique years from file names
unique_years <- unique(gsub(".*\\.(\\d{4})\\.\\d{2}\\.tif", "\\1", ond_rast))

# Create an empty list to store the raster objects for each year
raster_list <- list()


# Loop through unique years, read files, calculate mean, and clip for each year
for (year in unique_years) {
  # Filter files for the current year
  files_for_year <- grep(paste0(year, "\\."), ond_rast, value = TRUE)
  
  # Read raster files for the current year
  rasters_for_year <- lapply(files_for_year, raster)
  
  # Calculate mean for the rasters of the current year
  mean_raster <- mean(stack(rasters_for_year))
  
  # Clip the mean raster to a specific extent or region
  # Replace the extent coordinates with your desired values
  
  #mean_raster_clipped <- crop(mean_raster, kenya)
  
  # Add clipped mean raster to the list with the year as the identifier
  raster_list[[year]] <- mean_raster
}

# save the raster layers
#save the layers
for (i in 1:length( raster_list)) {
  # Extract the ith layer from the stack
  layer <- raster_list[[i]]
  
  # Define the output file name for this layer (you can customize the naming)
  output_filename <- paste("Chp_mean_ond_K", i, ".tif", sep = "")
  
  # Create the full output path
  output_path <- file.path(output, output_filename)
  
  # Save the layer as a TIFF file
  writeRaster(layer, filename = output_path, format = "GTiff")
} 

# Calculating anomalies ---------------------------------------------------


#load the clipped layers
setwd("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/Outputs/CHIRPS/MAM/Chapter 2/Chirps_MAM/kenya clipped")

#output folder
output <- ("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/Outputs/CHIRPS/MAM/Chapter 2/Chirps_MAM/zscore_mam")

raster_files <- list.files(pattern = ".tif$")

#stack the images
kmam_stack <- stack(raster_files)

#calculate the mean of the raster files
kch_mean <- calc(kmam_stack, fun = mean, na.rm = TRUE)

#save the mean
writeRaster(kch_mean, file.path(output, "mean_chirps_K_mam.tif"), format = "GTiff")

#calculate the standard deviation
kch_sd <- calc(kmam_stack, fun = sd, na.rm = TRUE)

#save the SD
writeRaster(kch_sd, file.path(output, "SD_chirps_K_mam.tif"), format = "GTiff")

#write to folder
Czcore_mam <- stack()

#calculate the zscore
for (i in 1:length(kmam_stack)) {
  chirp_layer <- kmam_stack[[i]] #extract the nth layer
  z_score <- (chirp_layer - kch_mean)/kch_sd
  Czcore_mam <- stack(Czcore_mam, z_score)
  
}

#save the layers
for (i in 1:length(Czcore_mam)) {
  # Extract the ith layer from the stack
  layer <- Czcore_mam[[i]]
  
  # Define the output file name for this layer (you can customize the naming)
  output_filename <- paste("Crp_Zscore_K_", i, ".tif", sep = "")
  
  # Create the full output path
  output_path <- file.path(output, output_filename)
  
  # Save the layer as a TIFF file
  writeRaster(layer, filename = output_path, format = "GTiff")
} 



