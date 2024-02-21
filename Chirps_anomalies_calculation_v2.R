#Betty Makena
# 11-08-2023
 
# Outline -----------------------------------------------------------------
#the code calculates the mean,sd and anomalies and saves  of chirps data 1982-2022
#loa libraries
library(raster)
library(sf)
library(ggplot2)
library(tidyverse)
library(rgdal)
library(dplyr)


# MAM session -------------------------------------------------------------
setwd("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/Outputs/CHIRPS/MAM/null_corrected")



raster_files <- list.files(pattern = ".tif$")
output <-("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/Outputs/CHIRPS/MAM/z-score_v2") 

#stack the images
mam_stack <- stack(raster_files)

#calculate the mean of the raster files
ch_mean <- calc(mam_stack, fun = mean, na.rm = TRUE)

#save the mean
writeRaster(ch_mean, file.path(output, "mean_chirps_mam.tif"), format = "GTiff")

#calculate the standard deviation
ch_sd <- calc(mam_stack, fun = sd, na.rm = TRUE)

#save the mean
writeRaster(ch_sd, file.path(output, "SD_chirps_mam.tif"), format = "GTiff")

#read the mean and SD. corrected in arcgis to remove negative values
ch_sd_mam <- raster("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/Outputs/CHIRPS/MAM/z-score/ch_sd_mam.tif")
ch_mean_mam <- raster("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/Outputs/CHIRPS/MAM/z-score/ch_mean_mam_v2.tif")

Czcore_mam <- stack()

#calculate the zscore
for (i in 1:length(mam_stack)) {
 chirp_layer <- mam_stack[[i]] #extract the nth layer
  z_score <- (chirp_layer - ch_mean)/ch_sd
  Czcore_mam <- stack(Czcore_mam, z_score)
  
}
 #save the files

for (i in 1:nlayers(Czcore_mam)) {
  # Extract the ith layer from the stack
  layer <- Czcore_mam[[i]]
  
  # Define the output file name for this layer (you can customize the naming)
  output_filename <- paste("Chirps_MAM_", i, ".tif", sep = "")
  
  # Create the full output path
  output_path <- file.path(output, output_filename)
  
  # Save the layer as a TIFF file
  writeRaster(layer, filename = output_path, format = "GTiff")
}


# ONDA Analysis -----------------------------------------------------------
setwd("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/Outputs/CHIRPS/OND/null_corrected")


raster_files <- list.files(pattern = ".tif$")
output <- ("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/Outputs/CHIRPS/OND/zscore_v2")
           
#stack the images
ond_stack <- stack(raster_files)

#calculate the mean of the raster files
ch_mean_o <- calc(ond_stack, fun = mean, na.rm = TRUE)

#save the mean
writeRaster(ch_mean, file.path(output, "mean_chirps_ond.tif"), format = "GTiff")

#calculate the standard deviation
ch_sd_o <- calc(ond_stack, fun = sd, na.rm = TRUE)

#save the mean
writeRaster(ch_sd, file.path(output, "SD_chirps_ond.tif"), format = "GTiff")

Czcore_ond <- stack()

#calculate the zscore
for (i in 1:length(ond_stack)) {
  chirp_layer <- ond_stack[[i]] #extract the nth layer
  z_score <- (chirp_layer - ch_mean_o)/ch_sd_o
  Czcore_ond <- stack(Czcore_ond, z_score)
  
}

#save the layers
for (i in 1:nlayers(Czcore_ond)) {
  # Extract the ith layer from the stack
  layer <- Czcore_ond[[i]]
  
  # Define the output file name for this layer (you can customize the naming)
  output_filename <- paste("Chirps_OND_", i, ".tif", sep = "")
  
  # Create the full output path
  output_path <- file.path(output, output_filename)
  
  # Save the layer as a TIFF file
  writeRaster(layer, filename = output_path, format = "GTiff")
}

