#Betty Makena
# 10-23-2023- Monday

library(raster)
library(sf)
library(ggplot2)
library(tidyverse)
library(rgdal)
library(dplyr)

data_dir <- ("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/NDVI/Raster_mam_tiff")


# Anomaly for MAM season --------------------------------------------------

#list all the files
raster_files <- list.files(pattern = ".tif")
ras_names <- names(raster_files)


for (i in raster_files) {
  raster_layer <- raster(file_path)
  raster_layers[[length(raster_layers) + 1]] <- raster_layer
}

# stack the raster files
stack_raster <- stack(raster_files)
stack_names <- names(stack_raster)

#calculate mean for the raster stak
#mean_ndvi <- mean(stack_raster, na.rm = TRUE)
#read the raster
mean_ndvi <- raster("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/Outputs/NDVI/MAM/mean_ndvi_mam.tif")
sd_ndvi <-  raster("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/Outputs/NDVI/MAM/sd_ndvi_mam.tif")

#calc allows calculation to each layer in the stack
#sd_ndvi <- calc(stack_raster, fun = sd, na.rm = TRUE)
#sd_ndvi <- sd(stack_raster, na.rm = TRUE)

#output folder
output <- ("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/Outputs/NDVI/MAM/zscore")

#save the mean
writeRaster(mean_ndvi, file.path(output, "mean_ndvi_mam.tif"), format = "GTiff")

#save the sd
writeRaster(sd_ndvi, file.path(output, "sd_ndvi_mam.tif"), format = "GTiff")

#calculate the z-score(anomailes)

#create an empty z-score stack
ndvi_mam_Z <- stack()

for (i in 1:length(stack_raster)) {
  ndvi_layer <- stack_raster[[i]] #extract the nth layer
  z_score <- (ndvi_layer - mean_ndvi)/sd_ndvi
  ndvi_mam_Z <- stack(ndvi_mam_Z, z_score)
  
}

z <- brick(ndvi_mam_Z)

#save every layer in the raster stack

# Get the names of the layers in the raster stack
layer_names <- names(ndvi_mam_Z)

# Loop through the layers in the raster stack
for (i in 1:nlayers(ndvi_mam_Z)) {
  # Extract the ith layer from the stack
  layer <- ndvi_mam_Z[[i]]
  
  # Extract the name of the layer from the layer names
  layer_name <- layer_names[i]
  
  # Define the output file name using the layer name
  output_filename <- paste(layer_name, ".tif", sep = "")
  
  # Create the full output path
  output_path <- file.path(output, output_filename)
  
  # Save the layer as a TIFF file
  writeRaster(layer, filename = output_path, format = "GTiff")
}



for (i in 1:nlayers(ndvi_mam_Z)) {
  # Extract the ith layer from the stack
  layer <- ndvi_mam_Z[[i]]
  
  # Define the output file name for this layer (you can customize the naming)
  output_filename <- paste("MAM_OND_", i, ".tif", sep = "")
  
  # Create the full output path
  output_path <- file.path(output, output_filename)
  
  # Save the layer as a TIFF file
  writeRaster(layer, filename = output_path, format = "GTiff")
}

writeRaster(ndvi_mam_Z, filename = output, format = "GTiff")



# Anomally for OND season -------------------------------------------------

setwd("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/NDVI/NDVI_OND")

#list all the files
raster_files <- list.files(pattern = ".tif")

#stack the images
ndvi_stack <- stack(raster_files)

# read ndvi and sd files
mean_ndvi_on <- raster("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/Outputs/NDVI/OND/zscore/mean_ndvi_ond.tif")

sd_ndvi_on <- raster("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/Outputs/NDVI/OND/zscore/mean_sd_ond.tif")

output <- ("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/Outputs/NDVI/OND/zscore")

Nzcore_ond <- stack()

#calculate the zscore
for (i in 1:length(ndvi_stack )) {
  ndvi_layer <- ndvi_stack[[i]] #extract the nth layer
  z_score <- (ndvi_layer - mean_ndvi_on)/sd_ndvi_on
  Nzcore_ond <- stack(Nzcore_ond, z_score)
  
}

#save the layers
for (i in 1:nlayers(Nzcore_ond )) {
  # Extract the ith layer from the stack
  layer <- Nzcore_ond [[i]]
  
  # Define the output file name for this layer (you can customize the naming)
  output_filename <- paste("NDVI_OND_", i, ".tif", sep = "")

  # Create the full output path
  output_path <- file.path(output, output_filename)
  
  # Save the layer as a TIFF file
  writeRaster(layer, filename = output_path, format = "GTiff")
}




