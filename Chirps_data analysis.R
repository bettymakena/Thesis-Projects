#Betty Makena
# 10-27-2023- Monday
library(raster)
library(sf)
library(ggplot2)
library(tidyverse)
library(rgdal)
#library(dplyr)

# Reading files -----------------------------------------------------------

setwd("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/Data/CHIRPS/zipped")

#selected_files <- raster_files[grep("\\.(03|04|05)\\.tif$", raster_files)]
raster_files <- list.files(pattern = ".tif")

# Filter file names to select those with "04" in the specified position
march <- raster_files[str_detect(raster_files, "^chirps-v2\\.0\\.\\d{4}\\.03\\.tif$")]

april <- raster_files[str_detect(raster_files, "^chirps-v2\\.0\\.\\d{4}\\.04\\.tif$")]

may <- raster_files[str_detect(raster_files, "^chirps-v2\\.0\\.\\d{4}\\.05\\.tif$")]

#join the three lists
mam <- c(march, april, may)

#first approach- stack, extact unique years, calculate mean based on unique years

#stack the  raster files

mam_stack <- stack(mam)

#create an empty mean stack
mean_chirps <- stack()


#loop through the files 3,4,5 for every year

for (i in 1:length(mam_stack)) {
  chirps_layer <- mam_stack[[i]] #extract the nth layer
  z_score <- (ndvi_layer - mean_ndvi)/sd_ndvi
  ndvi_mam_Z <- stack(ndvi_mam_Z, z_score)
  
}

for files 

# join  the three list
# group by subtring(files, )



selected <- substring(raster_files, 17,20) == "04"

raster_sfiles2 <- list.files(str_detect(raster_files, "^chirps-v2\\.0\\.\\d{4}\\.04\\.tif$"))

raster_files <- list.files(pattern = c(".tiff"))

raster_files <- list.files(substring())

files = list.files(pattern = '.*04.*\\.tif$')

student_files <- list.files(pattern = 'exams_*.csv')
#the months
selected_months <- c("03", "04", "05")

p <- "bettymakena"

for (i in data_dir) {
  
  
}

for (i in raster_files) {
  raster_layer <- raster(file_path)
  raster_layers[[length(raster_layers) + 1]] <- raster_layer
}

# List of years to process
years <- 1990:2022

#the months
selected_months <- c(3, 4, 5)

# Create an empty list to store the mean rasters
mean_rasters <- list()


for (year in years) {
  # Create a list to store the monthly rasters for this year
  monthly_rasters <- list()
  
  for (month in selected_months) {
    # Create the file name based on the year and month
    file_name <- paste("chirps-v2.0.", year, ".", sprintf("%02d", month), ".tif", sep = "")
    file_path <- file.path(data_dir, file_name)
    
    if (file.exists(file_path)) {
      # Read the monthly raster file
      monthly_raster <- raster(file_path)
      monthly_rasters[[month - 2]] <- monthly_raster
    }
  }
  
  # Stack the monthly rasters for this year
  year_stack <- stack(monthly_rasters)
  
  # Calculate the mean for these months
  mean_raster <- mean(year_stack, na.rm = TRUE)
  
  # Add the result to the list
  mean_rasters[[year - 1989]] <- mean_raster
}







