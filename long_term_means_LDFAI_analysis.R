#Betty Makena
# 11-26-2023- Sunday
library(raster)
library(sf)
library(ggplot2)
library(tidyverse)
library(rgdal)
library(dplyr)



# Reading files -----------------------------------------------------------
setwd("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/data_2024/Chapter 2/Long  term means/OND/Stand_Nov_Feb/Jan_feb/Means")
output <- ("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/data_2024/Chapter 2/Long  term means/OND/Stand_Nov_Feb/Jan_feb/Rstan_Jan_feb_83_23")


#selected_files <- raster_files[grep("\\.(03|04|05)\\.tif$", raster_files)]
ras_sep <- list.files(pattern = ".tif$")
 

sep_stack <- stack(ras_sep)

#average means for septs
#mean_sept
sept_avg_mean <- calc(sep_stack, fun = mean, na.rm = TRUE)

#save the mean
writeRaster(sept_avg_mean, file.path(output, "Avg_mean_sept.tif"), format = "GTiff")

diff_stck <- stack()

#Calculating the sum of squares
for (i in 1:length(sep_stack)){
  mnth <- sep_stack[[i]] #extract each month
  difn <- (mnth - sept_avg_mean) ^ 2
  diff_stck <- stack(diff_stck,difn ) 
}

#saving the differences files 
for (i in 1:nlayers(diff_stck)) {
  # Extract the ith layer from the stack
  layer <- diff_stck[[i]]
  
  # Define the output file name for this layer (you can customize the naming)
  output_filename <- paste("Dif_MAM_sept_", i, ".tif", sep = "")
  
  # Create the full output path
  output_path <- file.path(output, output_filename)
  
  # Save the layer as a TIFF file
  writeRaster(layer, filename = output_path, format = "GTiff")
}


# sum up all the layers in the difference folder
sum_diff <- calc(diff_stck, fun = sum, na.rm = TRUE)

#save the sum
writeRaster(sum_diff, file.path(output, "sum_diff_sept.tif"), format = "GTiff")

n = 1/42

# final formulae for the root of squares n=43
rs_sept <- sqrt(n* sum_diff)


#save the sum
writeRaster(rs_sept, file.path(output, "sum_square_sept.tif"), format = "GTiff")

stan_prep <- stack()

#Standardizing the precipitation( monthl prep( in sep_stack), ltm_ave_prep(sept_avg_mean), sum of squares(rs_sept)
for (i in 1:length(sep_stack)){
  mnth <- sep_stack[[i]] #extract each month
  R_prep <- (mnth - sept_avg_mean)/rs_sept
  stan_prep <- stack(stan_prep,R_prep)
  
}

# Save the final standardized layers(R)

#saving the differences files 
for (i in 1:nlayers(stan_prep)) {
  # Extract the ith layer from the stack
  layer <- stan_prep[[i]]
  
  # Define the output file name for this layer (you can customize the naming)
  output_filename <- paste("Rsrd_jan_feb_v2_", i, ".tif", sep = "")
  
  # Create the full output path
  output_path <- file.path(output, output_filename)
  
  # Save the layer as a TIFF file
  writeRaster(layer, filename = output_path, format = "GTiff")
}



# Calculating the LDFAI final Formulae ------------------------------------
#02-03-2024
#LDFAI= (Rnd - RSO)*(|Rnd| + |RSO|) * 1.8^(|Rnd| + |RSO|)
 
sep_oct <- ("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/data_2024/Chapter 2/Long  term means/OND/Standardized/Sept_Oct/Final_R")
lsep_O <- list.files(sep_oct,pattern = ".tif$")

nov_dec <- ("C:/Users/mbetty2/OneDrive - University of Nebraska-Lincoln/Projects/data_2024/Chapter 2/Long  term means/OND/Standardized/Nov_Dec/Final_R")
lno_de<- list.files(nov_dec,pattern = ".tif$")


# Difference (Rnd - RSO)
# Extract unique years from file names
years <- as.numeric(gsub("Rsrd_Sept_Oct_(\\d+)\\.tif", "\\1", lsep_O,lno_de))

diff <- stack()
ab_sum <- stack()

for (year in years){
  file_ND <- grep(paste0(year, "\\."), lno_de, value = TRUE)
  file_SO <- grep(paste0(year, "\\."), lsep_O, value = TRUE)
  ras_nd <- stack(lapply(file_ND, raster))
  ras_SO <- stack(lapply(file_SO, raster))
  ld_min <- (ras_nd - ras_SO)
  ld_sum <- abs(ras_nd) + abs(ras_SO)
  diff <- stack(diff, ld_min[[1]])
  ab_sum <- stack(ab_sum, ld_sum[[1]])
}


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

sep_stack <- stack(september)
dec <- stack(marriage)







# september ---------------------------------------------------------------

#september <- raster_files[str_detect(raster_files, "^chirps-v2\\.0\\.\\d{4}\\.09\\.tif$")]

sep_stack <- stack(september)

#mean_sept
sept_mean <- calc(sep_stack, fun = mean, na.rm = TRUE)

#save the mean
writeRaster(sept_mean, file.path(output, "Lmean_Chp_sept.tif"), format = "GTiff")

#calculate the standard deviation
sept_sd <- calc(sep_stack, fun = sd, na.rm = TRUE)

#save the mean
writeRaster(sept_sd, file.path(output, "LSD_Chp_sept.tif"), format = "GTiff")


# October -----------------------------------------------------------------
october <- raster_files[str_detect(raster_files, "^chirps-v2\\.0\\.\\d{4}\\.10\\.tif$")]

sep_stack <- stack(september)

#mean_sept
sept_mean <- calc(sep_stack, fun = mean, na.rm = TRUE)

#save the mean
writeRaster(sept_mean, file.path(output, "Lmean_Chp_sept.tif"), format = "GTiff")

#calculate the standard deviation
sept_sd <- calc(sep_stack, fun = sd, na.rm = TRUE)

#save the mean
writeRaster(sept_sd, file.path(output, "LSD_Chp_sept.tif"), format = "GTiff")







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







