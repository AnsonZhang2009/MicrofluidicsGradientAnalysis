#=============================================
#
# Prerequisites 
#
#=============================================

# If not already installed, install the image processing package "magick". 
# install.packages("magick")

library(magick)
library(dplyr)
library(ggplot2)

#=============================================
#
# Read Image
#
#=============================================

img <- image_read("Images/Untitled.png")

# Use the quantize function to merge pixels. 
palette <- image_quantize(img, max = 50, colorspace = "RGB")

print(palette)

info <- image_info(palette)

data <- image_data(palette, channels = "RGB")

# Reshape to get each pixel as a row with RGBA values. 
total_pixels <- info$width * info$height
pixel_matrix <- matrix(data, nrow = 3, ncol = total_pixels)
# Transpose matrix to show pixels as rows. 
pixel_matrix <- t(pixel_matrix)  

# This will count the number of unique combinations for the quantized data. 
unique_colors <- unique(pixel_matrix)
num_colors <- nrow(unique_colors)
print(num_colors)


