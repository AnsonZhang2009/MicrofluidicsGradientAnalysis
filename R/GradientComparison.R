#=============================================
#
# Prerequisites 
#
#=============================================

library(magick)
library(ggplot2)
library(dplyr)
library(stringr)
options(stringsAsFactors = FALSE)

#=============================================
#
# Data Input Functions
#
#=============================================

# Use command line input. 
args <- commandArgs(trailingOnly = TRUE)

# If not used in RStudio, the above option would be more convenient. 
args <- c("Images/") 

if (length(args) == 0) {
  stop("Number of files cannot be 0. ")
}

inputFiles <- function(paths) {
  
  validFiles <- c()
  
  for (path in paths) {
    
    if (file.exists(path)) {
      if (dir.exists(path)) {
        
        imageFiles <- list.files(path,
                                 pattern = "\\.(jpg|jpeg|png|tiff)",
                                 ignore.case = TRUE,
                                 full.names = TRUE)
        
        validFiles <- c(validFiles, imageFiles)
        
      } else {
        warning(paste("No path found, path: ", path))
      }
    } else {
      warning(paste("No path found, path:", path))
    }
  }
  
  if (length(validFiles) == 0) {
    stop("No files found. Please provide an accurate directory. ")
  }
  
  return(validFiles)
}

#=============================================
#
# Normalize and Convert to Grayscale
#
#=============================================

normalizeImage <- function(img, targetHeight = 500, targetWidth = 500) {
  
  # Convert to grayscale
  imgGray <- image_convert(img, colorspace = "gray")
  
  # Obtain the current dimensions of the image
  info <- image_info(imgGray)
  width0 <- info$width
  height0 <- info$height
  
  # Scale height and width. We would prefer the scaling with less magnitude because that means that we would lose less data. We only need either the height or the width to meet our requirement.  
  scaleW <- targetWidth / width0
  scaleH <- targetHeight / height0
  scaleFactor <- min(scaleW, scaleH)
  
  width1 <- round(width0 * scaleFactor)
  height1 <- round(height0 * scaleFactor)
  
  # Resize the image
  imgResized <- image_resize(imgGray, geometry = paste0(width1, "x", height1, "!"))
  
  # We will need to normalize the image to ensure that all images are not affected by the difference in the lighting
  imgNormalized <- image_normalize(imgResized)
  
  return(imgNormalized)
}

#=============================================
#
# Sobel Edge Detection Function
#
#=============================================

sobelGradient <- function(imgPath) {
  
  tryCatch({
    
    # Image input
    img <- image_read(imgPath)
    
    # Normalize the brightness of the image (check Magick docs for how this is done or if any changes might need to be made)
    imgNormalized <- normalizeImage(img)
    
    # Obtain a bit map of the normalized photo and convert it to a numeric matrix
    imgMatrix <- as.numeric(image_data(imgNormalized, channels = "gray"))
    
    # Obtain the dimensions of the image
    info <- image_info(imgNormalized)
    width <- info$width
    height <- info$height
    
    # Since imgMatrix is a 1d array, we rearrange it into a 2d array similar to a bit map holding information regarding every pixel
    img2D <- matrix(imgMatrix, nrow = height, ncol = width, byrow = TRUE)
    
    # Sobel kernels (as mentioned here https://www.mdpi.com/2313-433X/4/6/74)
    sobelX <- matrix(c(-1, 0, 1,
                       -2, 0, 2,
                       -1, 0, 1), nrow = 3, byrow = TRUE)
    
    sobelY <- matrix(c(-1, -2, -1,
                       0,  0,  0,
                       1,  2,  1), nrow = 3, byrow = TRUE)
    
    # Apply convolution
    gradientX <- applySobelKernel(img2D, sobelX)
    gradientY <- applySobelKernel(img2D, sobelY)
    
    # Calculate the gradient magnitude
    gradientMagnitude <- sqrt(gradientX^2 + gradientY^2)
    
    # Calculate gradient orientation η = arctan(Iy/Ix)
    gradientOrientation <- atan2(gradientY, gradientX)
    
    # Calculate overall gradient strength (mean of magnitude)
    gradientStrength <- mean(gradientMagnitude, na.rm = TRUE)
    
    # Apply non-maximum suppression
    edges <- nonMaximumSuppression(gradientMagnitude, gradientOrientation)
    
    # Count edge pixels above threshold
    threshold <- quantile(gradientMagnitude, 0.85, na.rm = TRUE)
    edgePixelCount <- sum(edges > threshold, na.rm = TRUE)
    
    # Save gradient magnitude image
    gradientImg <- image_read(as.raster(gradientMagnitude/max(gradientMagnitude, na.rm = TRUE)))
    directory <- tools::file_path_sans_ext(imgPath)
    filename <- paste0(directory, "_gradient.png")
    image_write(gradientImg, path = filename)
    cat("Saved gradient image:", filename, "\n")
    
    return(list(
      gradientStrength = gradientStrength,
      edgePixelCount = edgePixelCount,
      maxGradient = max(gradientMagnitude, na.rm = TRUE)
    ))
    
  }, error = function(e) {
    warning(paste0("Error processing ", imgPath, ": ", e$message))
    return(list(gradientStrength = NA, edgePixelCount = NA, maxGradient = NA))
  })
}

#=============================================
#
# Sobel Kernel Transformation
#
#=============================================

applySobelKernel <- function(imgMatrix, kernel) {
  
  rows <- nrow(imgMatrix)
  cols <- ncol(imgMatrix)
  result <- matrix(0, nrow = rows, ncol = cols)
  
  # Apply kernel to each pixel (excluding borders)
  for (i in 2:(rows-1)) {
    for (j in 2:(cols-1)) {
      # Extract 3x3 neighborhood
      neighborhood <- imgMatrix[(i-1):(i+1), (j-1):(j+1)]
      # Apply convolution
      result[i, j] <- sum(neighborhood * kernel)
    }
  }
  
  return(result)
}

nonMaximumSuppression <- function(gradientMag, gradientDir) {
  
  # Import data. 
  rows <- nrow(gradientMag)
  cols <- ncol(gradientMag)
  result <- matrix(0, nrow = rows, ncol = cols)
  
  for (i in 2:(rows-1)) {
    for (j in 2:(cols-1)) {
      
      angle <- gradientDir[i, j]
      
      # Determine neighboring pixels based on gradient direction
      if (abs(angle) < pi/8 || abs(angle) > 7*pi/8) {
        
        # Horizontal edge
        neighbors <- c(gradientMag[i, j-1], gradientMag[i, j+1])
        
      } else if (abs(angle - pi/2) < pi/8 || abs(angle + pi/2) < pi/8) {
        
        # Vertical edge
        neighbors <- c(gradientMag[i-1, j], gradientMag[i+1, j])
        
      } else if (angle > 0) {
        
        # Diagonal edge (+ slope)
        neighbors <- c(gradientMag[i-1, j-1], gradientMag[i+1, j+1])
        
      } else {
        
        # Diagonal edge (- slope)
        neighbors <- c(gradientMag[i-1, j+1], gradientMag[i+1, j-1])
        
      }
      
      # Retain pixel that has the maximum gradient
      if (gradientMag[i, j] > max(neighbors, na.rm = TRUE)) {
        result[i, j] <- gradientMag[i, j]
      }
    }
  }
  
  return(result)
}

#=============================================
#
# Obtain Results
#
#=============================================

filePaths <- inputFiles(args)

results <- data.frame(
  file = character(),
  gradientStrength = numeric(),
  edgePixelCount = numeric(),
  maxGradient = numeric()
)

for (i in seq_along(filePaths)) {
  
  gradientResults <- sobelGradient(filePaths[i])
  
  results <- rbind(results, data.frame(
    file = basename(filePaths[i]),
    gradientStrength = gradientResults$gradientStrength,
    edgePixelCount = gradientResults$edgePixelCount,
    maxGradient = gradientResults$maxGradient
  ))
}

#=============================================
#
# Extract Flow Rate
#
#=============================================

results$flowRate <- as.numeric(str_extract(results$file, "\\d+\\.?\\d*"))

#=============================================
#
# Save Plots
#
#=============================================

# A quick overview of the plots. In our presentation, the data will be plotted with prism. 
p1 <- ggplot(results, aes(x = flowRate, y = gradientStrength)) +
  geom_point(size = 3, alpha = 0.7, color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  labs(
    title = "Gradient Strength vs Flow Rate",
    x = "Flow Rate",
    y = "Average Gradient Magnitude"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  )

p2 <- ggplot(results, aes(x = flowRate, y = edgePixelCount)) +
  geom_point(size = 3, alpha = 0.7, color = "darkgreen") +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  labs(
    title = "Edge Pixel Count vs Flow Rate",
    x = "Flow Rate",
    y = "Number of Edge Pixels"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  )

p3 <- ggplot(results, aes(x = flowRate, y = maxGradient)) +
  geom_point(size = 3, alpha = 0.7, color = "purple") +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  labs(
    title = "Maximum Gradient vs Flow Rate",
    x = "Flow Rate",
    y = "Maximum Gradient Value"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  )

#=============================================
#
# Print Result
#
#=============================================

# Print plots
print(p1)
print(p2)
print(p3)

print(results)

#=============================================
#
# Save Result
#
#=============================================

# Write results data frame as a file
write.csv(results,"Data/results.csv", row.names = TRUE)

