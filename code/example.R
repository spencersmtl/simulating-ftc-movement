source("code/mapping.R")
source("code/model.R")

# Generate and process the landscape
# landscape <- create_hex_landscape()

# Load an image as landscape
landscape <- load_landscape("images/blobs.png",
                            scale = 0.5,# Set scale <1 to shrink image
                            cellsize = 80) # set cellsize <8 to make a finer hex grid overlay
visualise_landscape(landscape)

# Generate sparse dispersal kernel
K <- sparse_kernel(landscape, d_bar_x = 1)

# Initialize host distribution
H <- rep(0, nrow(landscape))
H[c(6)] <- 1

# Disperse
H_next <- K %*% H

# Visualize
landscape$H <- as.numeric(H_next)
ggplot(landscape) + 
  geom_sf(aes(fill = H)) + 
  scale_fill_viridis_c() + 
  theme_void()
