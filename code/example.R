source("code/mapping.R")
source("code/model.R")

# Generate and process the landscape
# landscape <- create_hex_landscape()

# Load an image
landscape <- load_landscape("images/blobs.png",scale = 0.5,cellsize = 8)
visualise_landscape(landscape)
centroids <- st_coordinates(st_centroid(landscape))
dist_mat <- as.matrix(dist(centroids))

# Generate sparse dispersal kernel
K <- sparse_kernel(landscape, d_bar_x = 1.5)

# Initialize host distribution
H <- rep(0, nrow(landscape))
H[c(1, 43)] <- 1

# Disperse
H_next <- K %*% H

# Visualize
landscape$H <- as.numeric(H_next)
ggplot(landscape) + geom_sf(aes(fill = H)) + scale_fill_viridis_c() + theme_void()
