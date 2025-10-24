source("code/mapping.R")
source("code/model.R")

# Generate and process the landscape
landscape <- create_hex_landscape()
centroids <- st_coordinates(st_centroid(landscape))
dist_mat <- as.matrix(dist(centroids))

# Generate sparse dispersal kernel
K <- sparse_kernel(landscape, d_bar_x = 1.5)

# Initialize host distribution
H <- rep(0, nrow(landscape))
H[c(1, 33, 38, 39, 43)] <- 1

# Disperse
H_next <- K %*% H

# Visualize
landscape$H <- as.numeric(H_next)
ggplot(landscape) + geom_sf(aes(fill = H)) + scale_fill_viridis_c()
