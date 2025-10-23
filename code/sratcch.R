source("code/mapping.R")
source("code/model.R")

# Generate the landscape
landscape <- create_hex_landscape()
centroids <- st_coordinates(st_centroid(landscape))
dist_mat <- as.matrix(dist(centroids))

# Dispersal
debugonce(sparse_kernel)
K <- sparse_kernel(landscape, d_bar_x = 1.5)


# Pick 7 central cells and mark as good habitat
center_ids <- c(33, 38, 39, 43, 48, 49, 53)  # adjust as needed
landscape$habitat <- ifelse(landscape$cell_id %in% center_ids, "good", "bad")
ggplot() + # visualise
  geom_sf(data = landscape, aes(fill = habitat), colour = "black", linewidth = 0) +
  scale_fill_manual(values = c("good" = "black", "bad" = "white")) +
  theme_void()


# 1. Dispersal kernel matrix
d_bar_x <- 1.5
K <- dispersal_kernel(dist_mat, d_bar_x)
K <- K / rowSums(K)  # normalize rows

# 3. Initial host distribution
H <- rep(0, nrow(landscape))
H[43] <- 1
H[1] <- 1

# 4. Apply dispersal
H_next <- K %*% H

# 5. Visualize
landscape$H <- as.numeric(H_next)
ggplot(landscape) + geom_sf(aes(fill = H)) + scale_fill_viridis_c()
