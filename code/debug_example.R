source("code/mapping.R")
source("code/model.R")

# Load an image as landscape
landscape <- load_landscape(
  "images/blobs.png",
  scale = 0.5,# Set scale <1 to shrink image
  cellsize = 10) |> # set cellsize <8 to make a finer hex grid overlay
  basic_habitat_quality(threshold = 0.5)

visualise_landscape(landscape)

# Generate sparse distance matrix
D <- compute_sparse_distance(landscape, 
                             avg_host_dispersal = 1, 
                             max_host_dispersal_mult = 1.1, 
                             normalize_rows = TRUE)

# Initialize host distribution
H <- rep(0, nrow(landscape))
H[c(5)] <- 1




# 
# # Disperse
# H_next <- D %*% H
# 
# # Visualize
# landscape$H <- as.numeric(H_next)
# ggplot(landscape) + 
#   geom_sf(aes(fill = H)) + 
#   scale_fill_viridis_c() + 
#   theme_void()
