source("code/mapping.R")
source("code/model.R")

# Load an image as landscape
landscape <- load_landscape(
  "images/blobs.png",
  scale = 0.5,# Set scale <1 to shrink image
  cellsize = 100) |> # set cellsize <8 to make a finer hex grid overlay
  basic_habitat_quality(threshold = 0.5)

# Generate sparse dispersal kernel
D <- compute_sparse_distance(
  landscape, 
  avg_host_dispersal = 1, 
  max_host_dispersal_mult = 1.1, 
  normalize_rows = TRUE)

compute_dispersal