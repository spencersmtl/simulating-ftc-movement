source("code/mapping.R")
source("code/model.R")

# EXTREMELY SMALL
landscape <- load_landscape(
  "images/blobs.png",
  scale = 0.5, # Set scale <1 to shrink image
  cellsize = 8# set cellsize <8 to make a finer hex grid overlay
) |> 
  basic_habitat_quality(threshold = 0.5) 
visualise_landscape(landscape)

# PRETTY SMALL
landscape <- load_landscape(
  "images/blobs.png",
  scale = 0.5, # Set scale <1 to shrink image
  cellsize = 80 # set cellsize <8 to make a finer hex grid overlay
) |> 
  basic_habitat_quality(threshold = 0.5) 
visualise_landscape(landscape)

# Generate sparse dispersal kernel
D <- compute_sparse_distance(
  landscape, 
  avg_host_dispersal = 1, 
  max_host_dispersal_mult = 1.1, 
  normalize_rows = TRUE)
D

