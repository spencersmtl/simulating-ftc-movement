source("code/mapping.R")
source("code/model.R")
{
  timesteps <- 2
  
  # load landscape
  landscape <- load_landscape(
    "images/blobs.png",
    scale = 0.5, # Set scale <1 to shrink image
    cellsize = 80 # set cellsize <8 to make a finer hex grid overlay
  )
  landscape <- basic_habitat_quality( # set habitat quality threshold
    landscape, 
    threshold = 0.5
  ) 
  
  # Generate sparse distance matrix
  D <- compute_sparse_distance(landscape, max_host_dispersal = 1)
  
  # Generate dispersal kernel
  K <- initialize_dispersal(
    D, 
    stay_prob = 0.5,
    scale = 1, 
    kernel_function = "negative_exp"
  )
  
  # Initialize density
  H <- matrix(NA_real_, nrow = nrow(landscape), ncol = timesteps)
  H[,1] <- c(0,0,0,1,0,0,0,0,0) # initial host population in cell 4
  
  # Disperse once
  H[,2] <- disperse( 
    kernel = K, 
    density = H[,1], 
    survival = 1 
  )
  
  visualise_landscape(landscape, density = H[,2])
}
#   row_sums <- rowSums(K)