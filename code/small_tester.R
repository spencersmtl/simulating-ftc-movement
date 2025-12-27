source("code/mapping.R")
source("code/model.R")
{
  timesteps <- 20
  max_host_dispersal <- 2 # number of adjacent cells
  beta <- 0.8 # habitat preference
  
  # load landscape
  landscape <- load_landscape("images/blobs.png", scale = 0.5, # Set scale <1 to shrink image
    cellsize = 80) # set cellsize <8 to make a finer hex grid overlay
  landscape <- basic_habitat_quality(landscape, threshold = 0.5) # set habitat quality threshold
  landscape$Q <- ifelse(landscape$type == "high", 1, 0)
  
  # Compute neighbours
  neighbs <- compute_neighbors(landscape, max_host_dispersal)
  
  # Compute neighbour distances
  dists <- compute_sparse_distance(landscape, neighbs, max_host_dispersal)
  
  # Generate dispersal kernel
  K <- initialize_dispersal(
    dists, # sparse distance matrix
    stay_prob = 0.5, # probability of staying in cell
    scale = 1, # drop-off scale
    kernel_function = "negative_exp",
    beta, # habitat preference strength
    Q = landscape$Q # habitat quality vector
  )
  
  # Initialize density
  H <- matrix(NA_real_, nrow = nrow(landscape), ncol = timesteps)
  H[,1] <- 1
  H[,1] <- H[,1] / sum(H[,1])
  
  # Disperse
  for (i in 2:timesteps) {
    H[,i] <- disperse(
      kernel = K,
      density = H[,i-1],
      survival = 1
    )
  }
  
  visualise_landscape(landscape, density = H[,20])
}
