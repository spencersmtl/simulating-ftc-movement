source("code/mapping.R")
source("code/model.R")
{
  timesteps <- 9
  max_host_dispersal <- 8 # number of adjacent cells
  scale <- 1 # drop-off strength
  stay_prob <- 0.5  # probability of staying in cell
  beta <- 1 # habitat preference
  survival <- 1 # dispersal survival rate
  
  # load landscape
  landscape <- load_landscape("images/nb-scape-bw-thick.png", scale = 0.5, # Set scale <1 to shrink image
    cellsize = 15) # set cellsize <8 to make a finer hex grid overlay
  landscape <- basic_habitat_quality(landscape, threshold = 0.5) # set habitat quality threshold
  landscape$Q <- ifelse(landscape$type == "high", 1, 0)
  
  # Compute neighbours
  neighbs <- compute_neighbors(landscape, max_host_dispersal)
  
  # Compute neighbour distances
  dists <- compute_sparse_distance(landscape, neighbs, max_host_dispersal)
  
  # Generate dispersal kernel
  K <- initialize_dispersal(
    dists, # sparse distance matrix
    scale, # drop-off scale
    stay_prob, # probability of staying in cell
    kernel_function = "negative_exp",
    beta, # habitat preference strength
    Q = landscape$Q # habitat quality vector
  )
  
  # Initialize density
  H <- matrix(0, nrow = nrow(landscape), ncol = timesteps)
  H[landscape$type == "high",1] = 1
  H[,1] <- H[,1] / sum(H[,1])
  
  # Disperse
  for (i in 2:timesteps) {
    H[,i] <- disperse(
      kernel = K,
      density = H[,i-1],
      survival
    )
  }
}

visualise_landscape(landscape)

visualise_landscape(landscape, density = H[,timesteps])

{
  dots <- dot_density_points(
    landscape, density = H[,4], dot_clutter = 5)
  visualise_lansdscape(landscape, dots = dots, dotsize = 1, show_legend = FALSE)
}

# Time-evolved with dot-density
{
  plots <- lapply(1:timesteps, function(t) {
    density <- H[, t]
    dots <- dot_density_points(landscape, density, dot_clutter = 3)
    visualise_landscape(landscape, dots = dots, show_legend = FALSE)
  })
  wrap_plots(plots, ncol = 3)  # adjust ncol for layout
}
