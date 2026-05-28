source("code/mapping.R")
source("code/model.R")
# load landscape
{
landscape <- load_landscape("images/blobs.png", scale = 0.5, # Set scale <1 to shrink image
                            cellsize = 70) # set cellsize <8 to make a finer hex grid overlay
landscape <- basic_habitat_quality(landscape, threshold = 0.5) # set habitat quality threshold
landscape$quality <- ifelse(landscape$type == "high", 1, 0)
quality <- landscape$quality
}
timesteps <- 10

{ # Hosts
max_dispersal <- 2 # number of adjacent cells, set to >=1
scale <- 1 # drop-off strength
stay_prob <- 0.9  # probability of staying in cell
beta <- 2 # habitat preference
survival <- 1 # dispersal survival rate
kernel_function <- "negative_exp" # use negative_exp, cauchy, or gaussian

neighbs <- compute_neighbors(landscape, max_dispersal) # Compute neighbours
dists <- compute_sparse_distance(landscape, neighbs, max_dispersal) # Compute neighbour distances
K <- initialize_dispersal(  # Generate dispersal kernel
  dists, scale, stay_prob, kernel_function, beta, quality) 
H <- matrix(0, nrow = nrow(landscape), ncol = timesteps)
H[landscape$type == "high",1] = 1
H[,1] <- H[,1] / sum(H[,1])
for (i in 2:timesteps) { # Disperse
  H[,i] <- disperse(kernel = K, density = H[,i-1], survival)}
}

{ # Parasitoids
  max_dispersal <- 2 # number of adjacent cells, set to >=1
  scale <- 1 # drop-off strength
  stay_prob <- 0.9  # probability of staying in cell
  beta <- 2 # habitat preference
  survival <- 1 # dispersal survival rate
  kernel_function <- "negative_exp" # use negative_exp, cauchy, or gaussian
  
  neighbs <- compute_neighbors(landscape, max_dispersal) # Compute neighbours
  dists <- compute_sparse_distance(landscape, neighbs, max_dispersal) # Compute neighbour distances
  K <- initialize_dispersal(  # Generate dispersal kernel
    dists, scale, stay_prob, kernel_function, beta, quality) 
  P <- matrix(0, nrow = nrow(landscape), ncol = timesteps)
  P[landscape$type == "high",1] = 1
  P[,1] <- P[,1] / sum(P[,1])
  for (i in 2:timesteps) { # Disperse
    P[,i] <- disperse(kernel = K, density = P[,i-1], survival)}
}

# Time-evolved with dot-density
{
  plots <- lapply(1:timesteps, function(t) {
    density <- H[, t]
    dots <- dot_density_points(landscape, density, dot_clutter = 10)
    visualise_landscape(landscape, dots = dots, show_legend = FALSE)
  })
  wrap_plots(plots, ncol = 4)  # adjust ncol for layout
}

# Time-evolved with numerical density
{
  plots <- lapply(1:timesteps, function(t) {
    visualise_landscape(landscape, density = H[, t], show_legend = FALSE)
  })
  wrap_plots(plots, ncol = 4)  # adjust ncol for layout
}
