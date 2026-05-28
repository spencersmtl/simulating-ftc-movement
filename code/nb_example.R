source("code/mapping.R")
source("code/model.R")

# load and visualize landscape
nb_landscape <- load_landscape(
  "images/nb-scape-bw-thick.png", 
  scale = 0.5, # Shrink or grow image
  cellsize = 10) # Lower gives finer detail
nb_landscape <- basic_habitat_quality(nb_landscape, threshold = 0.5)
nb_landscape$Q <- ifelse(nb_landscape$type == "high", 1, 0)

# Initialize density
timesteps <- 4
H <- matrix(0, nrow = nrow(nb_landscape), ncol = timesteps)
H[nb_landscape$type == "high",1] = 1
H[,1] <- H[,1] / sum(H[,1])

# Compute neighbours and distances
max_dispersal <- 12
neighbs <- compute_neighbors(nb_landscape, max_dispersal)
dists <- compute_sparse_distance(nb_landscape, neighbs, max_dispersal)

# Disperse
scale <- 2
stay_prob <- 0.5
beta <- 1
K <- initialize_dispersal(
  dists, # sparse distance matrix
  scale, # drop-off scale
  stay_prob, # probability of staying in cell
  kernel_function = "negative_exp",
  beta, # habitat preference strength
  Q = nb_landscape$Q # habitat quality vector
)
for (i in 2:timesteps) {
  H[,i] <- disperse(
    kernel = K,
    density = H[,i-1],
    survival = 1
  )
}

# just landscape
visualise_landscape(nb_landscape)
# landscape with initial density
dots <- dot_density_points(nb_landscape, density = H[,1], dot_clutter =5)
visualise_landscape(nb_landscape, dots = dots, dotsize = 1, show_legend = FALSE)
# after 1 round of dispersal
dots <- dot_density_points(nb_landscape, density = H[,2], dot_clutter = 5)
visualise_landscape(nb_landscape, dots = dots, dotsize = 1, show_legend = FALSE)
# after 2 rounds
dots <- dot_density_points(nb_landscape, density = H[,3], dot_clutter = 5)
visualise_landscape(nb_landscape, dots = dots, dotsize = 1, show_legend = FALSE)
# after 3 rounds
dots <- dot_density_points(nb_landscape, density = H[,4], dot_clutter = 5)
visualise_landscape(nb_landscape, dots = dots, dotsize = 1, show_legend = FALSE)
