library(sf)
library(nngeo)
library(Matrix)

# Host growth
host_growth <- function(H, P, lambda, a) {
  lambda * H * exp(-log(lambda) * H) * exp(-P * a)
}

# Parasitoid growth
parasitoid_growth <- function(H, P, lambda, a) {
  a * H * (1 - exp(-P)) * exp(-log(lambda) * H)
}

# Dispersal kernel
dispersal_kernel <- function(d_ij, d_bar_x) {
  (2 / (pi * d_bar_x^2)) * exp(-2 * d_ij / d_bar_x)
}

# Sparse distance matrix computation
compute_sparse_distance <- function(landscape, 
                                    avg_host_dispersal, 
                                    max_host_dispersal_mult = 1.1, 
                                    normalize_rows = TRUE) 
{
  
  n <- nrow(landscape) # number of cells
  centroids <- st_centroid(landscape)
  centroid_dist <- as.numeric(st_distance(centroids[1, ], centroids[2, ])) # distance between adjacent centroids (cellsize)
  avg_host_dispersal <- avg_host_dispersal * centroid_dist
  max_host_dispersal <- max_host_dispersal_mult * avg_host_dispersal # max dispersal distance
  
  # Get neighbor coordinates and distances
  neighbors <- st_nn(
    x = centroids, 
    y = centroids, 
    k = n,   # large number to get all within maxdist
    maxdist = max_host_dispersal, 
    returnDist = TRUE,   # we want distances
    progress = FALSE
  )
  origins <- rep(1:n, lengths(neighbors$nn)) # origin indices
  destinations <- unlist(neighbors$nn) # destination indices
  neighbor_distances <- as.numeric(unlist(neighbors$dist)) # distances
  browser()
  # assemble as sparse matrix
  D <- sparseMatrix(i = origins, 
                    j = destinations, 
                    x = neighbor_distances, 
                    dims = c(n, n))
  
  return(D)
}