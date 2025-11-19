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

# Sparse dispersal
sparse_kernel <- function(landscape, d_bar_x, maxdist_mult = 4, normalize_rows = TRUE) {
  n <- nrow(landscape) # number of cells
  maxdist <- maxdist_mult * d_bar_x # max dispersal distance
  
  # Get normalized centroids
  coords <- as.matrix(
    sf::st_set_geometry(landscape, NULL)[, c("centroid_x", "centroid_y")])
  
  # Assemble dispersal matrix
  dists <- as.matrix(dist(coords)) # Compute pairwise distance matrix
  dists[dists > maxdist] <- 0 # ignore distances beyond maxdist
  diag(dists) <- 0   # distance to self is zero
  weights <- dispersal_kernel(dists, d_bar_x)   # Compute kernel weights
  weights <- (weights + t(weights)) / 2  # Ensure symmetry
  weights <- weights / rowSums(weights)  # Normalize rows
  
  K <- Matrix(weights, sparse = TRUE)   # Convert to sparse matrix
  
  return(K)
}
  