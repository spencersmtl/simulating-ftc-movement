library(sf)
library(spdep)
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

# Sparsify dispersal
sparse_kernel <- function(landscape, d_bar_x, maxdist_mult = 4) {
  coords <- st_coordinates(st_centroid(landscape)) # Set coordinates to cell centroids
  n <- nrow(coords) # number of cells
  maxdist <- maxdist_mult * d_bar_x # max dispersal distance
  
  # Build sparse dispersal kernel. Ignore computing dispersal to cells further than maxdist
  nb <- dnearneigh(x = coords, d1 = 0, d2 = maxdist, longlat = FALSE)   # List of possible destinations for each cell
  
  rows <- rep(1:n, lengths(nb) + 1) # big vector of all origin cells. Repeat origin for each destination
  cols <- unlist( # big vector of all destinations for each origin cell. unlist turns list into one big vector
    lapply(
      1:n, # for each cell
      function(i) c(i, nb[[i]]))) # make pairs: `cell index` paired with `possible destination`
  
  diffcells <- coords[rows, , drop = FALSE] - coords[cols, , drop = FALSE] # differences between coordinates of possible paths
  dists <- sqrt(rowSums(diffcells^2))
  
  weights <- dispersal_kernel(dists, d_bar_x)   # compute kernel weights
  
  # assemble sparse matrix (rows = origin, cols = destination)
  K <- sparseMatrix(i = rows, j = cols, x = weights, dims = c(n, n))
  # Here normalize rows if you don't want absorbing conditions
  
  return(K)  
}

# Example usage:
# landscape <- create_hex_landscape()
# K_sparse <- build_sparse_kernel(landscape, d_bar_x = 1.5, maxdist_mult = 4)
# H <- numeric(nrow(landscape)); H[30] <- 1
# H_next <- as.numeric(K_sparse %*% H)
