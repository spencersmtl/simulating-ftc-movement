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

# Sparse distance matrix computation
compute_sparse_distance <- function(
    landscape, 
    max_host_dispersal = nrow(landscape), 
    normalize_rows = TRUE) 
{
  n <- nrow(landscape) # number of cells
  centroids <- st_centroid(landscape)
  centroid_dist <- as.numeric(st_distance(centroids[1, ], centroids[2, ])) # distance between adjacent centroids (cellsize)
  max_host_dispersal <- max_host_dispersal * centroid_dist # max dispersal distance

  # Get neighbor coordinates and distances
  neighbors <- st_nn(
    x = centroids, 
    y = centroids, 
    k = n,   # large number to get all within maxdist
    maxdist = max_host_dispersal+1e-9, # include maxdist 
    returnDist = TRUE,   # we want distances
    progress = FALSE
  )
  origins <- rep(1:n, lengths(neighbors$nn)) # origin indices
  destinations <- unlist(neighbors$nn) # destination indices
  neighbor_distances <- as.numeric(unlist(neighbors$dist)) # distances
  
  # assemble as sparse matrix
  sparse_distances <- sparseMatrix(i = origins, 
                                   j = destinations, 
                                   x = neighbor_distances, 
                                   dims = c(n, n))
  
  return(sparse_distances)
}

# build a normalized dispersal kernel (works with sparse matrices)
initialize_dispersal <- function(
    distances,
    scale,
    kernel_function = c("neg_exponential","gaussian","cauchy")) 
{
  kernel_function <- match.arg(kernel_function)
  kernel_fun <- switch(
    kernel_function,
    neg_exponential = function(d, s) exp(-d / s),
    gaussian        = function(d, s) exp(-(d^2) / (2 * s^2)),
    cauchy          = function(d, s) 1 / (1 + (d / s)^2)
  )
  
  # apply kernel to stored sparse entries
  K <- distances
  K@x <- kernel_fun(K@x, scale)
  
  # normalize rows
  row_totals <- Matrix::rowSums(K)
  row_totals[row_totals == 0] <- 1
  K <- Diagonal(x = 1 / row_totals) %*% K
  
  K
}

# single timestep dispersal (call every timestep)
disperse <- function(
    kernel, 
    density, 
    survival = 1) # Survival probability during dispersal
{
  as.numeric(kernel %*% density) * survival
}

# convenience: return a disperser closure
make_disperser <- function(distances, scale, kernel_function, survival = 1) {
  K <- initialize_dispersal(distances, scale, kernel_function)
  function(density) disperse(K, density, survival)
}
