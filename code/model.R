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
  
  # Prepare coordinates ands cale
  n <- nrow(landscape) # number of cells
  centroids <- st_centroid(landscape) # get centroids
  coords <- st_coordinates(centroids) # get centroid coordinates
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
    stay_prob,
    kernel_function = c("negative_exp","gaussian","cauchy")) 
{
  kernel_function <- match.arg(kernel_function)
  kernel_fun <- switch(
    kernel_function,
    negative_exp = function(d, s) exp(-d / s),
    gaussian     = function(d, s) exp(-(d^2) / (2 * s^2)),
    cauchy       = function(d, s) 1 / (1 + (d / s)^2)
  )
  
  # scale relative to minimum non-zero distance
  scale <- scale * min(distances@x[distances@x != 0]) 

  # apply kernel to stored sparse entries
  K <- distances
  K@x <- kernel_fun(K@x, scale)
  
  # normalize rows to sum to 1 and set diagonal to proportion that do not disperse
  diag(K) <- 0
  row_totals <- Matrix::rowSums(K)
  row_totals[row_totals == 0] <- 1 # avoid division by zero
  K <- Diagonal(x = (1-stay_prob) / row_totals) %*% K
  K <- K + Diagonal(nrow(K), stay_prob)
  
  K
}

# single timestep dispersal (call every timestep)
disperse <- function(
    kernel, 
    density, 
    survival = 1) # Survival probability during dispersal
{
  as.numeric(t(kernel) %*% density) * survival
}
