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
compute_sparse_distance <- function(landscape, 
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
  sparse_distance_matrix <- sparseMatrix(i = origins, 
                                         j = destinations, 
                                         x = neighbor_distances, 
                                         dims = c(n, n))
  
  return(sparse_distance_matrix)
}

compute_dispersal <- function(landscape,
                              compute_sparse_distance,
                              average_host_dispersal = 1,
)
{
  # Example kernel function: exponential (Laplace-like)
  kernel_fun <- function(d, scale) {
    exp(- d / scale)
  }
  # You can swap in other families: gaussian: exp(-(d^2)/(2*scale^2)),
  # cauchy/fat-tailed: 1 / (1 + (d/scale)^alpha), etc.
  
  initialize_dispersal <- function(D, density, scale,
                                   kernel_fun = kernel_fun,
                                   normalize_rows = TRUE,
                                   dispersal_survival = 1) {
    stopifnot(inherits(D, "dgCMatrix") || inherits(D, "dgTMatrix"),
              length(density) == nrow(D))
    
    # 1) convert distances in the non-zero entries to kernel values
    K <- D  # copy the sparse structure
    K@x <- kernel_fun(K@x, scale)  # only transforms non-zero entries
    
    # 2) (optional) normalize rows so that each row sums to 1 (row-stochastic)
    if (normalize_rows) {
      rs <- rowSums(K)             # efficient for sparse matrices
      rs[rs == 0] <- 1             # avoid division by zero (isolated cells)
      inv_rs <- 1 / rs
      K <- Diagonal(x = inv_rs) %*% K  # left-multiply scales rows
    }
    
    # 3) Apply dispersal (rows=origin convention)
    # If density is a numeric vector (length n), compute post-dispersal density.
    # For rows-as-origin, new_density = t(K) %*% density
    new_density <- as.numeric(t(K) %*% density) * dispersal_survival
    
    # 4) Quick mass check (should be ~ same sum if kernel is stochastic and survival=1)
    mass_before <- sum(density)
    mass_after <- sum(new_density)
    
    list(K = K,
         post_density = new_density,
         mass_before = mass_before,
         mass_after = mass_after)
  }
}
