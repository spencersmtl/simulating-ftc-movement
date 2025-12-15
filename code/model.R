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

  landscape <- landscape %>%
    mutate(
      x = as.integer(factor(round(coords[,1], 6))),
      y = as.integer(factor(round(coords[,2], 6)))
    )

  n_x <- max(landscape$x)
  n_y <- max(landscape$y)

  k <- ceiling(max_host_dispersal / centroid_dist)

  neighbors_list <- lapply(
    1:n,
    function(i) {
      origin_col <- landscape$x[i]
      origin_row <- landscape$y[i]

      candidate_cols <- ((origin_col + (-(k*2):(k*2)) - 1) %% n_x) + 1 # wrap in x direction
      candidate_rows <- origin_row + (-k:k) # don't wrap in y direction
      candidate_rows <- candidate_rows[
        candidate_rows >= 1 & candidate_rows <= n_y]

      candidate_cells <- as.vector(
        outer(candidate_rows, candidate_cols, Vectorize(
          function(y, x) {
            idx <- landscape$cell_id[landscape$y==y & landscape$x==x]
            if(length(idx)==0) NA else idx
          }
        ))
      )
      candidate_cells <- candidate_cells[!is.na(candidate_cells)]

      candidate_cells
    })

  dx <- coords[candidate_cells, 1] - coords[i, 1]
  dx <- dx - (max(coords[,1]) - min(coords[,1])) * round(dx / (max(coords[,1]) - min(coords[,1])))  # wrap in x
  dy <- coords[candidate_cells, 2] - coords[i, 2]  # no wrap in y

  dists <- sqrt(dx^2 + dy^2)

  # keep only neighbors within max_host_dispersal
  keep <- candidate_cells[dists <= max_host_dispersal]
  dists <- dists[dists <= max_host_dispersal]


  
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
