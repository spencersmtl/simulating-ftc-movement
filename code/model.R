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

compute_neighbors <- function(
    landscape,
    max_host_dispersal = nrow(landscape))
{
  # Prepare coordinates and scale
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
  
  neighbors_list <- lapply(1:n, function(i) {
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
    neighbors <- candidate_cells[!is.na(candidate_cells)]
    
    neighbors
  })
  neighbors_list
}

compute_sparse_distance <- function(
    landscape, 
    neighbors_list,
    max_host_dispersal = nrow(landscape)) 
{
  # Prepare coordinates and scale
  n <- nrow(landscape) # number of cells
  centroids <- st_centroid(landscape) # get centroids
  coords <- st_coordinates(centroids) # get centroid coordinates
  centroid_dist <- min(as.numeric(dist(coords))) # distance between adjacent centroids (cellsize)
  grid_width <- max(coords[,1]) - min(coords[,1]) + centroid_dist
  max_host_dispersal <- max_host_dispersal * centroid_dist
  # Row width
  y_coord_row1 <- coords[1, 2]
  x_coords_row1 <- coords[coords[,2] == y_coord_row1, 1]
  centroid_row_width <- max(x_coords_row1) - min(x_coords_row1)
  
  # Precompute total number of non-zero entries and set-up matrix
  nonzeros <- sum(lengths(neighbors_list))
  rows <- integer(nonzeros)
  cols <- integer(nonzeros)
  vals <- numeric(nonzeros)
  
  k <- 1 # hexdex (hex index) for distance matrix entries
  for (i in seq_len(n)) {
    neighbors <- neighbors_list[[i]]
    m <- length(neighbors)
    if (m == 0) next
    dx <- coords[neighbors, 1] - coords[i, 1]
    dy <- coords[neighbors, 2] - coords[i, 2]
    d <- sqrt(dx^2 + dy^2)
    
    # Split local vs wrapped
    local_mask <- d <= max_host_dispersal + 1
    wrapped_mask <- !local_mask
    
    if (length(wrapped_mask) > 0) {
      wrapped_neighbors <- neighbors[wrapped_mask]  # track neighbor indices
      
      # Copy coordinates of wrapped neighbors
      wrapped_coords <- coords[wrapped_neighbors, , drop = FALSE]
      
      # Shift wrapped hexes horizontally by grid width
      shifted_coords <- wrapped_coords
      flip <- sign(wrapped_coords[,1] - coords[i,1])
      shifted_coords[,1] <- wrapped_coords[,1] - 
        flip * (centroid_row_width + centroid_dist)
      
      # Compute Euclidean distances to ghost hexes
      dx_wrapshifted <- coords[i,1] - shifted_coords[,1]
      dy_wrapshifted <- shifted_coords[,2] - coords[i,2]
      d[wrapped_mask] <- sqrt(dx_wrapshifted^2 + dy_wrapshifted^2)
    }
    
    rows[k:(k+m-1)] <- i
    cols[k:(k+m-1)] <- neighbors
    vals[k:(k+m-1)] <- d
    
    k <- k + m
  }
  
  sparseMatrix(i = rows, 
               j = cols, 
               x = vals, 
               dims = c(n, n))
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
