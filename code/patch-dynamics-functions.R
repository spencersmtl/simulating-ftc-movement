library(tidyverse)
library(viridis)
library(png)
library(paletteer)

random_walk_step <- function(current_x, current_y) { # function to perform random walks for both FTC and fly
  dx <- sample(c(-8, -4, 0, 4, 8), 1) # x direction
  dy <- sample(c(-8, -4, 0, 4, 8), 1) # y direction
  return(c(current_x + dx, current_y + dy))
}

H_growth <- function(H, r, K, t, n) { # patch wise host growth term
  return(r * H[t - 1, n] * (1 - H[t - 1, n] / K))
}

H_shrink <- function(H, P, a, handle_t, t, n) { # patch wise host shrink term
  return((a * H[t - 1, n] * P[t - 1, n]) / (1 + a * handle_t * H[t - 1, n]))
}

H_delta <- function(H, P, r, a, K, handle_t, t, n) { # patch wise host change
  growth <- H_growth(H, r, K, t, n)
  shrink <- H_shrink(H, P, a, handle_t, t, n)
  return(growth - shrink)
}

P_growth <- function(H, P, e, a, handle_t, t, n) { # patch wise parasitoid growth term
  return((e * a * H[t - 1, n] * P[t - 1, n]) / (1 + a * handle_t * H[t - 1, n]))
}

P_shrink <- function(P, d, t, n) { # patch wise parasitoid shrink term
  return(d * P[t - 1, n])
}

P_delta <- function(H, P, e, a, d, handle_t, t, n) { # patch wise parasitoid change
  growth <- P_growth(H, P, e, a, handle_t, t, n)
  shrink <- P_shrink(P, d, t, n)
  return(growth - shrink)
}


for(i in 1:n_patches) # Compute CR for each patch
{
  if(patch_death_tracker[i]>20) {
    landscape_data$clusters[landscape_data$clusters == i] <- 0
    patch_occupied[i,] <- FALSE
  }
  if(patch_occupied[i,1] & patch_occupied[i,2]){ # If occupied by both FTC and fly
    
    cr_patches_r[t,i] = cr_patches_r[t-1,i] + a * cr_patches_r[t-1,i] *  # Resource AKA FTC dynamics
      (1 - cr_patches_r[t-1,i] / K) - b * cr_patches_r[t-1,i] / (A + cr_patches_r[t-1,i]) * cr_patches_c[t-1,i]
    cr_patches_c[t,i] = cr_patches_c[t-1] + e * # Consumer AKA Fly dynamics
      ((cr_patches_r[t-1] / (A + cr_patches_r[t-1])) - (B / (B + A))) * cr_patches_c[t-1] 
  }
  if(patch_occupied[i,1] & patch_occupied[i,2] == FALSE){ # If occupied by JUST FTC
    
    cr_patches_r[t,i] = cr_patches_r[t-1] + a * cr_patches_r[t-1] *  # Resource AKA FTC dynamics
      (1 - cr_patches_r[t-1] / K)
  }
  if(cr_patches_r[t,i]>80) { # Track if patch is going extinct
    patch_death_tracker[i] <- patch_death_tracker[i]+1} else {patch_death_tracker[i]=0}
}