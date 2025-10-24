library(tidyverse)
library(viridis)
library(png)
library(paletteer)

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

walk_step <- function(n_steps, step_length, cluster_id, landscape_data) {
  # Filter valid positions within the specified patch
  walker_path <- landscape_data %>%
    filter(clusters == cluster_id) %>%
    slice_sample(n = 1) %>%
    dplyr::select(x, y)
  
  for (k in 1:n_steps) {
    # Take step
    new_x <- walker_path[k,1] + sample(step_length, 1) 
    new_y <- walker_path[k,2] + sample(step_length, 1)
    if (new_x <= 0 | new_y <= 0) break
    else walker_path <- rbind(walker_path, c(new_x, new_y))
  }
  return(walker_path)
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