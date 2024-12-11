# Issue: need to keep track of which patch gets occupied
# Load necessary libraries
library(png)
library(viridis)
library(tidyverse)

# Load the landscape data
landscape_data <- readRDS("data/50x50_scale2a_simple_patches_data.rds")

# Initialize occupancy status of patch
landscape_data$resource_occupied <- landscape_data$clusters == 1
landscape_data$consumer_occupied <- landscape_data$clusters == 1

# Initial plot for fun
patch_map <- ggplot(landscape_data) +
  geom_raster(aes(x = x, y = -y, fill = as.factor(clusters))) +
  scale_fill_viridis(discrete = TRUE, option = "D", na.value = "white") +  # "D" is just one of the available options
  theme(legend.position = "none") +
  coord_equal()
patch_map

# Function to perform random walks for both FTC and fly
random_walk_step <- function(current_x, current_y) {
  dx <- sample(c(-8, -4, 0, 4, 8), 1) # x direction
  dy <- sample(c(-8, -4, 0, 4, 8), 1) # y direction
  return(c(current_x + dx, current_y + dy))
}

# Parameters
{
  n_patches = length(unique(na.omit(landscape_data$clusters))) # Number of patches omitting NAs
  t <- 1               # Initial time
  t_final <- 200       # number of time steps
  R0 <- 60             # Initial FTC in patch
  C0 <- 20              # Initial Fly in patch
  a <- 3               # Ability of FTC to avoid flies
  A <- 50              # FTC half-saturation constant ???
  b <- 5               # Fly consumption ability
  B <- 35              # Fly half-saturation constant
  e <- 0.5             # Extinction rate of the fly
  K <- 100               # FTC Carrying capacity
  
  # Initialize objects before loop
  cr_patches <- vector("list", n_patches) # cr series all patches
  cr_patch <- data.frame(R = numeric(t_final), # cr series patch 1
                         C = numeric(t_final))
  cr_patch$R[1] = R0 # Set initials
  cr_patch$C[1] = C0
  cr_patches[[1]] <- cr_patch
  
  per_t_ftc_walker_paths <- vector("list", 1) # All paths in all patches per time step
  per_t_fly_walker_paths <- vector("list", 1)
  
  all_patchwise_ftc_paths <- vector("list", t_final/5 - 2) # All paths in all patches for all time steps
  all_patchwise_fly_paths <- vector("list", t_final/5 - 2) 
  
  all_ftc_paths <- vector("list", n_patches)
  
  patch_occupied <- data.frame(R_occupied = rep(FALSE, n_patches), 
                               C_occupied = rep(FALSE, n_patches))
  patch_occupied$R_occupied[1] <- TRUE
  patch_occupied$C_occupied[1] <- TRUE#landscape_data[landscape_data$clusters == i,][1,]$resource_occupied
  
  patch_death_tracker <- rep(0,n_patches)
}

# Subset data for inhabited patch
current_cluster <- filter(landscape_data, clusters==t)

for(t in 1:(t_final-1)) # Time series loop containing CR and random walking
{
  for(i in 1:n_patches) # Compute CR for each patch
  { 
    if(patch_occupied[i,1] & patch_occupied[i,2]){ # If occupied by both FTC and fly
      
      cr_patch$R[t+1] = cr_patch$R[t] + a * cr_patch$R[t] *  # Resource AKA FTC dynamics
        (1 - cr_patch$R[t] / K) - b * cr_patch$R[t] / (A + cr_patch$R[t]) * cr_patch$C[t]
      cr_patch$C[t+1] = cr_patch$C[t] + e * # Consumer AKA Fly dynamics
        ((cr_patch$R[t] / (A + cr_patch$R[t])) - (B / (B + A))) * cr_patch$C[t] 
    }
    if(patch_occupied[i,1] & patch_occupied[i,2] == FALSE){ # If occupied by JUST FTC
      
      cr_patch$R[t+1] = cr_patch$R[t] + a * cr_patch$R[t] *  # Resource AKA FTC dynamics
        (1 - cr_patch$R[t] / K)
    }
    cr_patches[[i]] <- cr_patch
    if(cr_patch$R[t+1]>80) { # Track if patch is going extinct
      patch_death_tracker[i] <- patch_death_tracker+1} else patch_death_tracker=0
  }
  
  if(t %% 5 == 0) # every 5 time steps...
  {
    for(i in 1:n_patches) # Compute random walks for each patch
    {
      if(patch_occupied[i,1]) # FTC
      {
        patchwise_ftc_walkers <- floor(cr_patches[[i]]$R[t+1]/10) # Number of walkers is proportional to patch population size
        ftc_walker_paths <- vector("list", patchwise_ftc_walkers)
        
        start_positions_ftc <- landscape_data %>% # FTC starting positions
          filter(resource_occupied, clusters == i) %>%
          sample_n(patchwise_ftc_walkers, replace = TRUE) %>%
          select(x, y)
        
        for (j in 1:patchwise_ftc_walkers) # Simulate each ftc walker
        {
          # starting position for ftc_walker
          current_x <- start_positions_ftc$x[j]
          current_y <- start_positions_ftc$y[j]
          ftc_walker_path <- data.frame(x = current_x, y = current_y)
          
          for (k in 1:10) { # Take the steps
            new_coords <- random_walk_step(current_x, current_y)
            current_x <- new_coords[1]
            current_y <- new_coords[2]
            
            # Constrain the ftc_walker's position within the image bounds
            current_x <- max(min(current_x, max(landscape_data$x)), 1)
            current_y <- max(min(current_y, max(landscape_data$y)), 1)
            
            ftc_walker_path <- rbind(ftc_walker_path, c(current_x, current_y))
          }
          ftc_walker_paths[[j]] <- ftc_walker_path # Patch wise paths
        }
      }
      if(patch_occupied[i,2]) # Fly
      {
        patchwise_fly_walkers <- floor(cr_patches[[i]]$C[t+1]/10)
        fly_walker_paths <- vector("list", patchwise_fly_walkers)
        
        start_positions_fly <- landscape_data %>% # Fly starting positions
          filter(consumer_occupied, clusters == i) %>%
          sample_n(max(patchwise_fly_walkers, 1), replace = TRUE) %>%
          select(x, y)
        
        for (j in 1:max(patchwise_fly_walkers, 1)) # Simulate each fly walker
        {
          # starting position for fly_walker
          current_x <- start_positions_fly$x[j]
          current_y <- start_positions_fly$y[j]
          fly_walker_path <- data.frame(x = current_x, y = current_y)
          
          for (k in 1:20) { # Take the steps
            new_coords <- random_walk_step(current_x, current_y)
            current_x <- new_coords[1]
            current_y <- new_coords[2]
            
            # Constrain the fly_walker's position within the image bounds
            current_x <- max(min(current_x, max(landscape_data$x)), 1)
            current_y <- max(min(current_y, max(landscape_data$y)), 1)
            
            fly_walker_path <- rbind(fly_walker_path, c(current_x, current_y))
          }
          fly_walker_paths[[j]] <- fly_walker_path # Patch wise paths
        }
        
      }
      per_t_ftc_walker_paths[[i]] <- ftc_walker_paths # Store patch wise walker paths
      per_t_fly_walker_paths[[i]] <- fly_walker_paths
    }
    all_patchwise_ftc_paths[[t/5]] <- per_t_ftc_walker_paths
    all_patchwise_fly_paths[[t/5]] <- per_t_fly_walker_paths
  }
}

# Random walker plot
patch_map <- ggplot(landscape_data) +
  geom_raster(aes(x = x, y = -y, fill = as.factor(clusters))) +
  geom_path(data = do.call(rbind, all_patchwise_ftc_paths), 
            aes(x = x, y = -y), 
            color = "blue", linewidth = 1) +
  geom_path(data = do.call(rbind, all_patchwise_fly_paths), 
            aes(x = x, y = -y), 
            color = "red", linewidth = 1) +
  scale_fill_viridis(discrete = TRUE, option = "D", na.value = "white") +
  theme(legend.position = "none") +
  coord_equal()
patch_map

# patch 1 CR plot
ggplot(cr_patch, aes(x = as.numeric(row.names(cr_patch)))) +
  geom_path(aes(y = R), colour = "salmon") +
  geom_path(aes(y = C), colour = "purple") +
  labs(title = "Patch Dynamics", x = "Timestep", y = "Population size") +
  theme_minimal() +
  theme(legend.position="none")

# Iterate over patches
unique_patches <- unique(na.omit(landscape_data$clusters))