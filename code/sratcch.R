# Load necessary libraries
library(png)
library(viridis)
library(tidyverse)

# Load the landscape data
landscape_data <- readRDS("data/50x50_scale2a_simple_patches_data.rds")

# Initialize occupancy status of patch
landscape_data$resource_occupied <- landscape_data$clusters == 1
landscape_data$consumer_occupied <- landscape_data$clusters == 1

# Function to perform random walks for both FTC and fly
random_walk_step <- function(current_x, current_y) {
  dx <- sample(c(-8, -4, 0, 4, 8), 1) # x direction
  dy <- sample(c(-8, -4, 0, 4, 8), 1) # y direction
  return(c(current_x + dx, current_y + dy))
}

# Parameters
n_patches <- length(unique(na.omit(landscape_data$clusters))) # Number of patches omitting NAs
ftc_steps <- 10
fly_steps <- 20
t <- 1               # Initial time
t_final <- 200       # number of time steps
walker_interval <- t_final/5 # how often each patch sends out walkers
R0 <- 60             # Initial FTC in patch
C0 <- 20              # Initial Fly in patch
a <- 3               # Ability of FTC to avoid flies
A <- 50              # FTC half-saturation constant ???
b <- 5               # Fly consumption ability
B <- 35              # Fly half-saturation constant
e <- 0.5             # Extinction rate of the fly
K <- 100               # FTC Carrying capacity

cr_patches_r <- matrix(0, ncol = n_patches, nrow = t_final)
cr_patches_c <- matrix(0, ncol = n_patches, nrow = t_final)
cr_patches_r[1,1] = R0 # Set initials
cr_patches_c[1,1] = C0

patch_occupied <- data.frame(resource_occupied = rep(FALSE, n_patches), 
                             consumer_occupied = rep(FALSE, n_patches))
patch_occupied$resource_occupied[1] <- TRUE
patch_occupied$consumer_occupied[1] <- TRUE#landscape_data[landscape_data$clusters == i,][1,]$resource_occupied

patch_death_tracker <- rep(0,n_patches)

per_t_ftc_walker_paths <- vector("list", 1) # All paths in all patches per time step
per_t_fly_walker_paths <- vector("list", 1)

all_patchwise_ftc_paths <- vector("list", t_final/5 - 2) # All paths in all patches for all time steps
all_patchwise_fly_paths <- vector("list", t_final/5 - 2) 

all_ftc_paths <- vector("list", n_patches)

# Create the main list for all FTC and Fly walks
all_ftc_walks <- vector("list", n_patches)
all_fly_walks <- vector("list", n_patches)
for (n in 1:n_patches) {
  patch_walks <- vector("list", 200)
  names(patch_walks) <- paste0("walkers at time ", 1:200)
  all_ftc_walks[[n]] <- patch_walks
  all_fly_walks[[n]] <- patch_walks
}
names(all_ftc_walks) <- paste0("Patch_", 1:n_patches)
names(all_fly_walks) <- paste0("Patch_", 1:n_patches)




for(t in 2:(t_final)) # Time series loop containing CR and random walking
{
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
  if(t %% 5 == 0) # every 5 time steps...
  {
    for(i in 1:n_patches) # Compute random walks for each patch that is occupied
    {
      if(patch_occupied[i,1]) # FTC
      {
        num_ftc_walkers <- pmax(floor(cr_patches_r[t,i]/10),1) # Number of walkers is proportional to patch population size
        ftc_walkers <- vector("list", num_ftc_walkers) # Per patch, per step
        
        start_positions_ftc <- landscape_data %>% # FTC starting positions
          filter(resource_occupied, clusters == i) %>%
          sample_n(num_ftc_walkers, replace = TRUE) %>%
          select(x, y)
        
        for (j in 1:num_ftc_walkers) # Simulate each ftc walker
        {
          # starting position for ftc_walker
          current_x <- start_positions_ftc$x[j]
          current_y <- start_positions_ftc$y[j]
          ftc_walker_path <- data.frame(x = current_x, y = current_y)
          
          for (k in 1:ftc_steps) { # Take the steps
            new_coords <- random_walk_step(current_x, current_y)
            current_x <- new_coords[1]
            current_y <- new_coords[2]
            
            # Constrain the ftc_walker's position within the image bounds
            current_x <- max(min(current_x, max(landscape_data$x)), 1)
            current_y <- max(min(current_y, max(landscape_data$y)), 1)
            
            ftc_walker_path <- rbind(ftc_walker_path, c(current_x, current_y))
          }
          ftc_walkers[[j]] <- ftc_walker_path # Each walker gets a list entry
        }
        # store list of individual walker dataframes in a new list! Per patch
        all_ftc_walks[[i]][[t]] <- ftc_walkers
      }
    
      if(patch_occupied[i,2]) # FLY
      {
        num_fly_walkers <- pmax(floor(cr_patches_c[t,i]/10),1) # Number of walkers is proportional to patch population size
        fly_walkers <- vector("list", num_fly_walkers) # Per patch, per step
        
        start_positions_fly <- landscape_data %>% # fly starting positions
          filter(consumer_occupied, clusters == i) %>%
          sample_n(num_fly_walkers, replace = TRUE) %>%
          select(x, y)
        
        for (j in 1:num_fly_walkers) # Simulate each fly walker
        {
          # starting position for fly_walker
          current_x <- start_positions_fly$x[j]
          current_y <- start_positions_fly$y[j]
          fly_walker_path <- data.frame(x = current_x, y = current_y)
          
          for (k in 1:fly_steps) { # Take the steps
            new_coords <- random_walk_step(current_x, current_y)
            current_x <- new_coords[1]
            current_y <- new_coords[2]
            
            # Constrain the fly_walker's position within the image bounds
            current_x <- max(min(current_x, max(landscape_data$x)), 1)
            current_y <- max(min(current_y, max(landscape_data$y)), 1)
            
            fly_walker_path <- rbind(fly_walker_path, c(current_x, current_y))
          }
          fly_walkers[[j]] <- fly_walker_path # Each walker gets a list entry
        }
        # store list of individual walker dataframes in a new list! Per patch
        all_fly_walks[[i]][[t]] <- fly_walkers
      }
    }
  }
}

# Random walker plot
patch_map <- ggplot(landscape_data) +
  geom_raster(aes(x = x, y = -y, fill = as.factor(clusters))) +
  geom_path(data = do.call(rbind, all_ftc_walks[[1]][[5]]), 
            aes(x = x, y = -y), 
            color = "blue", linewidth = 1) +
  geom_path(data = do.call(rbind, all_fly_walks[[1]][[5]]), 
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

# Initial plot for fun
patch_map <- ggplot(landscape_data) +
  geom_raster(aes(x = x, y = -y, fill = as.factor(clusters))) +
  scale_fill_viridis(discrete = TRUE, option = "D", na.value = "white") +  # "D" is just one of the available options
  theme(legend.position = "none") +
  coord_equal()
patch_map

# Iterate over patches
unique_patches <- unique(na.omit(landscape_data$clusters))