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
R0 <- 10             # Initial FTC in patch
C0 <- 5              # Initial Fly in patch
a <- 3               # Ability of FTC to avoid flies
A <- 50              # FTC half-saturation constant ???
b <- 4               # Fly consumption ability
B <- 33              # Fly half-saturation constant
e <- 0.5             # Extinction rate of the fly
K <- 100             # FTC Carrying capacity

# Initialize objects before loop
cr_patches <- vector("list", n_patches) # cr series all patches
cr_patch1 <- data.frame(Resource = numeric(t_final), # cr series patch 1
                        Consumer = numeric(t_final))
cr_patch1$Resource[1] = R0 # Set initials
cr_patch1$Consumer[1] = C0
cr_patches[[1]] <- cr_patch1

all_ftc_paths <- vector("list", t_final/5 - 2) # FTC paths for plotting
all_fly_paths <- vector("list", t_final/5 - 2) # Fly paths for plotting
}

# Time series loop containing CR and random walking
for(t in 1:(t_final-1)) {
  # Patch CR dynamics ####
  
  # Subset data for the current patch
  current_cluster <- filter(landscape_data, clusters==t)
  
  # Calculate current values
  R1 <- cr_patch1$Resource # for readability
  C1 <- cr_patch1$Consumer # for readability
  R1[t+1] = R1[t] + a * R1[t] * (1 - R1[t] / K) - b * R1[t] / (A + R1[t]) * C1[t] # Resource AKA FTC dynamics
  C1[t+1] = C1[t] + e * ((R1[t] / (A + R1[t])) - (B / (B + A))) * C1[t] # Consumer AKA Fly dynamics
  cr_patch1$Resource[t+1] <- R1[t+1]
  cr_patch1$Consumer[t+1] <- C1[t+1]
  
  # Random walking ####
  if(t %% 5 == 0) {
    num_ftc_walkers <- floor(R1[t+1]/10) # Number of FTC walkers is proportional to patch population size
    num_fly_walkers <- floor(C1[t+1]/10)
    ftc_walker_paths <- vector("list", num_ftc_walkers)
    fly_walker_paths <- vector("list", num_fly_walkers)
    
    # Find starting positions within the blue patch
    start_positions_ftc <- landscape_data %>% # FTC starting positions
      filter(resource_occupied) %>%
      sample_n(num_ftc_walkers, replace = TRUE) %>%
      select(x, y)
    start_positions_fly <- landscape_data %>% # Fly starting positions
      filter(consumer_occupied) %>%
      sample_n(max(num_fly_walkers, 1), replace = TRUE) %>%
      select(x, y)
    
    # Simulate ftc_walkers
    for (i in 1:num_ftc_walkers) {
      # starting position for ftc_walker
      current_x <- start_positions_ftc$x[i]
      current_y <- start_positions_ftc$y[i]
      ftc_walker_path <- data.frame(x = current_x, y = current_y)
      
      for (j in 1:10) { # Take the steps
        new_coords <- random_walk_step(current_x, current_y)
        current_x <- new_coords[1]
        current_y <- new_coords[2]
        
        # Constrain the ftc_walker's position within the image bounds
        current_x <- max(min(current_x, max(landscape_data$x)), 1)
        current_y <- max(min(current_y, max(landscape_data$y)), 1)
        
        ftc_walker_path <- rbind(ftc_walker_path, c(current_x, current_y))
      }
      
      ftc_walker_paths[[i]] <- ftc_walker_path
    }
    all_ftc_paths[[t/5]] <- ftc_walker_paths[[i]]
    
    # Simulate fly_walkers
    for (i in 1:max(num_fly_walkers, 1)) {
      # starting position for fly_walker
      current_x <- start_positions_fly$x[i]
      current_y <- start_positions_fly$y[i]
      fly_walker_path <- data.frame(x = current_x, y = current_y)
      
      for (j in 1:20) { # Take the steps
        new_coords <- random_walk_step(current_x, current_y)
        current_x <- new_coords[1]
        current_y <- new_coords[2]
        
        # Constrain the fly_walker's position within the image bounds
        current_x <- max(min(current_x, max(landscape_data$x)), 1)
        current_y <- max(min(current_y, max(landscape_data$y)), 1)
        
        fly_walker_path <- rbind(fly_walker_path, c(current_x, current_y))
      }
      
      fly_walker_paths[[i]] <- fly_walker_path
    }
    all_fly_paths[[t/5]] <- fly_walker_paths[[i]]
  }
}

# Random walker plot
patch_map <- ggplot(landscape_data) +
  geom_raster(aes(x = x, y = -y, fill = as.factor(clusters))) +
  geom_path(data = do.call(rbind, all_ftc_paths), 
            aes(x = x, y = -y), 
            color = "blue", linewidth = 1) +
  geom_path(data = do.call(rbind, all_fly_paths), 
            aes(x = x, y = -y), 
            color = "red", linewidth = 1) +
  scale_fill_viridis(discrete = TRUE, option = "D", na.value = "white") +
  theme(legend.position = "none") +
  coord_equal()
patch_map

# patch 1 CR plot
ggplot(cr_patch1, aes(x = as.numeric(row.names(cr_patch1)))) +
  geom_path(aes(y = Resource), colour = "salmon") +
  geom_path(aes(y = Consumer), colour = "purple") +
  labs(title = "Patch Dynamics", x = "Timestep", y = "Population size") +
  theme_minimal() +
  theme(legend.position="none")

# Iterate over patches
unique_patches <- unique(na.omit(landscape_data$clusters))