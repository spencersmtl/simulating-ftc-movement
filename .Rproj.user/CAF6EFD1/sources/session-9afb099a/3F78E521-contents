# Load necessary libraries
library(png)
library(viridis)
library(tidyverse)

# Load the landscape data
landscape_data <- readRDS("data/50x50_scale2a_simple_patches_data.rds")

patch_map <- ggplot(landscape_data) +
  geom_raster(aes(x = x, y = -y, fill = as.factor(clusters))) +
  scale_fill_viridis(discrete = TRUE, option = "D", na.value = "white") +  # "D" is just one of the available options
  theme(legend.position = "none") +
  coord_equal()

# Initialize occupancy status of patch
landscape_data$resource_occupied <- landscape_data$clusters == 1
landscape_data$consumer_occupied <- landscape_data$clusters == 1

# Parameters (example values; adjust as needed)

{
n_patches = length(unique(na.omit(landscape_data$clusters)))
t <- 1               # Initial time
t_final <- 200       # number of time steps
R0 <- 10
C0 <- 5
a <- 3             # Avoidability of the resource
A <- 50              # Resource half-saturation constant
b <- 4             # Consumer consumption ability
B <- 33              # Consumer half-saturation constant
e <- 0.5             # Extinction rate of the consumer
K <- 100             # Carrying capacity
}

cr_patch1 <- data.frame(Resource = numeric(t_final), # consumer-resource dynamics time series for patch 1
                        Consumer = numeric(t_final))
cr_patch1$Resource[1] = R0
cr_patch1$Consumer[1] = C0

for(t in 1:(t_final-1)) {
  # Subset data for the current patch
  current_cluster <- filter(landscape_data, clusters==t)
  
  # Calculate current values
  R1 <- cr_patch1$Resource # resource time series data for patch 1. Done for readability
  C1 <- cr_patch1$Consumer # consumer time series data for patch 1. Done for readability
  
  R1[t+1] = R1[t] + a * R1[t] * (1 - R1[t] / K) - b * R1[t] / (A + R1[t]) * C1[t]
  C1[t+1] = C1[t] + e * ((R1[t] / (A + R1[t])) - (B / (B + A))) * C1[t]
  cr_patch1$Resource[t+1] <- R1[t+1]
  cr_patch1$Consumer[t+1] <- C1[t+1]
  
  
}

# Create the plot
ggplot(cr_patch1, aes(x = as.numeric(row.names(cr_patch1)))) +
  geom_path(aes(y = Resource), colour = "salmon") +
  geom_path(aes(y = Consumer), colour = "purple") +
  labs(title = "Patch Dynamics", x = "Timestep", y = "Population size") +
  theme_minimal() +
  theme(legend.position="none")

# Iterate over patches
unique_patches <- unique(na.omit(landscape_data$clusters))




# Function to perform random walk
random_walk_step <- function(current_x, current_y) {
  # Random step in x and y directions
  dx <- sample(c(-5, -2, 0, 2, 5), 1)
  dy <- sample(c(-5, -2, 0, 2, 5), 1)
  
  return(c(current_x + dx, current_y + dy))
}

# Add random walkers
num_walkers <- 1
walker_paths <- vector("list", num_walkers)

# Find starting positions within the blue patch
start_positions <- landscape_data %>%
  filter(resource_occupied) %>%
  sample_n(num_walkers, replace = TRUE) %>%
  select(x, y)

# Simulate random walks
for (i in 1:num_walkers) {
  current_x <- start_positions$x[i]
  current_y <- start_positions$y[i]
  
  walker_path <- data.frame(x = current_x, y = current_y)
  
  for (j in 1:5) { # Number of steps
    new_coords <- random_walk_step(current_x, current_y)
    current_x <- new_coords[1]
    current_y <- new_coords[2]
    
    # Constrain the walker's position within the image bounds
    current_x <- max(min(current_x, max(landscape_data$x)), 1)
    current_y <- max(min(current_y, max(landscape_data$y)), 1)
    
    walker_path <- rbind(walker_path, c(current_x, current_y))
  }
  
  walker_paths[[i]] <- walker_path
}

# Overlay random walker paths on the image
patch_map <- ggplot(landscape_data) +
  geom_raster(aes(x = x, y = -y, fill = as.factor(clusters))) +
  geom_path(data = do.call(rbind, walker_paths), 
            aes(x = x, y = -y), 
            color = "blue", linewidth = 1) +
  scale_fill_viridis(discrete = TRUE, option = "D", na.value = "white") +
  theme(legend.position = "none") +
  coord_equal()
patch_map
