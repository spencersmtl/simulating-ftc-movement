# Load the landscape data
landscape_data <- readRDS("data/50x50_scale2a_simple_patches_data.rds")
landscape_data$resource_occupied <- landscape_data$clusters == 1 # initial occupancy status
landscape_data$consumer_occupied <- landscape_data$clusters == 1

# Parameters
n_patches <- length(unique(na.omit(landscape_data$clusters))) # Number of patches omitting NAs

H_n_steps <- 1
P_n_steps <- 5
H_step_range <- seq(from = -20, to = 20, 5)
P_step_range <- seq(from = -20, to = 20, 5)

t <- 1               # Initial time
t_final <- 200       # number of time steps

walker_interval <- t_final/5 # how often each patch sends out walkers
H0 <- 60             # Initial host population in patch
P0 <- 20             # Initial parasitoid population in patch
r <- 3               # Intrinsic host growth rate
a <- 1               # Parasitoid attack rate
handle_t <- 1        # Handling time
e <- 0.5             # Parasitoid efficiency
d <- 1               # Parasitoid mortality
K <- 100             # Host carrying capacity

H <- matrix(nrow = t_final, ncol = n_patches)
P <- matrix(nrow = t_final, ncol = n_patches)
H[1,1] <- H0
P[1,1] <- P0

H_path <- walk_step(n_steps = H_n_steps, step_length = H_step_range, 
                    cluster_id, landscape_data)





patch_occupied <- data.frame(resource_occupied = rep(FALSE, n_patches), # occupied tag
                             consumer_occupied = rep(FALSE, n_patches))
patch_occupied$resource_occupied[1] <- TRUE
patch_occupied$consumer_occupied[1] <- TRUE

patch_death_tracker <- rep(0,n_patches) # patch extinction tracker

# list for all FTC and Fly walk lists
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


for (t in 2:t_final) {
  for (n in 1:n_patches) {
    # Update host population
    H[t, n] <- H[t - 1, n] + H_delta(H, P, r, a, K, handle_t, t, n)
    
    # Update parasitoid population
    P[t, n] <- P[t - 1, n] + P_delta(H, P, e, a, d, handle_t, t, n)
  }
}



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
  geom_path(data = do.call(rbind, all_ftc_walks[[1]][[10]]), 
            aes(x = x, y = -y), 
            color = "salmon", linewidth = 1) +
  geom_path(data = do.call(rbind, all_fly_walks[[1]][[10]]), 
            aes(x = x, y = -y), 
            color = "purple", linewidth = 1) +
  scale_fill_viridis(discrete = TRUE, option = "D", na.value = "white") +
  theme(legend.position = "none") +
  coord_equal()
patch_map

cr_both_1 <- as.data.frame(cbind(cr_patches_r[,1],cr_patches_c[,1]))
# patch 1 CR plot
ggplot(cr_both_1, aes(x = as.numeric(row.names(cr_both_1)))) +
  geom_path(aes(y = cr_both_1[,1]), colour = "salmon") +
  geom_path(aes(y = cr_both_1[,2]), colour = "purple") +
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
