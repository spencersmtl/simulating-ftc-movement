# Load the landscape data
landscape_data <- readRDS("data/50x50_scale2a_simple_patches_data.rds")
landscape_data$resource_occupied <- landscape_data$clusters == 1 # initial occupancy status
landscape_data$consumer_occupied <- landscape_data$clusters == 1

# General Parameters
n_patches <- length(unique(na.omit(landscape_data$clusters))) # Number of patches omitting NAs
t <- 1               # Initial time
t_final <- 200       # number of time steps
walker_interval <- t_final/5 # how often each patch sends out walkers
R0 <- 50             # Initial FTC in patch
C0 <- 10             # Initial Fly in patch
a <- 2.8               # FTC Growth rate
A <- 52              # FTC half-saturation constant ???
b <- 6               # Fly consumption ability
m <- 0.4              # Mortality
e <- 0.5             # Fly extinction rate
K <- 100               # FTC Carrying capacity

cr_patches_r <- matrix(0, ncol = n_patches, nrow = t_final)
cr_patches_c <- matrix(0, ncol = n_patches, nrow = t_final)
cr_patches_r[1,1] = R0 # Set initials
cr_patches_c[1,1] = C0

patch_occupied <- data.frame(resource_occupied = rep(FALSE, n_patches), # occupied tag
                             consumer_occupied = rep(FALSE, n_patches))
patch_occupied$resource_occupied[1] <- TRUE
patch_occupied$consumer_occupied[1] <- TRUE

patch_death_tracker <- rep(0,n_patches) # patch extinction tracker

for(t in 2:(t_final)) # Time series loop containing CR and random walking
{
  for(i in 1:n_patches) # Compute CR for each patch
  {
    
    if(patch_death_tracker[i]>20) # Patch death!
    { 
      landscape_data$clusters[landscape_data$clusters == i] <- 0
      patch_occupied[i,] <- FALSE
    }
    
    if(patch_occupied[i,1] & patch_occupied[i,2]) # Occupied by both
    { 
      cr_patches_r[t,i] = pmax(
        cr_patches_r[t-1,i] + 
          a * cr_patches_r[t-1,i] * (1 - cr_patches_r[t-1,i] / K) - 
          b * cr_patches_r[t-1,i] / (A + cr_patches_r[t-1,i]) * cr_patches_c[t-1,i],0
      )
      cr_patches_c[t,i] = pmax(
        cr_patches_c[t-1,i] + 
          e * cr_patches_c[t-1,i] * 
          ((cr_patches_r[t-1,i] / (A + cr_patches_r[t-1,i])) - m),0
      )
    }
    
    if(patch_occupied[i,1] & patch_occupied[i,2] == FALSE) # Occupied by just FTC
    { 
      cr_patches_r[t,i] = pmax(
        cr_patches_r[t-1,i] + 
          a * cr_patches_r[t-1,i] * (1 - cr_patches_r[t-1,i] / K),0)
    }

    if(cr_patches_r[t,i]>80)# Track if patch is going extinct
    { 
      patch_death_tracker[i] <- patch_death_tracker[i]+1} else {patch_death_tracker[i]=0}
  }
}

cr_both_1 <- as.data.frame(cbind(cr_patches_r[,1],cr_patches_c[,1]))
# patch 1 CR plot
ggplot(cr_both_1, aes(x = as.numeric(row.names(cr_both_1)))) +
  geom_path(aes(y = cr_both_1[,1]), colour = "salmon") +
  geom_path(aes(y = cr_both_1[,2]), colour = "purple") +
  labs(title = "Patch Dynamics", x = "Timestep", y = "Population size") +
  theme_minimal() +
  theme(legend.position="none")
