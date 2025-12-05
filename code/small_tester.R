source("code/mapping.R")
source("code/model.R")

average_host_dispersal <- 1

# load and visualize landscape
landscape <- load_landscape(
  "images/blobs.png",
  scale = 0.5, # Set scale <1 to shrink image
  cellsize = 80 # set cellsize <8 to make a finer hex grid overlay
) |> 
  basic_habitat_quality(threshold = 0.5) 
visualise_landscape(landscape)

# Generate sparse dispersal kernel
D <- compute_sparse_distance(landscape, max_host_dispersal = 1)



H <- matrix(nrow = nrow(landscape), ncol = 20)
H[,1] <- ifelse(landscape$type == "high", 1, 0)

for(t in 2:20){
  
  H[,t]
}


head(H)
