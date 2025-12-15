source("code/mapping.R")
source("code/model.R")

# load and visualize landscape
landscape <- load_landscape(
  "images/nb-scape-bw-thick.png",
  scale = 0.5, # Set scale <1 to shrink image
  cellsize = 8 # set cellsize <8 to make a finer hex grid overlay
) |> 
  basic_habitat_quality(threshold = 0.5) 
visualise_landscape(landscape)

# Generate sparse dispersal kernel
D <- compute_sparse_distance(
  landscape, 
  avg_host_dispersal = 4, 
  max_host_dispersal_mult = 1.1, 
  normalize_rows = TRUE)



H <- matrix(nrow = nrow(landscape), ncol = 20)
H[,1] <- ifelse(landscape$type == "high", 1, 0)

for(t in 2:20){
  
  H[,t]
}


head(H)
