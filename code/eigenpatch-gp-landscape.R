# Load functions ####
library(eigenmove)
library(gridExtra)
library(ggplotify)

# Parameters ####
{
  clustering_method = c("hclust") # options: "hclust", "DBSCAN", "OPTICS"
  n_eigs = 10 # Number of eigenvectors/eigenvalues to calculate
  tau = 40 # Interaction time
  landscape_img = "images/50x50_scale2a_simple.png" # Use any png. Images larger than 100x100 may take a very long time to compute

  # hclust parameters
  n_clust = 4 # Number of patches to define.
  
  # DBSCAN parameters
  minPts = 4 # Defines minimum size of cluster.
  eps = 20 # DBSCAN radius of neighbourhood around each point. See DBSCAN literature for more information
  eps_cl = 100 # OPTICS radius of reachability around each point. See OPTICS literature for more information
  
  # Movement model parameters: c(high, mid, low) habitat quality
  step_length = c(0.5, 0.5, 2) # Higher numbers increase distance traveled per step
  speed = c(0.5, 0.5, 2) # Higher numbers increase likelihood of cell inhabitant taking a step
  pref_strength = c(4, 2, 1) # Higher numbers increase likelihood of traveling to corresponding habitat quality
}

# Run instance ####
{
  
  gp_landscape = image_to_dataframe(landscape_img,
                                    1) # Scale image larger or smaller. Value of 1 maintains size.
  
  gp_movement_matrix = calculate_dispersal(landscape = gp_landscape,
                                           step_length = step_length,
                                           speed = speed, 
                                           pref_strength = pref_strength)
  
  gp_dist <- calculate_kinetic_distances(gp_movement_matrix,
                                         tau = tau,
                                         d = n_eigs,
                                         scale_by_density = TRUE,
                                         progress_bar = TRUE)
  
  if(clustering_method == "hclust") {
    gp_landscape <- calculate_clusters(cluster_type = "hclust",
                                       landscape = gp_landscape,
                                       out = gp_dist,
                                       min_dens = 1/nrow(gp_landscape),
                                       n_clust = n_clust,
                                       method = "ward.D2")
    
    patch_map <- ggplot(filter(gp_landscape))+
      geom_raster(aes(x=x, y=-y, fill = as.factor(clusters)))+
      labs(caption=paste0("Patch structure: n=",n_clust,
                          "\nDepth of calculation: n_eigs=",n_eigs,
                          "\nInteraction time: tau=",tau))+
      theme_void() +
      theme(legend.position = "none") +
      coord_equal()
    if(n_clust > 20) { # better colours for larger number of patches
      patch_map = patch_map + 
        scale_fill_paletteer_d("ggsci::default_igv", na.value = "grey95")}
    if (n_clust <= 20) { # better colours for smaller number of patches
      patch_map = patch_map + 
        scale_fill_paletteer_d("ggsci::category20_d3", na.value = "grey95")}
    print(patch_map)
  }
  
  if(clustering_method == "DBSCAN") {
    gp_landscape <- calculate_clusters(cluster_type = "DBSCAN",
                                       landscape = gp_landscape,
                                       out = gp_dist,
                                       min_dens = 1/nrow(gp_landscape),
                                       n_clust = n_clust,
                                       eps = eps,
                                       minPts = minPts)
    patch_map <- ggplot(filter(gp_landscape))+
      geom_raster(aes(x=x, y=-y, fill = as.factor(clusters)))+
      scale_fill_paletteer_d("ggsci::default_igv", na.value = "grey95")+
      theme_void() +
      theme(legend.position = "none") +
      coord_equal()
    print(patch_map)
  }
  
  if(clustering_method == "OPTICS") {
    gp_landscape <- calculate_clusters(cluster_type = "OPTICS",
                                       landscape = gp_landscape,
                                       out = gp_dist,
                                       min_dens = 1/nrow(gp_landscape),
                                       n_clust = n_clust,
                                       minPts = minPts,
                                       eps_cl = eps_cl)
    gp_landscape$landscape$clusters <- replace(gp_landscape$landscape$clusters, gp_landscape$landscape$clusters == 0, NA)
    
    patch_map <- ggplot(filter(gp_landscape$landscape))+
      geom_raster(aes(x=x, y=-y, fill = as.factor(clusters)))+
      scale_fill_paletteer_d("ggsci::default_igv", na.value = "grey95")+
      theme_void() +
      theme(legend.position = "none") +
      coord_equal() + 
      labs(subtitle = paste0("To adjust patch structure: 
                            \nset eps_cl based on reachability plot"))
    grid.arrange(as.grob(~plot(gp_landscape$reachability)), patch_map, ncol = 2)
  }
}
