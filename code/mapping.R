library(sf) 
library(ggplot2)
library(imager)
library(dplyr)

load_landscape <- function(png, scale = 1, cellsize = 8) {
  # Load and scale
  img <- imager::load.image(png) # load image
  if (!is.numeric(scale) || scale <= 0) stop("scale must be > 0") # basic scale requirement
  if (scale != 1) img <- imager::imresize(img, scale) # rescale image if specified. Higher is bigger
  browser()
  
  # Turn into dataframe and process
  pixel_df <- as.data.frame(img) |> 
    group_by(x, y) |> 
    summarize(value = mean(value, na.rm = TRUE), # get mean RGBO of each pixel
              .groups = "drop") |>
    mutate(y = max(y) - y + 1) # flip image right side up
  
  # Turn into hex grid
  point_geometry <- st_as_sf(pixel_df, coords = c("x", "y")) # points as sf pixel cartesian coordinates
  bbox_sfc <- st_as_sfc(st_bbox(point_geometry)) # pixel bounding box
  bbox_exp <- st_buffer(bbox_sfc, dist = 1e-5) # tiny buffer to avoid edge issues
  hex_grid <- st_make_grid(   # Make a hex grid over the bounding box
    bbox_exp, 
    cellsize = cellsize,
    square = FALSE,
    what = "polygons")
  inside <- st_within(st_centroid(hex_grid), bbox_sfc, sparse = FALSE)   # specify only hexes fully inside the image area
  hex_grid <- hex_grid[inside, ] # keep specified hexes
  hex_sf <- st_sf(cell_id = seq_along(hex_grid), geometry = hex_grid) # Convert to sf and assign cell IDs
  pts_to_hex <- st_join(point_geometry, hex_sf, join = st_intersects) # assign all pixels to hexes
  
  # Add centroids
  centroids <- st_coordinates(st_centroid(hex_sf))
  hex_sf$centroid_x <- centroids[,1]
  hex_sf$centroid_y <- centroids[,2]
  
  # Make sure scale/cellsize is appropriate
  na_count <- sum(is.na(pts_to_hex$cell_id)) # Count how many pixels were not assigned to a hex
  if (na_count/nrow(pts_to_hex) > 0.2) {warning(paste(
    "cellsize is large for this image: ", 
    round(na_count/nrow(pts_to_hex)*100, 1), 
    "% of pixels not assigned to hexes, consider decreasing cellsize."))}
  pixel_per_hex_count <- tabulate(pts_to_hex$cell_id, nbins = nrow(hex_sf)) # Count how many pixels were assigned to each hex
  empty_frac <- mean(pixel_per_hex_count == 0) # Count how many empty hexes there are
  if (empty_frac > 0.01) {stop(paste(
    "cellsize too small for this image: ",  
    empty_frac * 100, 
    "% of hexes contain no pixels. Increase cellsize or increase scale factor."))}
  
  return(hex_sf)
}



create_hex_landscape <- function() {
  bounding_box <- st_bbox(c(xmin = 1, ymin = 1, xmax = 8, ymax = 8)) # Set bounds
  
  hex_grid <- st_make_grid(
    st_as_sfc(bounding_box), # Convert bounding polygon to sfc object
    cellsize = 1,        # Size of each cell
    square = FALSE,      # Cells are hexagons instead of squares
    what = "polygons"    # We want the actual polygon cells
  )
  
  # Convert to sf object and add cell IDs
  hex_sf <- st_sf(cell_id = seq_along(hex_grid), geometry = hex_grid)
  
  return(hex_sf)
}

visualise_landscape <- function(landscape){
  ggplot() +
    geom_sf(data = landscape, 
            colour = "black",     # Cell border colour
            aes(fill = type),            # No fill
            linewidth = 0) +    # Border line width (integers)
    scale_fill_manual(values = c("low" = "lightgrey",
                                 "high" = "darkgreen")) +
    labs(title = "Toy Landscape",
         x = "Longitude", 
         y = "Latitude") +
    theme_void() # Switch to theme_minimal if you want to see labels
}
