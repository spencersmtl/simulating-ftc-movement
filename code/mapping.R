library(sf) 
library(ggplot2)
library(imager)
library(dplyr)

load_landscape <- function(png, scale = 1, cellsize = 8) {
  # Load and scale
  img <- imager::load.image(png) # load image
  if (!is.numeric(scale) || scale <= 0) stop("scale must be > 0") # basic scale requirement
  if (scale != 1) img <- imager::imresize(img, scale) # rescale image if specified. Higher is bigger
  
  # Turn into dataframe and process
  df_px <- as.data.frame(img) |> 
    group_by(x, y) |> 
    summarize(value = mean(value, na.rm = TRUE), # get mean RGBO of each pixel
              .groups = "drop") |>
    mutate(y = max(y) - y + 1) # flip image right side up
  
  # Turn into hex grid
  pts <- st_as_sf(df_px, coords = c("x", "y"), crs = NA) # points as sf pixel cartesian coordinates
  bbox_sfc <- st_as_sfc(st_bbox(pts)) # pixel bounding box
  bbox_exp <- st_buffer(bbox_sfc, dist = 1e-5) # tiny buffer to avoid edge issues
  hex_grid <- st_make_grid(   # Make a hex grid over the bounding box
    bbox_exp, 
    cellsize = cellsize,
    square = FALSE,
    what = "polygons")
  inside <- st_within(st_centroid(hex_grid), bbox_sfc, sparse = FALSE)   # specify only hexes fully inside the image area
  hex_grid_full <- hex_grid[inside, ] # keep specified hexes
  hex_sf <- st_sf(cell_id = seq_along(hex_grid_full), geometry = hex_grid_full) # Convert to sf and assign cell IDs
  pts_to_hex <- st_join(pts, hex_sf, join = st_intersects) # assign all pixels to hexes
  
  # Make sure scale/cellsize is appropriate
  counts <- tabulate(pts_to_hex$cell_id, nbins = nrow(hex_sf)) # Count how many pixels were assigned to each hex
  empty_frac <- mean(counts == 0) # Count how many empty hexes there are
  if (empty_frac > 0.01) {stop(sprintf(
    "cellsize too small for this image: %.1f%% of hexes contain no pixels. Increase cellsize or increase scale factor.",
    empty_frac * 100))}
  
  # Quantify hexes based on pixel values for optional modeling use later (e.g., habitat quality)
  hex_vals <- pts_to_hex |> # hex_vals is an sf of aggregated points with geometries
    group_by(cell_id) |> 
    summarize(value = mean(value, na.rm = TRUE), .groups = "drop") # get mean value of hexes
  hex_vals_df <- sf::st_set_geometry(hex_vals, NULL)   # plain df to use left_join
  hex_sf <- dplyr::left_join(hex_sf, hex_vals_df, by = "cell_id") # add values to hex_sf
  hex_sf$value[is.na(hex_sf$value)] <- 1   # set NAs to base habitat type
  
  # Set habitat qualities
  hex_sf$type <- dplyr::case_when( # set habitat types
    hex_sf$value > 0.5 ~ "low", # add more mean RGBO cases here
    hex_sf$value <= 0.5 ~ "high") |>
    factor(levels = c("low", "high")) 
  
  # Compute centroids of cells
  cents <- st_coordinates(st_centroid(hex_sf))
  hex_sf$centroid_x <- cents[,1]
  hex_sf$centroid_y <- cents[,2]
  
  # Compute raw spacing between any adjacent hexes
  neighbors <- st_touches(hex_sf)
  i <- 1
  j <- neighbors[[i]][1]
  raw_adj <- sqrt(sum((cents[i, ] - cents[j, ])^2))
  
  # Normalize centroid coordinates so adjacent spacing = 1
  hex_sf$centroid_x <- hex_sf$centroid_x / raw_adj
  hex_sf$centroid_y <- hex_sf$centroid_y / raw_adj
  
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
