library(sf)
library(ggplot2)

create_hex_landscape <- function() {
  bbox <- st_bbox(c(xmin = 1, ymin = 1, xmax = 8, ymax = 8)) # Set bounds
  
  hex_grid <- st_make_grid(
    st_as_sfc(bbox), # Convert bounding polygon to sfc object
    cellsize = 1,        # Size of each cell
    square = FALSE,      # Cells are hexagons instead of squares
    what = "polygons"    # We want the actual polygon cells
  )
  
  # Convert to sf object and add cell IDs
  hex_sf <- st_sf(cell_id = 1:length(hex_grid), geometry = hex_grid)
  
  return(hex_sf)
}

visualise_landscape <- function(landscape){
  ggplot() +
    geom_sf(data = landscape, 
            colour = "black",     # Cell border colour
            fill = NA,            # No fill
            linewidth = 0) +    # Border line width (integers)
    labs(title = "Toy Landscape",
         x = "Longitude", 
         y = "Latitude") +
    theme_void() # Switch to theme_minimal if you want to see labels
}
