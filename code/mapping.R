library(sf) 
library(ggplot2)
library(imager)
library(dplyr)
library(patchwork)

load_landscape <- function(png, scale = 1, cellsize = 8) {
  # Load and scale
  img <- imager::load.image(png) # load image
  if (!is.numeric(scale) || scale <= 0) stop("scale must be > 0") # basic scale requirement
  if (scale != 1) img <- imager::imresize(img, scale) # rescale image if specified. Higher is bigger
  
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
  
  # Add simple hex values for habitat modeling
  hex_vals <- pts_to_hex |> # hex_vals is an sf of aggregated points with geometries
    filter(!is.na(cell_id)) |>
    group_by(cell_id) |> 
    summarize(value = mean(value, na.rm = TRUE), .groups = "drop") # get mean value of hexes
  hex_sf <- st_join(hex_sf, hex_vals["value"]) # add values to hex_sf
  
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

basic_habitat_quality <- function(landscape, threshold = 0.5){
  # Set habitat qualities
  landscape$type <- case_when( # set habitat types
    landscape$value > threshold ~ "low", # add more mean RGBO cases here
    landscape$value <= threshold ~ "high") |>
    factor(levels = c("low", "high")) 
  return(landscape)
}

dot_density_points <- function(landscape, 
                               density, 
                               dot_value = NULL, 
                               seed = NULL,
                               dot_clutter = 10) {
  if (!is.null(seed)) set.seed(seed)
  
  landscape$density <- density

  # auto-scale dot_value based on mean density
  if (is.null(dot_value)) {
    mean_density <- mean(density, na.rm = TRUE)
    if (mean_density == 0) stop("Mean density is zero, cannot scale dots automatically.")
    dot_value <- mean_density / dot_clutter
  }
  
  # compute number of dots per hex
  landscape <- landscape |>
    dplyr::mutate(n_dots = round(density / dot_value))
  
  # check for dots
  if (all(landscape$n_dots <= 0)) {
    stop("No dots to plot: all densities are zero or dot_value is too high.")
  }
  
  # keep hexes with >0 dots
  landscape_with_dots <- landscape |>
    dplyr::filter(n_dots > 0)
  
  # generate dots
  dots <- sf::st_sample(
    landscape_with_dots,
    size = landscape_with_dots$n_dots,
    type = "random"
  ) |>
    sf::st_as_sf() |>
    dplyr::mutate(dummy = 1)
  
  dots
}

visualise_landscape <- function(landscape, 
                                density = NULL, 
                                dots = NULL,
                                dotsize = 0.6,
                                show_legend = TRUE) {
  p <- ggplot() +
    geom_sf(
      data = landscape,
      aes(fill = type),
      colour = "black",
      linewidth = 0
    ) +
    scale_fill_manual(values = c(
      "low"  = "lightgrey",
      "high" = "turquoise4"
    )) +
    theme_void()
  
  if (!is.null(density)) {
    p <- p +
      geom_sf_text(
        data = landscape,
        aes(label = round(density, 2)),
        size = 3
      )
  }
  
  if (!is.null(dots)) {
    p <- p +
      geom_sf(
        data = dots,
        colour = "black",
        size = dotsize,
        alpha = 0.7
      )
  }
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  p
}
