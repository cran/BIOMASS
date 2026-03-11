## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, echo = TRUE,
  comment = "#>", fig.align = "center"
)
require(BIOMASS)
require(knitr)
require(ggplot2)

## -----------------------------------------------------------------------------
data("NouraguesPlot201")

kable(head(NouraguesPlot201), digits = 5, row.names = FALSE, caption = "Head of NouraguesPlot201")

## -----------------------------------------------------------------------------
data("NouraguesCoords")

kable(head(NouraguesCoords), digits = 5, row.names = FALSE, caption = "Head of NouraguesCoords")

## -----------------------------------------------------------------------------
data("NouraguesTrees")

kable(head(NouraguesTrees), digits = 3, row.names = FALSE, caption = "Head of the table trees")

## ----check_plot_trust_GPS, dpi=200, message=FALSE-----------------------------
check_plot_trust_GPS <- check_plot_coord(
  corner_data = NouraguesPlot201,
  longlat = c("Long", "Lat"),  # or proj_coord = c("Xutm", "Yutm"), 
  rel_coord = c("Xfield", "Yfield"),
  trust_GPS_corners = T,
  draw_plot = TRUE,
  max_dist = 10, rm_outliers = TRUE)

## ----dpi=200------------------------------------------------------------------
degraded_corner_coord <- NouraguesPlot201[c(1:2,11:12,21:22,31:32),]

check_plot_trust_field <- check_plot_coord(
  corner_data = degraded_corner_coord,
  longlat = c("Long", "Lat"),  # or proj_coord = c("Xutm", "Yutm"), 
  rel_coord = c("Xfield", "Yfield"),
  trust_GPS_corners = FALSE,
  draw_plot = TRUE, rm_outliers = FALSE)

## -----------------------------------------------------------------------------
kable(check_plot_trust_GPS$corner_coord, row.names = FALSE, caption = "Reference corner coordinates")

## ----eval=FALSE---------------------------------------------------------------
# sf::st_write(check_plot_trust_GPS$polygon, "your_directory/plot201.shp")

## -----------------------------------------------------------------------------
plot201Trees <- NouraguesTrees[NouraguesTrees$Plot==201,]

check_plot_trust_GPS <- check_plot_coord(
  corner_data = NouraguesPlot201,
  longlat = c("Long", "Lat"), rel_coord = c("Xfield", "Yfield"),
  trust_GPS_corners = TRUE,
  tree_data = plot201Trees, tree_coords = c("Xfield","Yfield"))

## -----------------------------------------------------------------------------
kable(head(check_plot_trust_GPS$tree_data[,-c(5,6,7)]), digits = 3, row.names = FALSE, caption = "Head of the $tree_data output")

## -----------------------------------------------------------------------------
plot_to_change <- check_plot_trust_GPS$plot_design
plot_to_change <- plot_to_change + ggtitle("A custom title")
plot_to_change

## -----------------------------------------------------------------------------
# Load internal CHM raster
nouraguesRaster <- terra::rast(system.file("extdata", "NouraguesRaster.tif",package = "BIOMASS", mustWork = TRUE))

check_plot_trust_GPS <- check_plot_coord(
  corner_data = NouraguesPlot201,
  longlat = c("Long", "Lat"), rel_coord = c("Xfield", "Yfield"),
  trust_GPS_corners = TRUE,
  tree_data = plot201Trees, tree_coords = c("Xfield","Yfield"),
  prop_tree = "D", threshold_tree = 20, # Display tree diameters >= 20
  ref_raster = nouraguesRaster )

## -----------------------------------------------------------------------------
multiple_checks <- check_plot_coord(
  corner_data = NouraguesCoords, # NouraguesCoords contains 4 plots
  proj_coord = c("Xutm", "Yutm"), rel_coord = c("Xfield", "Yfield"),
  trust_GPS_corners = TRUE, 
  plot_ID = "Plot",
  tree_data = NouraguesTrees, tree_coords = c("Xfield","Yfield"), 
  prop_tree = "D", tree_plot_ID = "Plot",
  ref_raster = nouraguesRaster)

## ----divide_plot--------------------------------------------------------------
subplots <- divide_plot(
  corner_data = check_plot_trust_GPS$corner_coord,
  rel_coord = c("x_rel","y_rel"),
  proj_coord = c("x_proj","y_proj"),
  grid_size = 25 # or c(25,25)
  )

kable(head(subplots$sub_corner_coord, 10), digits = 1, row.names = FALSE, caption = "Head of the divide_plot()$sub_corner_coord output.")

## -----------------------------------------------------------------------------
subplots <- divide_plot(
  corner_data = check_plot_trust_GPS$corner_coord, 
  rel_coord = c("x_rel","y_rel"), proj_coord = c("x_proj","y_proj"),
  grid_size = c(40,45), 
  grid_tol = 0.3, # by default =0.1, ie, if more than 10% of the plot is not covered by the grid, it will returned an error
  origin = c(10,5)
  )

## ----imperfect_cuts_visualisation, echo=FALSE, fig.show="hold", out.width="50%", warning=FALSE----
non_centred <- divide_plot(
  corner_data = check_plot_trust_GPS$corner_coord, 
  rel_coord = c("x_rel","y_rel"), proj_coord = c("x_proj","y_proj"),
  grid_size = c(40,45), 
  grid_tol = 0.3)

ggplot(data = subplots$sub_corner_coord, mapping = aes(x=x_proj, y=y_proj)) + 
  geom_point(data = check_plot_trust_GPS$corner_coord, 
             mapping = aes(x=x_proj, y=y_proj),
             shape = 15, size = 2) + 
  geom_point(col="red") + 
  coord_equal() + 
  theme_bw() + 
  labs(title = "subplot divisions with origin at (10,5)")

ggplot(data = non_centred$sub_corner_coord, mapping = aes(x=x_proj, y=y_proj)) + 
  geom_point(data = check_plot_trust_GPS$corner_coord, 
             mapping = aes(x=x_proj, y=y_proj),
             shape = 15, size = 2.5) + 
  geom_point(col="red") + 
  coord_equal() + 
  theme_bw() + 
  labs(title = "subplot divisions with origin at (0,0) by default")


## ----divide_plot_trees--------------------------------------------------------
# Add AGB predictions (calculated in Vignette BIOMASS) to plot201Trees
AGB_data <- readRDS("saved_data/NouraguesTreesAGB.rds")
plot201Trees <- merge(plot201Trees , AGB_data[c("Xfield","Yfield","D","AGB")], sort=FALSE)

subplots <- divide_plot(
  corner_data = check_plot_trust_GPS$corner_coord, 
  rel_coord = c("x_rel","y_rel"),
  proj_coord = c("x_proj","y_proj"),
  grid_size = 25, # or c(25,25)
  tree_data = plot201Trees, tree_coords = c("Xfield","Yfield")
  )

## -----------------------------------------------------------------------------
kable(head(subplots$tree_data[,-c(2,3,4)]), digits = 1, row.names = FALSE, caption = "Head of the divide_plot()$tree_data returns")

## ----divide_multiple_plots----------------------------------------------------
multiple_subplots <- divide_plot(
  corner_data = NouraguesCoords,
  rel_coord = c("Xfield","Yfield"), proj_coord = c("Xutm","Yutm"), corner_plot_ID = "Plot",
  grid_size = 25, 
  tree_data = NouraguesTrees, tree_coords = c("Xfield","Yfield"), tree_plot_ID = "Plot"
)

## ----divide_sd_coord----------------------------------------------------------
sd_coord_subplots <- divide_plot(
  corner_data = check_plot_trust_GPS$corner_coord,
  rel_coord = c("x_rel","y_rel"),
  proj_coord = c("x_proj","y_proj"),
  grid_size = 25, # or c(25,25)
  tree_data = plot201Trees, tree_coords = c("Xfield","Yfield"),
  sd_coord = check_plot_trust_GPS$sd_coord, n = 50
  )

## ----subplot_summary----------------------------------------------------------
subplot_metric <- subplot_summary(
  subplots = subplots,
  value = "AGB", # AGB was added before applying divide_plot()
  per_ha = TRUE) 

## ----subplot_summary_quantile-------------------------------------------------
subplot_metric <- subplot_summary(
  subplots = subplots,
  value = "AGB",
  fun = quantile, probs = 0.5, # yes, it is the median
  per_ha = FALSE)

## ----eval=FALSE---------------------------------------------------------------
# # Set the CRS of the polygons
# subplot_polygons <- sf::st_set_crs(
#   subplot_metric$polygon ,
#   value = "EPSG:2972") # EPSG:2972 (corresponding to UTM Zone 22N) is the UTM coordinate system of Nouragues
# 
# # Save the polygons in a shapefile
# sf::st_write(subplot_polygons, "your_directory/subplots_201.shp")

## ----subplot_summary_display_trees--------------------------------------------
multiple_subplot_metric <- subplot_summary(
  subplots = multiple_subplots, draw_plot = FALSE,
  value = "D", fun = mean, per_ha = FALSE)

## ----include=FALSE------------------------------------------------------------
error_prop <- readRDS("saved_data/error_prop.rds")

## -----------------------------------------------------------------------------
subplot_AGBD <- subplot_summary(
  subplots = sd_coord_subplots,
  AGB_simu = error_prop$AGB_simu # error_prop has been created in the previous vignette
)

## -----------------------------------------------------------------------------
raster_summary <- subplot_summary(
  subplots = subplots,
  ref_raster = nouraguesRaster, raster_fun = median)

## -----------------------------------------------------------------------------
subplot_metric <- subplot_summary(
  subplots = sd_coord_subplots,
  value = "D", fun = sd, per_ha = TRUE,
  AGB_simu = error_prop$AGB_simu, 
  ref_raster = nouraguesRaster # by default, the associated function is the mean function
  )

## ----customize_plot-----------------------------------------------------------
subplot_metric <- subplot_summary(subplots = subplots,
                                  value = "AGB") 

custom_plot <- subplot_metric$plot_design
# Change the title and legend:
custom_plot + 
  labs(title = "Nouragues plot" , fill="Sum of AGB per hectare")
# Display trees with diameter as size and transparency (and a smaller legend on the right): 
custom_plot + 
  geom_point(data=check_plot_trust_GPS$tree_data, mapping = aes(x = x_proj, y = y_proj, size = D, alpha= D), shape=1,) +
  labs(fill = "Sum of AGB per hectare") +
  guides(alpha = guide_legend(title = "Diameter (cm)"),
         size = guide_legend(title = "Diameter (cm)")) + 
  theme(legend.position = "right", legend.key.size = unit(0.5, 'cm'))

