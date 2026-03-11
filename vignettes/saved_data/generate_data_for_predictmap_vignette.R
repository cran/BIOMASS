###draft vignette

# load data
data("NouraguesTrees")
data("NouraguesCoords")
nouraguesRaster <- terra::rast(system.file("extdata", "NouraguesRaster.tif",package = "BIOMASS", mustWork = TRUE))

#1) compute AGBD
# # Height
# data("NouraguesHD")
# brm_model <- modelHD(
#   D = NouraguesHD$D, H = NouraguesHD$H,
#   method = "log2",
#   bayesian = TRUE, useCache = TRUE)
# 
# # WD
# Taxo <- readRDS(file = "vignettes/saved_data/Taxo_vignette.rds")
# NouraguesTrees$GenusCorrected <- Taxo$genusAccepted
# NouraguesTrees$SpeciesCorrected <- Taxo$speciesAccepted
# NouraguesTrees$family <- Taxo$familyAccepted
# wood_densities <- getWoodDensity(
#   genus = NouraguesTrees$GenusCorrected,
#   species = NouraguesTrees$SpeciesCorrected,
#   family = NouraguesTrees$family,
#   stand = NouraguesTrees$Plot # for unidentified or non-documented trees in the reference database
# )
# NouraguesTrees$WD <- wood_densities$meanWD
# 
# error_prop_4plots <- AGBmonteCarlo(
#   D = NouraguesTrees$D, WD = NouraguesTrees$WD, 
#   HDmodel = brm_model, 
#   Dpropag = "chave2004",
#   errWD = wood_densities$sdWD)
# # keep only 200 iterations per tree
# error_prop_4plots$AGB_simu <- error_prop_4plots$AGB_simu[,1:50]
# saveRDS(error_prop_4plots, file = "vignettes/saved_data/error_prop4plots.rds")

error_prop4plots <- readRDS(file = "vignettes/saved_data/error_prop4plots.rds")

#2) spatialize AGBD

# divide plots into subplots
multiple_subplots <- divide_plot(
  grid_size = 25, 
  corner_data = NouraguesCoords,
  rel_coord = c("Xfield","Yfield"), proj_coord = c("Xutm","Yutm"), corner_plot_ID = "Plot",
  tree_data = NouraguesTrees, tree_coords = c("Xfield","Yfield"), tree_plot_ID = "Plot"
)

# multiple_checks <- check_plot_coord(
#   corner_data = NouraguesCoords, # NouraguesCoords contains 4 plots
#   proj_coord = c("Xutm", "Yutm"), rel_coord = c("Xfield", "Yfield"),
#   trust_GPS_corners = TRUE,
#   plot_ID = "Plot",
#   tree_data = NouraguesTrees, tree_coords = c("Xfield","Yfield"),
#   prop_tree = "D", tree_plot_ID = "Plot",
#   ref_raster = nouraguesRaster)

# compute AGBD estimates and their uncertainty per subplot
subplot_AGBD <- subplot_summary(
  subplots = multiple_subplots,
  AGB_simu = error_prop_4plots$AGB_simu, draw = F
)

#3) spatialize LiDAR metric
# quick plot to visualise plot corners in the landscape
terra::plot(nouraguesRaster)
points(NouraguesCoords$Xutm, NouraguesCoords$Yutm, col ="red", pch = 20)

# get CHM median values for each suplot
raster_summary <- subplot_summary(
  subplots = multiple_subplots,
  ref_raster = nouraguesRaster, raster_fun = median, na.rm = T)
chm_subplot <- raster_summary$tree_summary %>%
  rename(raster_metric = z2012_median)

#4) calibrate model
# gather data for agbd-chm model calibration
agbd_subplot <- subplot_AGBD$long_AGB_simu
data_for_calibration <- agbd_subplot %>%
  left_join(chm_subplot, by = "subplot_ID") %>%
  arrange(subplot_ID)

## calibrate model
# model_cal <- calibrate_model(long_AGB_simu = data_for_calibration, nb_rep = 50, useCache = T,
#                             plot_model = TRUE, chains = 4, thin = 20, iter = 2500,
#                             warmup = 500, cores = 4)
# saveRDS(model_cal, file = "vignettes/saved_data/nouragues_model_cal.rds")

model_cal <- readRDS(file = "vignettes/saved_data/nouragues_model_cal.rds")

# check
summary(model_cal)

#5) predict AGBD map
map_agbd <- predict_map(fit_brms = model_cal,
                         pred_raster = nouraguesRaster,
                         grid_size = 25,
                         raster_fun = median,
                         n_cores = 4,
                         n_post_draws = 100,
                         alignment_raster = NULL,
                         plot_maps = T)
#saveRDS(map_agbd, file = "vignettes/saved_data/nouragues_map_pred.rds")
