# Angetter et al. (2011) -----------------
# https://doi.org/10.1111/j.1095-8312.2011.01780.x

# 1 Preparation ---------------------

# 1.1 Data sources ==================
# Source of occurence data:
# Gbif
# Reference code:
# GBIF.org (29 March 2023) GBIF Occurrence Download  https://doi.org/10.15468/dl.jh2rpq

# add vertnet reference

# Source of climate data:
# Worldclim, 30 arcsec, bio: 1, 2, 5, 8, 9, 10, 12, 13, 14, 18.
# link = https://www.worldclim.org/data/worldclim21.html

# original

# 1.2 Package library =============
# Packages always used:
library(tidyverse)
library(terra)
library(INLA)
library(pROC)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(ade4)
library(ecospat)

library(hypervolume)

# library(tidylog)
# library(sp)
# # inla.upgrade() # for the stable version
# library(rgdal)
# library(vegan)
# library(ggfortify)
# library(mapproj)
# library(ade4)
# library(ecospat)
# # library(dismo)

# Packages used, but not loaded:
# library(raster)
# library(spatialEco)
# library(isocat)

# Packages used for data prep:
library(readxl)
library(here)

## 1.3 Function files ========
source(paste0(here(), "/Joe_functions.R"))
source(paste0(here(), "/Isak_functions.R"))

#Note new occurrence data ----------
# range(-104.76345 9.42033,-68.8774 9.42033,-68.8774 43.83465,-104.76345 43.83465,-104.76345 9.42033)
# time (1911 - 2011)



## 1.4 Import occurrence ========
anolis_occurence <- read_xlsx(path = paste0(here(), "/data/anolis_sagrei_v2.xlsx"))


anolis_occurence_clean <- anolis_occurence |>
  mutate(Latitude = as.double(decimalLatitude), Longitude = as.double(decimalLongitude)) |>
  select(Latitude, Longitude, occurrenceStatus) |>
  replace(
    is.character(anolis_occurence$occurrenceStatus) == TRUE,
    anolis_occurence$occurrenceStatus <- 1
  )

Vertnet <- read.csv(file = "data/VertNet_Reptilia_Oct2015/VertNet_Reptilia_Oct2015.csv")

anolis_occurence2 <- Vertnet |> 
  filter(scientificname == "Anolis sagrei", year < 2011)

anolis_occurence_clean2 <- anolis_occurence2 |>
  mutate(Latitude = as.double(decimallatitude), Longitude = as.double(decimallongitude)) |>
  select(Latitude, Longitude, occurrenceStatus) |>
  replace(is.character(anolis_occurence2$occurrenceStatus) == TRUE, anolis_occurence2$occurrenceStatus <- 1)

anolis_occurence_final <- full_join(anolis_occurence_clean, anolis_occurence_clean2)

## 1.5 Import climate data ========
bio1 <- rast(paste0(here(), "/data/wc2.1_30s_bio/wc2.1_30s_bio_1.tif"))
bio2 <- rast(paste0(here(), "/data/wc2.1_30s_bio/wc2.1_30s_bio_2.tif"))
bio5 <- rast(paste0(here(), "/data/wc2.1_30s_bio/wc2.1_30s_bio_5.tif"))
bio8 <- rast(paste0(here(), "/data/wc2.1_30s_bio/wc2.1_30s_bio_8.tif"))
bio9 <- rast(paste0(here(), "/data/wc2.1_30s_bio/wc2.1_30s_bio_9.tif"))
bio10 <- rast(paste0(here(), "/data/wc2.1_30s_bio/wc2.1_30s_bio_10.tif"))
bio12 <- rast(paste0(here(), "/data/wc2.1_30s_bio/wc2.1_30s_bio_12.tif"))
bio13 <- rast(paste0(here(), "/data/wc2.1_30s_bio/wc2.1_30s_bio_13.tif"))
bio14 <- rast(paste0(here(), "/data/wc2.1_30s_bio/wc2.1_30s_bio_14.tif"))
bio18 <- rast(paste0(here(), "/data/wc2.1_30s_bio/wc2.1_30s_bio_18.tif"))


climfull_x <- terra::`add<-`(bio1, c(bio2, bio5, bio8, bio9, bio10, bio12, bio13, bio14, bio18))

# Make a list of names for each climate variable (make sure to have the correct order)
clim_names <- c("Annual Mean Temperature (C)",
                "Mean Diurnal Range (Mean of monthly (max temp - min temp, C))",
                "Max Temperature of Warmest Month (C)", "Mean Temperature of Wettest Quarter (C)",
                "Mean Temperature of Driest Quarter (C)", "Mean Temperature of Warmest Quarter (C)",
                "Annual Precipitation (mm)", "Precipitation of Wettest Month (mm)",
                "Precipitation of Driest Month (mm)", "Precipitation of Warmest Quarter (mm)")

variable_codes <- c("wc2.1_30s_bio_1", "wc2.1_30s_bio_02",
                    "wc2.1_30s_bio_5", "wc2.1_30s_bio_8",
                    "wc2.1_30s_bio_9", "wc2.1_30s_bio_10",
                    "wc2.1_30s_bio_12", "wc2.1_30s_bio_13",
                    "wc2.1_30s_bio_14", "wc2.1_30s_bio_18")

## 1.6 Define expanded prediction range ==============
# range(-104.76345 ,-68.8774, 9.42033, 43.83465)
model_range_extended <- c(-110, -67.5, 8, 44)

climate_extended_prediction_range <- define_range_and_aggregation(climfull, model_range_extended,
  aggregation_factor = 1
) |>
  raster::stack() |>
  as("SpatialGridDataFrame")


## 1.7 Native population ==============

model_range_native <- c(-86, -73.7, 19.4, 27.5)

# clean occurrence data
anolis_native <- anolis_occurence_final |>
  filter(Longitude > model_range_native[1], Longitude < model_range_native[2],
    Latitude > model_range_native[3], Latitude < model_range_native[4]
  )

native_exclusion_range1 <- c(-84, -79.5, 24, 27.5)
native_exclusion_range2 <- c(-81, -79, 19.4, 20)

anolis_native_exclusion1 <- anolis_occurence_final |>
  filter(Longitude > native_exclusion_range1[1], Longitude < native_exclusion_range1[2],
         Latitude > native_exclusion_range1[3], Latitude < native_exclusion_range1[4])

anolis_native_exclusion2 <- anolis_occurence_final |>
  filter(Longitude > native_exclusion_range2[1], Longitude < native_exclusion_range2[2],
         Latitude > native_exclusion_range2[3], Latitude < native_exclusion_range2[4])

anolis_native_merge <- anti_join(anolis_native, anolis_native_exclusion1)

anolis_native_final <- anti_join(anolis_native_merge, anolis_native_exclusion2)

# Prepare data
climate_native <- define_range_and_aggregation(climfull_x, model_range_native,
  aggregation_factor = 1
)

climate_native_exclusion1 <- define_range_and_aggregation(climfull_x, native_exclusion_range1)

climate_native_exclusion_null1 <- setValues(climate_native_exclusion1, values = NA)

climate_native_merge <- terra::merge(climate_native_exclusion_null1, climate_native)

climate_native_exclusion2 <- define_range_and_aggregation(climfull_x, native_exclusion_range2)

climate_native_exclusion_null2 <- setValues(climate_native_exclusion2, values = NA)

climate_native_final <- terra::merge(climate_native_exclusion_null2, climate_native_merge)

#4
climate_grid_native <- make_climate_grid(climate_native)

occurence_grid_native <- make_occurence_grid(anolis_native_final, climate_native)


## 1.8 Invasive population ==============

model_range_invasive <- c(-99.6, -75.6, 24.2, 37.9)

# clean occurrence data
anolis_invasive <- anolis_occurence_final |>
  filter(
    Longitude > model_range_invasive[1], Longitude < model_range_invasive[2],
    Latitude > model_range_invasive[3], Latitude < model_range_invasive[4]
  )

anolis_invasive_exclusion_range <- c(-79.5, -72, 24.2, 28)

anolis_native_exclusion <- anolis_occurence_final |>
  filter(Longitude > anolis_invasive_exclusion_range[1], Longitude < anolis_invasive_exclusion_range[2],
         Latitude > anolis_invasive_exclusion_range[3], Latitude < anolis_invasive_exclusion_range[4])


anolis_invasive_final <- anti_join(anolis_invasive, anolis_native_exclusion)

# (max(anolis_invasive_final$Longitude) - min(anolis_invasive_final$Longitude))*0.05
# (max(anolis_invasive_final$Latitude) - min(anolis_invasive_final$Latitude))*0.05
# summary(anolis_invasive_final)
# 
# coastline <- ne_countries(continent = c("south america", "north america"), returnclass = "sf", scale = 50)
# 
# ggplot() +
#   geom_sf(data = coastline, fill = "white") +
#   geom_point(data = anolis_invasive_final, aes(Longitude, Latitude, color = "native"))+
#   xlim(-99.6, -75.6)+
#   ylim(24.2, 37.9)


# Prepare data
climate_invasive <- define_range_and_aggregation(climfull_x, model_range_invasive,
  aggregation_factor = 1
)
#9
climate_invasive_exclusion <- define_range_and_aggregation(climfull_x, anolis_invasive_exclusion_range)

climate_invasive_exclusion_null <- setValues(climate_invasive_exclusion, values = NA)

climate_invasive_new <- terra::merge(climate_invasive_exclusion_null, climate_invasive)

# climate_invasive_new <- raster::aggregate(climate_invasive_new, fact = 9)

climate_grid_invasive <- make_climate_grid(climate_invasive_new)

occurence_grid_invasive <- make_occurence_grid(anolis_invasive, climate_invasive_new)


## 1.9 Plot occurrences ==============

# coastline <- ne_countries(continent = "north america", returnclass = "sf", scale = 50) |>
#   st_crop(c(xmin = -105, xmax = -70, ymin = 10, ymax = 38))
# 
# ggplot() +
#   geom_sf(data = coastline, fill = "white") +
#   geom_point(data = anolis_invasive, aes(Longitude, Latitude, color = "invasive")) +
#   geom_point(data = anolis_native, aes(Longitude, Latitude, color = "native")) +
#   scale_x_continuous(expand = c(0, 0), limits = c(-102, -70)) +
#   scale_y_continuous(expand = c(0, 0), limits = c(13, 35)) +
#   labs(title = "Species occurences (1990-2020)", x = NULL, y = NULL, color = "Population") +
#   theme(panel.background = element_rect(fill = "lightblue", colour = "lightblue"))

# 2 Model --------------

## 2.1 Native ============

mod_native <- sdmINLASpatRandEff(occurence_grid_native, climate_grid_native,
  meshParameters = list(max.n.strict = -1),
  createGeoTIFF = FALSE, outFolder = getwd(),
  myName = "Anolis_sagrei_Native"
)

## 2.1 Invasive ============
mod_invasive <- sdmINLASpatRandEff(occurence_grid_invasive, climate_grid_invasive,
  meshParameters = list(max.n.strict = -1),
  createGeoTIFF = FALSE,
  outFolder = getwd(),
  myName = "Anolis_sagrei_Invasive"
)

# 3 Model Results -------------
# 
# ## 3.1 Native Results =========
# # Retrive model results
# obj_native <- readRDS(file = "Anolis_sagrei_Native.rds")
# 
# ### 3.1.1 Mean occurrence probability map #########
# native_occurence_estimate <- as.data.frame(obj_native$spatialPredictions[1])
# 
# ggplot(native_occurence_estimate) +
#   geom_tile(aes(s1, s2, fill = meanEst)) +
#   coord_sf() +
#   scale_fill_viridis_c(
#     limits = c(0, 0.35), option = "magma",
#     begin = 0.12, end = 1, name = "Occurence probability"
#   ) +
#   scale_x_continuous(expand = c(0, 0), limits = c(model_range_native[1], model_range_native[2])) +
#   scale_y_continuous(expand = c(0, 0), limits = c(model_range_native[3], model_range_native[4])) +
#   labs(title = "Native occurrence prediction", x = NULL, y = NULL) +
#   theme_bw()
# 
# 
# ### 3.1.2 Spatial effect map #########
# native_spatial_effect <- as.data.frame(obj_native$spatialPredictions[7])
# 
# ggplot(native_spatial_effect, aes(s1, s2)) +
#   geom_tile(aes(fill = meanSpatialPred)) +
#   scale_fill_viridis_c(
#     option = "magma",
#     end = 0.9, name = "Spatial effect"
#   ) +
#   scale_x_continuous(expand = c(0, 0), limits = c(model_range_native[1], model_range_native[2])) +
#   scale_y_continuous(expand = c(0, 0), limits = c(model_range_native[3], model_range_native[4])) +
#   coord_sf() +
#   labs(title = "Spatial effect in native range", x = NULL, y = NULL) +
#   theme_bw()
# 
# ### 3.1.3 Predicted occurrence probability in an expanded range #########
# predict_native <- predictSD(obj_native, climate_extended_prediction_range,
#   origClimateData = climate_grid_native
# )
# native_predicted_estimate <- as.data.frame(predict_native[1])
# 
# ggplot(native_predicted_estimate, aes(s1, s2)) +
#   geom_tile(aes(fill = OccurrenceProb)) +
#   coord_sf() +
#   scale_fill_viridis_c(
#     limits = c(0, 1), option = "magma",
#     begin = 0.12, end = 1, name = "Suitability"
#   ) +
#   scale_x_continuous(expand = c(0, 0), limits = c(model_range_extended[1], model_range_extended[2])) +
#   scale_y_continuous(expand = c(0, 0), limits = c(model_range_extended[3], model_range_extended[4])) +
#   labs(title = "Native niche prediction", x = NULL, y = NULL) +
#   theme_bw()
# 
# ### 3.1.4 Individual climate response curves ########
# preds_native <- obj_native$responsePredictions
# 
# pars_native <- getPars(obj_native, climate_grid_native)
# 
# # Plot individual response curves
# map2(preds_native, clim_names, plot_response_curve)
# 
# # Calculate & then plot standardized individual response curves
# integral_native_bio5 <- calculate_integral_response_curve(
#   climate_grid_native, "wc2.1_30s_bio_5",
#   preds_native$wc2.1_30s_bio_5$covarVal,
#   pars_native
# )
# 
# #NOTE: Need help with function ------------
# integral_native <- map2(.x = variable_names, .y = preds_native,
#                         calculate_integral_response_curve,
#                         data = climate_grid_native, pars = pars_native,
#                         x_range = .x$covarVal)
# 
# plot_integral_response_curve(integral_native_bio5, "Name of the climate variable")
# 
# 
# ### 3.1.5 AUC ########
# observations_native <- occurence_grid_native$occurrence
# predictions_native <- obj_native$spatialPredictions$meanEst
# roc_object_native <- roc(observations_native, predictions_native)
# anolis_AUC_native <- auc(roc_object_native)
# anolis_AUC_native
# 
# 
# ## 3.2 Invasive Results =========
# # Retrive model results
# obj_invasive <- readRDS(file = "Anolis_sagrei_invasive.rds")
# 
# ### 3.2.1 Mean occurrence probability map #########
# invasive_occurence_estimate <- as.data.frame(obj_invasive$spatialPredictions[1])
# 
# ggplot(invasive_occurence_estimate) +
#   geom_tile(aes(s1, s2, fill = meanEst)) +
#   coord_sf() +
#   scale_fill_viridis_c(
#     limits = c(0, 1), option = "magma",
#     begin = 0.12, end = 1, name = "Occurence probability"
#   ) +
#   scale_x_continuous(expand = c(0, 0), limits = c(model_range_invasive[1], model_range_invasive[2])) +
#   scale_y_continuous(expand = c(0, 0), limits = c(model_range_invasive[3], model_range_invasive[4])) +
#   labs(title = "invasive occurrence prediction", x = NULL, y = NULL) +
#   theme_bw()
# 
# 
# ### 3.2.2 Spatial effect map #########
# invasive_spatial_effect <- as.data.frame(obj_invasive$spatialPredictions[7])
# 
# ggplot(invasive_spatial_effect, aes(s1, s2)) +
#   geom_tile(aes(fill = meanSpatialPred)) +
#   scale_fill_viridis_c(
#     option = "magma",
#     end = 0.9, name = "Spatial effect"
#   ) +
#   scale_x_continuous(expand = c(0, 0), limits = c(model_range_invasive[1], model_range_invasive[2])) +
#   scale_y_continuous(expand = c(0, 0), limits = c(model_range_invasive[3], model_range_invasive[4])) +
#   coord_sf() +
#   labs(title = "Spatial effect in invasive range", x = NULL, y = NULL) +
#   theme_bw()
# 
# ### 3.2.3 Predicted occurrence probability in an expanded range #########
# predict_invasive <- predictSD(obj_invasive, climate_extended_prediction_range,
#   origClimateData = climate_grid_invasive
# )
# invasive_predicted_estimate <- as.data.frame(predict_invasive[1])
# 
# ggplot(invasive_predicted_estimate, aes(s1, s2)) +
#   geom_tile(aes(fill = OccurrenceProb)) +
#   coord_sf() +
#   scale_fill_viridis_c(
#     limits = c(0, 1), option = "magma",
#     begin = 0.12, end = 1, name = "Suitability"
#   ) +
#   scale_x_continuous(expand = c(0, 0), limits = c(model_range_extended[1], model_range_extended[2])) +
#   scale_y_continuous(expand = c(0, 0), limits = c(model_range_extended[3], model_range_extended[4])) +
#   labs(title = "invasive niche prediction", x = NULL, y = NULL) +
#   theme_bw()
# 
# ### 3.2.4 Individual climate response curves ########
# preds_invasive <- obj_invasive$responsePredictions
# 
# pars_invasive <- getPars(obj_invasive, climate_grid_invasive)
# 
# # Plot individual response curves
# map2(preds_invasive, clim_names, plot_response_curve)
# 
# 
# # Calculate & then plot standardized individual response curves
# integral_invasive_bio5 <- calculate_integral_response_curve(
#   climate_grid_invasive, "wc2.1_30s_bio_5",
#   preds_invasive$wc2.1_30s_bio_5$covarVal,
#   pars_invasive
# )
# 
# plot_integral_response_curve(integral_invasive_bio5, "Name of the climate variable")
# 
# ### 3.2.5 AUC ########
# observations_invasive <- occurence_grid_invasive$occurrence
# predictions_invasive <- obj_invasive$spatialPredictions$meanEst
# roc_object_invasive <- roc(observations_invasive, predictions_invasive)
# anolis_AUC_invasive <- auc(roc_object_invasive)
# anolis_AUC_invasive
# 
# 
# # 4 Niche comparison ---------
# 
# 
# ## 4.1 Univariate ==========
# 
# 
# 
# compare_response_curves(
#   preds_native$wc2.1_30s_bio_1,
#   preds_invasive$wc2.1_30s_bio_1,
#   "a"
# )
# 
# pmap(c(preds_native, preds_invasive, clim_names), compare_response_curves)
# 
# 
# # PRIORITY ------------
# purrr::map2
# preds_native$wc2.1_30s_bio_1[1]
# 
# y <- map2(preds_native, )
# 
# map(preds_native, plot_response_curve, "a")
# map2(preds_native, clim_names, plot_response_curve)
# 
# 
# # Pars for native and invasive
# pars_native <- getPars(obj_native, climate_grid_native)
# pars_invasive <- getPars(obj_invasive, climate_grid_invasive)
# 
# # Standardised curves
# 
# compare_response_curves(integral_native_bio5, integral_invasive_bio5)
# 
# 
# # ## 4.2 Multivariate ===========
# # 
# # climate_background <- as.data.frame(climate_extended_prediction_range) |>
# #   select(-s1, -s2) |>
# #   rename(
# #     "Annual Mean Temperature" = wc2.1_30s_bio_1,
# #     "Mean Diurnal Range" = wc2.1_30s_bio_2,
# #     "Max Temperature of Warmest Month" = wc2.1_30s_bio_5,
# #     "Mean Temperature of Wettest Quarter" = wc2.1_30s_bio_8,
# #     "Mean Temperature of Driest Quarter" = wc2.1_30s_bio_9,
# #     "Mean Temperature of Warmest Quarter" = wc2.1_30s_bio_10,
# #     "Annual Precipitation" = wc2.1_30s_bio_12,
# #     "Precipitation of Wettest Month" = wc2.1_30s_bio_13,
# #     "Precipitation of Driest Month" = wc2.1_30s_bio_14,
# #     "Precipitation of Warmest Quarter" = wc2.1_30s_bio_18
# #   )
# # 
# # occurence_native_data_frame <- as.data.frame(occurence_grid_native)
# # 
# # # Data frame showing each climate variable at each point in the native area
# # native_data_frame <- as.data.frame(climate_grid_native) |>
# #   left_join(occurence_native_data_frame) |>
# #   select(-s1, -s2)
# # 
# # # Add predicted occurrence probability from the native model to extended climate area data frame
# # native_prediction_data_frame <- as.data.frame(climate_extended_prediction_range) |>
# #   left_join(native_predicted_estimate) |>
# #   select(-s1, -s2)
# # 
# # occurence_invasive_data_frame <- as.data.frame(occurence_grid_invasive)
# # 
# # # Data frame showing each climate variable at each point in the invasive area
# # invasive_data_frame <- as.data.frame(climate_grid_invasive) |>
# #   left_join(occurence_invasive_data_frame) |>
# #   select(-s1, -s2)
# # 
# # # Add predicted occurrence probability from the native model to extended climate area data frame
# # invasive_prediction_data_frame <- as.data.frame(climate_extended_prediction_range) |>
# #   left_join(invasive_predicted_estimate) |>
# #   select(-s1, -s2)
# # 
# # # Make a pca of the climate in the extended prediction area
# # pca.env <- dudi.pca(climate_background, scannf = F, nf = 2)
# # 
# # # Plots the contribution of each
# # ecospat.plot.contrib(contrib = pca.env$co, eigen = pca.env$eig)
# # 
# # ## Convert
# # # PCA scores for the whole study area
# # scores.globclim <- pca.env$li
# # # PCA scores for the species native distribution
# # scores.sp.nat <- suprow(pca.env, native_data_frame[which(native_data_frame[, 11] == 1), 1:10])$li
# # # PCA scores for the species invasive distribution
# # scores.sp.inv <- suprow(pca.env, invasive_data_frame[which(invasive_data_frame[, 11] == 1), 1:10])$li
# # # PCA scores for the whole native study area
# # scores.clim.nat <- suprow(pca.env, native_data_frame[, 1:10])$li
# # # PCA scores for the whole invaded study area
# # scores.clim.inv <- suprow(pca.env, invasive_data_frame[, 1:10])$li
# # # PCA scores for the occurrence prediction based on the native model
# # scores.pred.nat <- suprow(pca.env, native_prediction_data_frame[, 1:10])$li
# # # PCA scores for the occurrence prediction based on the invasive model
# # scores.pred.inv <- suprow(pca.env, invasive_prediction_data_frame[, 1:10])$li
# # 
# # # Assigns the occurrence probability from the native model to coordinates in the pca
# # clim_space_native <- bind_cols(scores.pred.nat, native_predicted_estimate) |>
# #   select(-s1, -s2)
# # # Assigns the occurrence probability from the invasive model to coordinates in the pca
# # clim_space_invasive <- bind_cols(scores.pred.inv, invasive_predicted_estimate) |>
# #   select(-s1, -s2)
# # 
# # 
# # # Gridding the native niche
# # grid.clim.nat <- ecospat.grid.clim.dyn(
# #   glob = scores.globclim,
# #   glob1 = scores.clim.nat,
# #   sp = scores.sp.nat, R = 100,
# #   th.sp = 0
# # )
# # 
# # # Gridding the invasive niche
# # grid.clim.inv <- ecospat.grid.clim.dyn(
# #   glob = scores.globclim,
# #   glob1 = scores.clim.inv,
# #   sp = scores.sp.inv, R = 100,
# #   th.sp = 0
# # )
# # 
# # 
# # # Plots the available climate and occurrences in the native and invasive ranges in a pca
# # ecospat.plot.niche.dyn(grid.clim.nat, grid.clim.inv,
# #   quant = 0.25, interest = 1,
# #   title = "Niche Overlap", name.axis1 = "PC1 = 40.85%",
# #   name.axis2 = "PC2 = 25.64%"
# # )
# # ecospat.shift.centroids(scores.sp.nat, scores.sp.inv, scores.clim.nat, scores.clim.inv)
# # 
# # # Save the coordinates of occurrence probabilities in the pca at several thresholds
# # # No threshold saves all points from the extended prediction area
# # coords_100 <- clim_space_native |>
# #   select(-OccurrenceProb) |>
# #   as.matrix()
# # 
# # coords_native_75 <- clim_space_native |>
# #   filter(OccurrenceProb > 0.25) |>
# #   select(-OccurrenceProb) |>
# #   as.matrix()
# # 
# # coords_native_50 <- clim_space_native |>
# #   filter(OccurrenceProb > 0.5) |>
# #   select(-OccurrenceProb) |>
# #   as.matrix()
# # 
# # coords_native_25 <- clim_space_native |>
# #   filter(OccurrenceProb > 0.75) |>
# #   select(-OccurrenceProb) |>
# #   as.matrix()
# # 
# # 
# # coords_invasive_75 <- clim_space_invasive |>
# #   filter(OccurrenceProb > 0.25) |>
# #   select(-OccurrenceProb) |>
# #   as.matrix()
# # 
# # coords_invasive_50 <- clim_space_invasive |>
# #   filter(OccurrenceProb > 0.5) |>
# #   select(-OccurrenceProb) |>
# #   as.matrix()
# # 
# # coords_invasive_25 <- clim_space_invasive |>
# #   filter(OccurrenceProb > 0.75) |>
# #   select(-OccurrenceProb) |>
# #   as.matrix()
# # 
# # # Make hulls
# # hull_100 <- inla.nonconvex.hull(coords_100, convex = -0.02, concave = -0.15, resolution = 300)
# # hull_native_75 <- inla.nonconvex.hull(coords_native_75, convex = -0.02, concave = -0.15, resolution = 300)
# # hull_native_50 <- inla.nonconvex.hull(coords_native_50, convex = -0.02, concave = -0.15, resolution = 300)
# # hull_native_25 <- inla.nonconvex.hull(coords_native_25, convex = -0.02, concave = -0.15, resolution = 300)
# # 
# # hull_invasive_75 <- inla.nonconvex.hull(coords_invasive_75, convex = -0.02, concave = -0.15, resolution = 300)
# # hull_invasive_50 <- inla.nonconvex.hull(coords_invasive_50, convex = -0.02, concave = -0.15, resolution = 300)
# # hull_invasive_25 <- inla.nonconvex.hull(coords_invasive_25, convex = -0.02, concave = -0.15, resolution = 300)
# # 
# # ecospat.plot.niche(grid.clim.nat)
# # lines(hull_100$loc, col = "red")
# # lines(hull_native_75$loc, col = "red")
# # lines(hull_native_50$loc, col = "red")
# # lines(hull_native_25$loc, col = "red")
# # # text() to ad values to the lines
# # # polygon() for fill
# # 
# # 
# # ecospat.plot.niche(grid.clim.inv)
# # lines(hull_invasive_75$loc, col = "blue")
# # 
# # plot(hull_100$loc,
# #   type = "lines",
# #   xlab = "PC1 = 40.85%", ylab = "PC2 = 25.64%", main = "Predicted niches in climate space"
# # )
# # polygon(hull_invasive_75$loc, col = adjustcolor("blue", alpha.f = 0.1))
# # polygon(hull_invasive_50$loc, col = adjustcolor("blue", alpha.f = 0.3))
# # polygon(hull_invasive_25$loc, col = adjustcolor("blue", alpha.f = 0.5))
# # polygon(hull_native_75$loc, col = adjustcolor("red", alpha.f = 0.1))
# # polygon(hull_native_50$loc, col = adjustcolor("red", alpha.f = 0.3))
# # polygon(hull_native_25$loc, col = adjustcolor("red", alpha.f = 0.5))
# # points(scores.sp.inv, cex = 0.7, pch = 2, col = "red")
# # points(scores.sp.nat, cex = 0.7, pch = 1, col = "blue")
# # # look at legend()
# # 
# # 
# # 
# # ## 4.3 Overlap index =============
# # y <- as.data.frame(c(climate_extended_prediction_range$wc2.1_30s_bio_1,
# #                      climate_extended_prediction_range$wc2.1_30s_bio_2,
# #                      climate_extended_prediction_range$wc2.1_30s_bio_5,
# #                      predict_native))
# # 
# # 
# # z <- as.data.frame(c(climate_extended_prediction_range$wc2.1_30s_bio_1,
# #                      climate_extended_prediction_range$wc2.1_30s_bio_2,
# #                      climate_extended_prediction_range$wc2.1_30s_bio_3))
# # 
# # 
# # z <- as.data.frame(climate_extended_prediction_range)
# # 
# # z <- z |> 
# #   select(-s1, -s2)
# # #Make toy data ------------------
# # #standardised mean/sd
# # test_df <- as.data.frame(climate_grid_native) |> 
# #   select(wc2.1_30s_bio_1, wc2.1_30s_bio_2, wc2.1_30s_bio_5)
# # test_hyper <- hypervolume(test_df)
# # 
# # #
# # comp_predict_native <- predictSD(obj_native, climate_grid_native,
# #                                  origClimateData = climate_grid_native
# # )
# # comp_native_predicted_estimate <- as.data.frame(comp_predict_native[1])
# # #
# # comp_df_native <- as.data.frame(climate_grid_native) |> 
# #   left_join(comp_native_predicted_estimate) |> 
# #   select(wc2.1_30s_bio_1, wc2.1_30s_bio_2, wc2.1_30s_bio_5, OccurrenceProb)
# # 
# # comp_df_test <- as.data.frame(climate_grid_native) |> 
# #   select(wc2.1_30s_bio_1, wc2.1_30s_bio_2, wc2.1_30s_bio_5) |> 
# #   mutate(OccurrenceProb = rnorm(5433, 0.5, 0.1))
# # 
# # 
# # comp_hyper_native <- hypervolume(comp_df_native)
# # comp_niche_native <- hypervolume_to_data_frame(comp_hyper_native)
# # 
# # comp_hyper_test <- hypervolume(comp_df_test)
# # comp_niche_test <- hypervolume_to_data_frame(comp_hyper_test)
# # 
# # #
# # comp_predict_native <- predictSD(obj_native, climate_grid_native,
# # origClimateData = climate_grid_native
# # )
# # comp_native_predicted_estimate <- as.data.frame(comp_predict_native[1])
# # #
# # 
# # x <- hypervolume(y)
# # w <- hypervolume(z)
# # 
# # A <- raster::raster(predict_invasive)
# # B <- raster::raster(predict_native)
# # 
# # isocat::schoenersD(A / raster::cellStats(A, stat = "sum"), B / raster::cellStats(B, stat = "sum"))
# # 
# # a <- rast(predict_invasive)
# # b <- rast(predict_native)
# # 
# # spatialEco::overlap(a, b)
