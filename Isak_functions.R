
#' @title Crop the range of a climate SpatRaster or raster map and aggregate if necessary
#' 
#' @description 
#' This function crops a SpatRaster or raster object to a range of coordinates and aggregates the cells if specified
#' 
#' @param climate_data A SpatRaster or raster object
#' @param range_border A list of four values defining the border of the cropped object.
#'List must have values in order: West, East, South, North.
#' @param aggregation_factor A value determining what factor to aggregate the cells by.
#' Default is set to 1 causing no aggregation
#' WARNING! Aggregation is not intended for final analysis. The model code requires a lot of memory and will crash if not provided with enough. Aggregation is intended for trial runs while working while computers with insufficient memory and to reduce time. It is unlikely that a consumer grade computer will be able to complete an unaggregated analysis given the the size and resolution of study areas form most of invasive niche shift studies.
define_range_and_aggregation <- function(climate_data, range_border, aggregation_factor = 1) {
  clim <- climate_data |> 
    crop(raster::extent(range_border))
  
  clim <- raster::aggregate(clim, fact = aggregation_factor)
  
  return(clim)
}



#' @title Turn a SpatRaster into a SpatialGridDataFrame
#' 
#' @description 
#' This function takes a SpatRaster and converts it into a raster stack that in turn is made into a SpatialGridDataFrame
#' 
#' @param clim A SpatRaster
make_climate_grid <- function(clim) {
  clim_raster <- clim |> 
    raster::stack() |> 
    as("SpatialGridDataFrame")
  
  return(clim_raster)
}



#' @title Crop the range of a climate SpatRaster or raster map and aggregate if necessary
#' 
#' @description 
#' This function crops a SpatRaster or raster object to a range of coordinates and aggregates the cells if specified
#' 
#' @param occurence_data A data frame containing 3 columns:
#'  Longitude, Latitude, occurrenceStatus (which should have values of 1 for present)
#' @param clim A SpatRaster of the climate data containing all coordinates of the occurence data
#' 
#'  @return A raster of occurences with cell values of 1 for present and 0 for absence. Resolution is identical to the resolution of @param clim
make_occurence_grid <- function(occurence_data, clim) {
  grid <- occurence_data|> 
    vect(geom = c("Longitude", "Latitude")) |> 
    rasterize(clim[[1]], background = 0, field = "occurrenceStatus")
  names(grid) <- "occurrence"
  
  grid_raster <- grid |> 
    raster::raster() |> 
    as("SpatialGridDataFrame")
  names(grid_raster) <- "occurrence"
  
  return(grid_raster)
}




#' @title Plot a climate response curve with confidence interval
#' 
#' @param data A data frame containing 4 columns:
#'  covarVal (the climate gradient), meanEst (mean response),
#'  lowerEst (lower confidence) & upperEst (upper confidence)
#' @param clim_variable_name A string containing the name of the climate variable
#'  (will be displayed on the x-axis)
#' @param y_limit A range that defines the y-axis (used to investigate low effect curves)
plot_response_curve <- function(data, clim_variable_name = "Name of the climate variable", y_limit = c(0,1)) {

    ggplot(data, aes(covarVal, meanEst))+
      geom_line(linewidth = 0.6)+
    geom_line(aes(covarVal, lowerEst), linetype = "dashed")+
    geom_line(aes(covarVal, upperEst), linetype = "dashed")+
    ylim(y_limit)+
    labs(y = "Probability of occurence", x = clim_variable_name)+
    theme_bw()
}



# Inverse Logit
invLogit <- function(x) exp(x)/(1+ exp(x))


#' @title Plot a standardised climate response curve
#' @description The response is standardised by dividing by the area under the curve. 
#' The plot is meant to show the trend of the response, the actual response values are not valuable.
#' @param data A SpatialGridDataFrame of the climate used in the model
#' @param clim_variable The name of a climate variable used in @param data
#' @param x_range A numeric list containing the range of the climate variable
#' @param pars The climate parameters extracted from the model
#'  using the getPars() function available in the Joe_functions.R file
calculate_integral_response_curve <- function(data, clim_variable, x_range, pars) {
  
  xx_st <- (x_range - mean(data@data[[clim_variable]], na.rm = TRUE)) / sd(data@data[[clim_variable]],
                                                                      na.rm = TRUE)
  fn <- function(x, pars, selVar) invLogit(pars["Intercept"] +
                                             x * pars[selVar] + x * x * pars[paste0(selVar, "_quadratic")])
  
  intVal <- integrate(fn, min(xx_st), max(xx_st), pars = pars, selVar = clim_variable)$value
  
  response <- data.frame(covarVal = x_range,
                         meanEst = invLogit(pars["Intercept"] + xx_st * pars[clim_variable] +
                                                               xx_st * xx_st * pars[paste0(clim_variable,
                                                                                           "_quadratic")]) / intVal)
}
# lower = min(c(preds_native$wc2.1_30s_bio_1$covarVal), preds_invasive$wc2.1_30s_bio_1$covarVal)
# 
# response_functions <- function(pars, variable_code) {
#   fn <- function(x, pars, variable_code) invLogit(pars["Intercept"] +
#                                              x * pars[variable_code] +
#                                              x * x * pars[paste0(variable_code, "_quadratic")])
#   i <- integrate(fn, -Inf, Inf, pars = pars, variable_code = variable_code)
#   
# }
# 
# fn <- function(x, pars, variable_code) invLogit(pars["Intercept"] +
#                                                   x * pars[variable_code] +
#                                                   x * x * pars[paste0(variable_code, "_quadratic")])
# 
# i <- integrate(fn, -Inf, Inf, pars = pars_native, variable_code = "wc2.1_30s_bio_5")
# 
# jk <- fn(x, pars = pars_native, "wc2.1_30s_bio_1")
# 
# r_list <- map(variable_codes, response_functions(pars = pars_native))
# 
# x <- 0:20
# 
# f <- function(x) { 
#   y <- function(x) invLogit(0 + 1.3*x - 0.15*x^2)
#   s <- integrate(y, -Inf, Inf)
#   y(x)/s$value
# }
# 
# 
# plot(x, f(x))
# 
# g <- function(x) {
#   y <- function(x) invLogit(-4 + 1.5*x - 0.08*x^2)
#   s <- integrate(y, -Inf, Inf)
#   y(x)/s$value
# }
# plot(x, g(x))
# 
# 
# 
# h <- function(x) pmin(f(x), g(x))
# 
# o <- integrate(h, -Inf, Inf)
# o$value
# 
# 
# f2 <- function(x) { 
#   y <- function(x) invLogit(-4 + 1.6*x - 0.08*x^2)
#   s <- integrate(y, -Inf, Inf)
#   y(x)/s$value
# }
# 
# 
# plot(x, f2(x))
# 
# g2 <- function(x) {
#   y <- function(x) invLogit(-4 + 1.5*x - 0.08*x^2)
#   s <- integrate(y, -Inf, Inf)
#   y(x)/s$value
# }
# plot(x, g2(x))
# 
# 
# 
# h2 <- function(x) pmin(f2(x), g2(x))
# 
# o2 <- integrate(h2, -Inf, Inf)
# o2$value
# 
# 
# f3 <- function(x) { 
#   y <- function(x) invLogit(0 + 0.5*x - 0.18*x^2)
#   s <- integrate(y, -Inf, Inf)
#   y(x)/s$value
# }
# 
# 
# plot(x, f3(x))
# 
# g3 <- function(x) {
#   y <- function(x) invLogit(-4 + 0.48*x - 0.017*x^2)
#   s <- integrate(y, -Inf, Inf)
#   y(x)/s$value
# }
# plot(x, g3(x))
# 
# 
# 
# h3 <- function(x) pmin(f3(x), g3(x))
# 
# o3 <- integrate(h3, -Inf, Inf)
# o3$value

#' @title Plot a standardised climate response curve
#' @description The response is standardised by dividing by the area under the curve. 
#' The plot is meant to show the trend of the response, the actual response values are not valuable.
#' @param data A data frame containing 2 columns: covarVal & meanEst
#' @param clim_variable_name A string containing the name of the climate variable
#'  (will be displayed on the x-axis)
#' @param y_limit A range that defines the y-axis
plot_integral_response_curve <- function(data, clim_variable_name = "Name of the climate variable", y_limit = c(0,1)) {
  
  ggplot(data, aes(covarVal, meanEst))+
    geom_line()+
    ylim(y_limit)+
    labs(y = "Probability of occurence (standardised)", x = clim_variable_name)+
    xlim(min(data$covarVal), max(data$covarVal))+
    theme_bw()
}




compare_response_curves <- function(native_data, invasive_data, clim_variable_name = "Name of the climate variable") {
  
  ggplot(native_data)+
    geom_area(aes(covarVal, meanEst, fill = "native"))+
    geom_vline(xintercept = min(native_data$covarVal), col = "blue")+
    geom_vline(xintercept = max(native_data$covarVal), col = "blue")+
    geom_area(data = invasive_data, alpha = 0.7,
              aes(covarVal, meanEst, fill = "invasive")) +
    geom_vline(xintercept = min(invasive_data$covarVal), col = "red")+
    geom_vline(xintercept = max(invasive_data$covarVal), col = "red")+
    labs(x = clim_variable_name, y = "Occurence Probability",
         fill ="Population") +
    theme_bw()
  
}




