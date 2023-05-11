#' Plot surface maps
#' @param variable A matrix of basis function values.
#' @param year A matrix of basis function values.
#' @param profJulFracofYearAggr A matrix of basis function values.
#' @param profYearAggr A matrix of basis function values.
#' @param profLatAggr A matrix of basis function values.
#' @param profLongAggr A matrix of basis function values.
#' @param profTempAggr A matrix of basis function values.
#' @param profPsalAggr A matrix of basis function values.
#' @import ggplot2

plotSurfaceMap <- function(variable = "temp", year = NA, profJulFracofYearAggr, profYearAggr, profLatAggr, profLongAggr, profTempAggr, profPsalAggr){

  if(is.na(year)){
    relevant_indices <- which(profJulFracofYearAggr > 31 & profJulFracofYearAggr < 60)
  }else{
    relevant_indices <- which(profYearAggr == as.numeric(year) & profJulFracofYearAggr > 31 & profJulFracofYearAggr < 60)
  }

  lat <- profLatAggr[relevant_indices]
  long <- profLongAggr[relevant_indices]
  long_p <- ifelse(long > 180, long - 360, long)
  data_full <- list(profTempAggr, profPsalAggr)[[which(c('temp', 'psal') == variable)]]
  data_sub <- sapply(data_full[relevant_indices], function(x) x[[1]][1])

  if(variable == 'temp'){
    outlier_index <- which(abs(data_sub) > 1)
  }else{
    outlier_index <- which(abs(data_sub) > 0.5)
  }

  lat <- lat[-outlier_index]
  long <- long[-outlier_index]
  long_p <- long_p[-outlier_index]
  data_sub <- data_sub[-outlier_index]

  zlim_range_temp <- c(-1, 1)
  zlim_range_sal <- c(-0.5, 0.5)

  zlim_range <- list(zlim_range_temp, zlim_range_sal)[[which(c('temp', 'psal') == variable)]]

  plot_title <- list('Surface temperature residuals for profile (deg C)', 'Surface salinity residuals for profile (PSU)')[[which(c('temp', 'psal') == variable)]]

  shapefile <- map_data('world')

  ggplot() +
    geom_polygon(data = shapefile, aes(x = shapefile$long, y = shapefile$lat, group = shapefile$group), fill = 'white', color = 'black', linewidth = .2) +
    geom_point(data = data.frame(lat, long_p, data_sub), aes(x = long_p, y = lat, color = data_sub), size = .005, shape = 20) +
    geom_vline(xintercept = c(-175, -30), linetype="dashed",
               color = "#63666A", size=0.5) +
    labs(color = plot_title, x = 'Longitude', y = 'Latitude') +
    scale_color_gradientn(colours = colorRamps::matlab.like(10), limits = zlim_range)  +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.position = 'top',
          legend.title = element_text(size = 13),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold.italic"))
}
