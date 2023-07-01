h_new <- function(scale_horizontal_space, scale_vertical_space,
              lat1d, lon1d, pres1, lat2d, lon2d, pres2, radius){
  return(scale_horizontal_space^2 * calculateDistance(lat1d, lon1d, lat2d, lon2d, radius)^2 +
           scale_vertical_space^2 * (pres1 - pres2)^2)
}

h1 <- function(scale_horizontal_space, lat1d, lon1d, lat2d, lon2d, radius){

  lat1r = deg2rad(lat1d)
  lon1r = deg2rad(lon1d)
  lat2r = deg2rad(lat2d)
  lon2r = deg2rad(lon2d)
  L = lat1r - lat2r
  l = lon1r - lon2r

  con = 4 * scale_horizontal_space^2 * radius^2

  return (con * (sin(L / 2) * cos(L / 2) - sin(lat1r) * cos(lat2r) * (sin(l / 2))^2))
}

h3 <- function(scale_horizontal_space, lat1d, lon1d, lat2d, lon2d, radius){

  lat1r = deg2rad(lat1d)
  lon1r = deg2rad(lon1d)
  lat2r = deg2rad(lat2d)
  lon2r = deg2rad(lon2d)
  l = lon1r - lon2r

  con = 4 * scale_horizontal_space^2 * radius^2

  return(con * sin(l / 2) * cos(l / 2))
}

h33 <- function(scale_horizontal_space, lat1d, lon1d, lat2d, lon2d, radius){

  lat1r = deg2rad(lat1d)
  lon1r = deg2rad(lon1d)
  lat2r = deg2rad(lat2d)
  lon2r = deg2rad(lon2d)
  l = lon1r - lon2r

  con = 2 * scale_horizontal_space^2 * radius^2

  return(con * ((cos(l / 2))^2 - (sin(l / 2))^2))
}

h12 <- function(scale_horizontal_space, lat1d, lon1d, lat2d, lon2d, radius){

  lat1r = deg2rad(lat1d)
  lon1r = deg2rad(lon1d)
  lat2r = deg2rad(lat2d)
  lon2r = deg2rad(lon2d)
  L = lat1r - lat2r
  l = lon1r - lon2r

  con = 4 * scale_horizontal_space^2 * radius^2

  return (con * (-(cos(L / 2))^2 / 2 + (sin(L / 2))^2 / 2 + sin(lat1r) * sin(lat2r) * (sin(l / 2))^2))
}

h13 <- function(scale_horizontal_space, lat1d, lon1d, lat2d, lon2d, radius){

  lat1r = deg2rad(lat1d)
  lon1r = deg2rad(lon1d)
  lat2r = deg2rad(lat2d)
  lon2r = deg2rad(lon2d)
  l = lon1r - lon2r

  con = 4 * scale_horizontal_space^2 * radius^2

  return (-con * sin(lat1r) * sin(l / 2) * cos(l / 2))
}

h4 <- function(scale_vertical_space, pres1, pres2){
  return(2 * scale_vertical_space^2 * (pres1 - pres2))
}

h44 <- function(scale_vertical_space){
  return(2 * scale_vertical_space^2)
}

C1 <- function(stationary_param, nonstationary_param1, nonstationary_param2, lat1d, lon1d, pres1, lat2d, lon2d, pres2, radius){

  a1 = stationary_param[3]
  b1 = stationary_param[4]
  d1 = stationary_param[5]
  a2 = stationary_param[6]
  b2 = stationary_param[7]
  d2 = stationary_param[8]

  c1 = nonstationary_param1
  c2 = nonstationary_param2

  H = h_new(stationary_param[1], stationary_param[2], lat1d, lon1d, pres1, lat2d, lon2d, pres2, radius);
  H1 = h1(stationary_param[1], lat1d, lon1d, lat2d, lon2d, radius);
  H2 = h1(stationary_param[1], lat2d, lon1d, lat1d, lon2d, radius);
  H3 = h3(stationary_param[1], lat1d, lon1d, lat2d, lon2d, radius);
  H4 = h4(stationary_param[2], pres1, pres2);

  return(0.25 * (a1 * a2 * H1 * H2 - b1 * b2 * H3^2 - c1 * c2 * H4^2 - a1 * b2 * H1 * H3
                 + a2 * b1 * H2 * H3 - a1 * c2 * H1 * H4 + a2 * c1 * H2 * H4
                 - b1 * c2 * H3 * H4 - b2 * c1 * H3 * H4) + H * d1 * d2)
}

C2 <- function(stationary_param, nonstationary_param1, nonstationary_param2, lat1d, lon1d, pres1, lat2d, lon2d, pres2, radius){

  a1 = stationary_param[3]
  b1 = stationary_param[4]
  d1 = stationary_param[5]
  a2 = stationary_param[6]
  b2 = stationary_param[7]
  d2 = stationary_param[8]

  nu = stationary_param[9]

  c1 = nonstationary_param1
  c2 = nonstationary_param2

  H12 = h12(stationary_param[1], lat1d, lon1d, lat2d, lon2d, radius)
  H13 = h13(stationary_param[1], lat1d, lon1d, lat2d, lon2d, radius)
  H23 = h13(stationary_param[1], lat2d, lon1d, lat1d, lon2d, radius)
  H33 = h33(stationary_param[1], lat1d, lon1d, lat2d, lon2d, radius)
  H44 = h44(stationary_param[2])

  return(-0.5 * (a1 * a2 * H12 - b1 * b2 * H33 - c1 * c2 * H44 - a1 * b2 * H13 + a2 * b1 * H23) + 2 * nu * d1 * d2)
}

uni_differential <- function(PARAM, fd_eval_mat_loc1, fd_eval_mat_loc2, LAT1D, LON1D, PRES1, LAT2D, LON2D, PRES2, radius){

  sigma_square = PARAM[1]
  SCALE_HORIZONTAL_SPACE = PARAM[2]
  SCALE_VERTICAL_SPACE = PARAM[3]
  smoothness = PARAM[4]

  a1 <- PARAM[5]
  b1 <- PARAM[6]
  d1 <- PARAM[7]

  a2 <- PARAM[8]
  b2 <- PARAM[9]
  d2 <- PARAM[10]

  c1 <- fd_eval_mat_loc1
  c2 <- fd_eval_mat_loc2

  con = 2^(smoothness - 1) * gamma(smoothness)
  con = 1.0 / con
  con = sigma_square * con

  expr <- sqrt(h_new(SCALE_HORIZONTAL_SPACE, SCALE_VERTICAL_SPACE, LAT1D, LON1D, PRES1, LAT2D, LON2D, PRES2, radius))

  STATIONARY_PARAM <- c(SCALE_HORIZONTAL_SPACE, SCALE_VERTICAL_SPACE, a1, b1, d1, a2, b2, d2, smoothness)

  f <- expr^(smoothness - 1) * besselK(expr, smoothness - 1)
  f_prime <- expr^(smoothness - 2) * besselK(expr, smoothness - 2)

  C1_val <- C1(stationary_param = STATIONARY_PARAM, nonstationary_param1 = c1, nonstationary_param2 = c2,
               lat1d = LAT1D, lon1d = LON1D, pres1 = PRES1, lat2d = LAT2D, lon2d = LON2D, pres2 = PRES2, radius)
  C2_val <- C2(stationary_param = STATIONARY_PARAM, nonstationary_param1 = c1, nonstationary_param2 = c2,
               lat1d = LAT1D, lon1d = LON1D, pres1 = PRES1, lat2d = LAT2D, lon2d = LON2D, pres2 = PRES2, radius)

  val <- con * (C1_val * f_prime + C2_val * f + d1 * d2 * (expr^2 * f_prime + 2 * (smoothness - 1) * f))
  diag(val) <- con * (diag(C1_val) + diag(C2_val)) + sigma_square * d1 * d2

  return(val)

}

#' Compute the bivariate differential operator cross-covariance function
#'
#' @description
#' \code{cov_bi_differential} evaluates the nonstationary spatial cross-covariance
#' function model based on the differential operators approach in 3D of the form:
#' \deqn{C_{ij}(L_1, L_2, l_1 - l_2, p_1, p_2) = K_{ij}^{1} \mathcal{M}_{\nu_{ij} - 1} \{ {h(L_1, L_2, l_1 - l_2, p_1 - p_2)}^{1/2} \}}
#' \deqn{\quad \quad \quad \quad \quad \quad + K_{ij}^{2} \mathcal{M}_{\nu_{ij}} \{ {h(L_1, L_2, l_1 - l_2, p_1 - p_2)}^{1/2} \},}
#' for different pairs of locations \eqn{(L_1, l_1, p_1)} and \eqn{(L_2, l_2, p_2)},
#' where \eqn{L} represents the latitude, \eqn{l} the longitude, and \eqn{p} the pressure coordinates, respectively.
#' The forms of \eqn{K_{ij}^{1}} and \eqn{K_{ij}^{2}} can be found in Appendix of Salvana, M. L., & Jun, M. (2022)
#' and \eqn{h(L_1, L_2, l_1 - l_2, p_1 - p_2)} is a distance function of the form:
#' \deqn{h(L_1, L_2, l_1 - l_2, p_1 - p_2) = {a_{h}^2 ch^2(L_1, L_2, l_1 - l_2) + a_{v}^2 (p_1 - p_2)^2},}
#' where \eqn{a_{h}} and \eqn{a_{v}} are the scale parameters in the horizontal and vertical directions,
#' respectively, and \eqn{ch(L_1, L_2, l_1 - l_2)} is the chordal distance with the following formula:
#' \deqn{ch(L_1, L_2, l_1 - l_2) = 2 R \left\{ \sin^2 \left( \frac{L_1 - L_2}{2} \right) + \cos L_1 \cos L_2 \sin^2 \left( \frac{l_1 - l_2}{2} \right) \right\}^{1/2}.}
#' Here \eqn{R} is the radius of the sphere. Note that for global processes,
#' the relevant sphere is the Earth with \eqn{R=6,371} km.
#'
#' @usage cov_bi_differential(location, beta, scale_horizontal, scale_vertical,
#' a1, b1, c1, d1, a2, b2, c2, d2, radius)
#'
#' @param location An \eqn{n \times 3} matrix of coordinates.
#' @param beta A numeric constant indicating the colocated correlation parameter.
#' @param scale_horizontal A numeric constant indicating the horizontal scale parameter.
#' @param scale_vertical A numeric constant indicating the vertical scale parameter.
#' @param a1 A numeric constant indicating the anisotropy in latitude parameter associated with variable 1.
#' @param b1 A numeric constant indicating the anisotropy in longitude parameter associated with variable 1.
#' @param c1 A numeric vector indicating the nonstationary parameters with depth associated with variable 1.
#' @param d1 A numeric constant indicating the variance parameter from the fully isotropic component associated with variable 1.
#' @param a2 A numeric constant indicating the anisotropy in latitude parameter associated with variable 2.
#' @param b2 A numeric constant indicating the anisotropy in longitude parameter associated with variable 2.
#' @param c2 A numeric vector indicating the nonstationary parameters with depth associated with variable 2.
#' @param d2 A numeric constant indicating the variance parameter from the fully isotropic component associated with variable 2.
#' @param radius A numeric constant indicating the radius of the sphere.
#'
#' @return A cross-covariance matrix of dimension \eqn{2 n \times 2n}.
#'
#' @author Mary Lai Salvana \email{yourlainess@gmail.com}
#'
#' @references Salvana, M. L., & Jun, M. (2022). 3D Bivariate Spatial Modelling of Argo Ocean Temperature and Salinity Profiles. \emph{arXiv preprint arXiv:2210.11611}.
#'
#' @examples
#'
#' library(dplyr)
#'
#' x <- seq(0, 1, length.out = 10)
#' y <- seq(0, 1, length.out = 10)
#'
#' loc2d <- expand.grid(x, y) %>% as.matrix()
#' depth <- seq(0, 1, length.out = 10)
#' loc3d <- cbind(rep(loc2d[, 1], each = length(depth)), rep(loc2d[, 2], each = length(depth)), depth)
#'
#' earthRadiusKm = 6371
#'
#' BETA = 0.5
#' SCALE_HORIZONTAL = 0.03
#' SCALE_VERTICAL = 0.3
#' A1 = A2 = 0.00001
#' B1 = B2 = 0.00001
#' C1 = sin((loc3d[1:10, 3] + 0.1) * pi / 0.5)
#' C2 = cos((loc3d[1:10, 3] + 0.1) * pi / 0.5)
#' D1 = D2 = 0
#'
#' cov_mat <- cov_bi_differential(location = loc3d, beta = BETA,
#'                                scale_horizontal = SCALE_HORIZONTAL,
#'                                scale_vertical = SCALE_VERTICAL,
#'                                a1 = A1, b1 = B1, c1 = C1, d1 = D1,
#'                                a2 = A2, b2 = B2, c2 = C2, d2 = D2,
#'                                radius = earthRadiusKm)
#'
#'
#' @export
cov_bi_differential <- function(location, beta, scale_horizontal, scale_vertical, a1, b1, c1, d1, a2, b2, c2, d2, radius){

  LAT1D <- matrix(location[, 2], nrow(location), nrow(location), byrow = F)
  LON1D <- matrix(location[, 1], nrow(location), nrow(location), byrow = F)
  PRES1 <- matrix(location[, 3], nrow(location), nrow(location), byrow = F)
  LAT2D <- matrix(location[, 2], nrow(location), nrow(location), byrow = T)
  LON2D <- matrix(location[, 1], nrow(location), nrow(location), byrow = T)
  PRES2 <- matrix(location[, 3], nrow(location), nrow(location), byrow = T)

  fd_eval_mat_loc1 <- matrix(c1, nrow(location), nrow(location), byrow = F)
  fd_eval_mat_loc2 <- matrix(c1, nrow(location), nrow(location), byrow = T)

  fd_eval2_mat_loc1 <- matrix(c2, nrow(location), nrow(location), byrow = F)
  fd_eval2_mat_loc2 <- matrix(c2, nrow(location), nrow(location), byrow = T)

  PARAM <- c(1, scale_horizontal, scale_vertical, 2, a1, b1, d1, a1, b1, d1)

  cov_val <- uni_differential(PARAM, fd_eval_mat_loc1, fd_eval_mat_loc2, LAT1D, LON1D, PRES1, LAT2D, LON2D, PRES2, radius)

  PARAM <- c(1, scale_horizontal, scale_vertical, 2, a2, b2, d2, a2, b2, d2)

  cov_val2 <- uni_differential(PARAM, fd_eval2_mat_loc1, fd_eval2_mat_loc2, LAT1D, LON1D, PRES1, LAT2D, LON2D, PRES2, radius)

  PARAM <- c(beta, scale_horizontal, scale_vertical, 2, a1, b1, d1, a2, b2, d2)

  cov_val3 <- uni_differential(PARAM, fd_eval_mat_loc1, fd_eval2_mat_loc2, LAT1D, LON1D, PRES1, LAT2D, LON2D, PRES2, radius)

  Sigma <- rbind(cbind(cov_val, cov_val3), cbind(t(cov_val3), cov_val2))

  return(Sigma)

}

#' Compute the bivariate differential operator cross-covariance function in parallel
#'
#' @description
#' \code{cov_bi_differential_parallel} evaluates the nonstationary spatial cross-covariance
#' function model based on the differential operators approach in 3D of the form:
#' \deqn{C_{ij}(L_1, L_2, l_1 - l_2, p_1, p_2) = K_{ij}^{1} \mathcal{M}_{\nu_{ij} - 1} \{ {h(L_1, L_2, l_1 - l_2, p_1 - p_2)}^{1/2} \}}
#' \deqn{\quad \quad \quad \quad \quad \quad + K_{ij}^{2} \mathcal{M}_{\nu_{ij}} \{ {h(L_1, L_2, l_1 - l_2, p_1 - p_2)}^{1/2} \},}
#' for different pairs of locations \eqn{(L_1, l_1, p_1)} and \eqn{(L_2, l_2, p_2)},
#' where \eqn{L} represents the latitude, \eqn{l} the longitude, and \eqn{p} the pressure coordinates, respectively.
#' The forms of \eqn{K_{ij}^{1}} and \eqn{K_{ij}^{2}} can be found in Appendix of Salvana, M. L., & Jun, M. (2022)
#' and \eqn{h(L_1, L_2, l_1 - l_2, p_1 - p_2)} is a distance function of the form:
#' \deqn{h(L_1, L_2, l_1 - l_2, p_1 - p_2) = {a_{h}^2 ch^2(L_1, L_2, l_1 - l_2) + a_{v}^2 (p_1 - p_2)^2},}
#' where \eqn{a_{h}} and \eqn{a_{v}} are the scale parameters in the horizontal and vertical directions,
#' respectively, and \eqn{ch(L_1, L_2, l_1 - l_2)} is the chordal distance with the following formula:
#' \deqn{ch(L_1, L_2, l_1 - l_2) = 2 R \left\{ \sin^2 \left( \frac{L_1 - L_2}{2} \right) + \cos L_1 \cos L_2 \sin^2 \left( \frac{l_1 - l_2}{2} \right) \right\}^{1/2}.}
#' Here \eqn{R} is the radius of the sphere. Note that for global processes,
#' the relevant sphere is the Earth with \eqn{R=6,371} km.
#'
#' @usage cov_bi_differential_parallel(location, beta, scale_horizontal, scale_vertical,
#' a1, b1, c1, d1, a2, b2, c2, d2, radius, num_processors)
#'
#' @param location An \eqn{n \times 3} matrix of coordinates.
#' @param beta A numeric constant indicating the colocated correlation parameter.
#' @param scale_horizontal A numeric constant indicating the horizontal scale parameter.
#' @param scale_vertical A numeric constant indicating the vertical scale parameter.
#' @param a1 A numeric constant indicating the anisotropy in latitude parameter associated with variable 1.
#' @param b1 A numeric constant indicating the anisotropy in longitude parameter associated with variable 1.
#' @param c1 A numeric vector indicating the nonstationary parameters with depth associated with variable 1.
#' @param d1 A numeric constant indicating the variance parameter from the fully isotropic component associated with variable 1.
#' @param a2 A numeric constant indicating the anisotropy in latitude parameter associated with variable 2.
#' @param b2 A numeric constant indicating the anisotropy in longitude parameter associated with variable 2.
#' @param c2 A numeric vector indicating the nonstationary parameters with depth associated with variable 2.
#' @param d2 A numeric constant indicating the variance parameter from the fully isotropic component associated with variable 2.
#' @param radius A numeric constant indicating the radius of the sphere.
#' @param num_processors A numeric constant indicating the number of available processors; must be a perfect square.
#'
#' @import pbdBASE
#'
#' @return A cross-covariance matrix of dimension \eqn{2 n \times 2n}.
#'
#' @author Mary Lai Salvana \email{yourlainess@gmail.com}
#'
#' @references Salvana, M. L., & Jun, M. (2022). 3D Bivariate Spatial Modelling of Argo Ocean Temperature and Salinity Profiles. \emph{arXiv preprint arXiv:2210.11611}.
#'
#' @examples
#'
#' library(dplyr)
#'
#' x <- seq(0, 1, length.out = 10)
#' y <- seq(0, 1, length.out = 10)
#'
#' loc2d <- expand.grid(x, y) %>% as.matrix()
#' depth <- seq(0, 1, length.out = 10)
#' loc3d <- cbind(rep(loc2d[, 1], each = length(depth)), rep(loc2d[, 2], each = length(depth)), depth)
#'
#' earthRadiusKm = 6371
#'
#' BETA = 0.5
#' SCALE_HORIZONTAL = 0.03
#' SCALE_VERTICAL = 0.3
#' A1 = A2 = 0.00001
#' B1 = B2 = 0.00001
#' C1 = sin((loc3d[1:10, 3] + 0.1) * pi / 0.5)
#' C2 = cos((loc3d[1:10, 3] + 0.1) * pi / 0.5)
#' D1 = D2 = 0
#'
#' cov_mat <- cov_bi_differential_parallel(location = loc3d, beta = BETA,
#'                                scale_horizontal = SCALE_HORIZONTAL,
#'                                scale_vertical = SCALE_VERTICAL,
#'                                a1 = A1, b1 = B1, c1 = C1, d1 = D1,
#'                                a2 = A2, b2 = B2, c2 = C2, d2 = D2,
#'                                radius = earthRadiusKm,
#'                                num_processors = 4)
#'
#'
#' @export
cov_bi_differential_parallel <- function(location, beta, scale_horizontal, scale_vertical, a1, b1, c1, d1, a2, b2, c2, d2, radius, num_processors){

  init.grid(NPROW = sqrt(num_processors), NPCOL = sqrt(num_processors))
  bldim <- c(sqrt(num_processors), sqrt(num_processors))

  LAT1D <- matrix(location[, 2], nrow(location), nrow(location), byrow = F)
  LON1D <- matrix(location[, 1], nrow(location), nrow(location), byrow = F)
  PRES1 <- matrix(location[, 3], nrow(location), nrow(location), byrow = F)
  LAT2D <- matrix(location[, 2], nrow(location), nrow(location), byrow = T)
  LON2D <- matrix(location[, 1], nrow(location), nrow(location), byrow = T)
  PRES2 <- matrix(location[, 3], nrow(location), nrow(location), byrow = T)

  fd_eval_mat_loc1 <- matrix(c1, nrow(location), nrow(location), byrow = F)
  fd_eval_mat_loc2 <- matrix(c1, nrow(location), nrow(location), byrow = T)

  fd_eval2_mat_loc1 <- matrix(c2, nrow(location), nrow(location), byrow = F)
  fd_eval2_mat_loc2 <- matrix(c2, nrow(location), nrow(location), byrow = T)

  PARAM <- c(1, scale_horizontal, scale_vertical, 2, a1, b1, d1, a1, b1, d1)

  cov_val <- uni_differential(PARAM, fd_eval_mat_loc1, fd_eval_mat_loc2, LAT1D, LON1D, PRES1, LAT2D, LON2D, PRES2, radius)

  PARAM <- c(1, scale_horizontal, scale_vertical, 2, a2, b2, d2, a2, b2, d2)

  cov_val2 <- uni_differential(PARAM, fd_eval2_mat_loc1, fd_eval2_mat_loc2, LAT1D, LON1D, PRES1, LAT2D, LON2D, PRES2, radius)

  PARAM <- c(beta, scale_horizontal, scale_vertical, 2, a1, b1, d1, a2, b2, d2)

  cov_val3 <- uni_differential(PARAM, fd_eval_mat_loc1, fd_eval2_mat_loc2, LAT1D, LON1D, PRES1, LAT2D, LON2D, PRES2, radius)

  Sigma <- rbind(cbind(cov_val, cov_val3), cbind(t(cov_val3), cov_val2))

  return(Sigma)

  finalize()

}

#' Step 1: WLS estimation of the parameters of the bivariate differential operator cross-covariance function
#' @param empirical_values A matrix of empirical values
#' @param location An nx3 matrix of coordinates.
#' @param init_beta A number for colocated correlation parameter.
#' @param init_beta_fix An indicator whether colocated correlation parameter should be estimated.
#' @param init_scale_horizontal A number for the horizontal scale parameter.
#' @param init_scale_vertical A number for the vertical scale parameter.
#' @param init_scale_horizontal_fix An indicator whether horizontal scale parameter should be estimated.
#' @param init_scale_vertical_fix An indicator whether vertical scale parameter should be estimated.
#' @param init_a1 A number for the anisotropy parameter in Latitude associated with variable 1.
#' @param init_b1 A number for the anisotropy parameter in longitude associated with variable 1.
#' @param init_c1 A number for vector for the nonstationarity parameter in depth associated with variable 1.
#' @param init_d1 A number for the variance parameter of the fully isotropic associated with variable 1.
#' @param init_a1_fix An indicator whether a1 parameter should be estimated.
#' @param init_b1_fix An indicator whether b1 parameter should be estimated.
#' @param init_a2 A number for the anisotropy parameter in Latitude associated with variable 2.
#' @param init_b2 A number for the anisotropy parameter in longitude associated with variable 2.
#' @param init_c2 A number for vector for the nonstationarity parameter in depth associated with variable 2.
#' @param init_d2 A number for the variance parameter of the fully isotropic associated with variable 2.
#' @param init_a2_fix An indicator whether a2 parameter should be estimated.
#' @param init_b2_fix An indicator whether b2 parameter should be estimated.
#' @param radius A number for the radius of the sphere.
#' @param basis1 A matrix of basis function values for variable 1.
#' @param basis2 A matrix of basis function values for variable 2.
#' @param splines_degree A number indicating the degree of the splines.
#' @param knots1 A vector of knot locations for variable 1.
#' @param knots2 A vector of knot locations for variable 2.
#' @import stats
#' @import mvtnorm
#' @return A vector of estimated parameter values.
#' @export
est_bi_differential_wls <- function(empirical_values, location, init_beta, init_scale_horizontal, init_scale_vertical, init_scale_horizontal_fix = F, init_scale_vertical_fix = F, init_a1, init_b1, init_c1, init_d1 = NULL, init_a1_fix = F, init_b1_fix = F, init_a2, init_b2, init_c2, init_d2 = NULL, init_a2_fix = F, init_b2_fix = F, init_beta_fix = F, radius, basis1, basis2, splines_degree = 4, knots1, knots2){

  NEGLOGLIK <- function(theta, empirical_values, basis1, basis2){

    print(paste("theta = c(", paste(round(theta, 8), collapse=","), ")", sep = ''))

    nb1 <- ncol(basis1)
    nb2 <- ncol(basis2)

    if(!init_beta_fix){
      BETA <- theta[1]

      if(init_a1_fix){
        if(init_scale_horizontal_fix & init_scale_vertical_fix){

          if(theta[1] < 0 | theta[1] > 1){
            return(Inf)
          }

          SCALE_HORIZONTAL <- init_scale_horizontal
          SCALE_VERTICAL <- init_scale_vertical
          A1 <- init_a1
          B1 <- init_b1
          C1_coef <- theta[1 + 1:nb1]

          if(is.null(init_d1) | is.null(init_d2)){
            D1 <- 0
            A2 <- init_a2
            B2 <- init_b2
            C2_coef <- theta[nb1 + 1 + 1:nb2]
            D2 <- 0
          }else{
            D1 <- theta[nb1 + 2]
            A2 <- init_a2
            B2 <- init_b2
            C2_coef <- theta[nb1 + 2 + 1:nb2]
            D2 <- theta[nb1 + nb2 + 3]
          }
        }else{

          if(theta[1] < 0 | theta[1] > 1){
            return(Inf)
          }

          SCALE_HORIZONTAL <- exp(theta[2])
          SCALE_VERTICAL <- exp(theta[3])
          A1 <- init_a1
          B1 <- init_b1
          C1_coef <- theta[3 + 1:nb1]

          if(is.null(init_d1) | is.null(init_d2)){
            D1 <- 0
            A2 <- init_a2
            B2 <- init_b2
            C2_coef <- theta[nb1 + 3 + 1:nb2]
            D2 <- 0
          }else{
            D1 <- theta[nb1 + 4]
            A2 <- init_a2
            B2 <- init_b2
            C2_coef <- theta[nb1 + 4 + 1:nb2]
            D2 <- theta[nb1 + nb2 + 5]
          }
        }
      }else{
        if(init_scale_horizontal_fix & init_scale_vertical_fix){

          if(theta[1] < 0 | theta[1] > 1 | theta[2] < 0){
            return(Inf)
          }

          SCALE_HORIZONTAL <- init_scale_horizontal
          SCALE_VERTICAL <- init_scale_vertical
          A1 <- theta[2]
          B1 <- theta[3]
          C1_coef <- theta[3 + 1:nb1]

          if(is.null(init_d1) | is.null(init_d2)){
            D1 <- 0
            A2 <- theta[nb1 + 4]
            B2 <- theta[nb1 + 5]
            C2_coef <- theta[nb1 + 5 + 1:nb2]
            D2 <- 0
          }else{
            D1 <- theta[nb1 + 4]
            A2 <- theta[nb1 + 5]
            B2 <- theta[nb1 + 6]
            C2_coef <- theta[nb1 + 6 + 1:nb2]
            D2 <- theta[nb1 + nb2 + 7]
          }
        }else{

          if(theta[1] < 0 | theta[1] > 1 | theta[4] < 0){
            return(Inf)
          }

          SCALE_HORIZONTAL <- exp(theta[2])
          SCALE_VERTICAL <- exp(theta[3])
          A1 <- theta[4]
          B1 <- theta[5]
          C1_coef <- theta[5 + 1:nb1]

          if(is.null(init_d1) | is.null(init_d2)){
            D1 <- 0
            A2 <- theta[nb1 + 6]
            B2 <- theta[nb1 + 7]
            C2_coef <- theta[nb1 + 7 + 1:nb2]
            D2 <- 0
          }else{
            D1 <- theta[nb1 + 6]
            A2 <- theta[nb1 + 7]
            B2 <- theta[nb1 + 8]
            C2_coef <- theta[nb1 + 8 + 1:nb2]
            D2 <- theta[nb1 + nb2 + 9]
          }
        }
      }

    }else{

      BETA <- init_beta

      if(init_a1_fix){
        if(init_scale_horizontal_fix & init_scale_vertical_fix){

          SCALE_HORIZONTAL <- init_scale_horizontal
          SCALE_VERTICAL <- init_scale_vertical
          A1 <- init_a1
          B1 <- init_b1
          C1_coef <- theta[1:nb1]

          if(is.null(init_d1) | is.null(init_d2)){
            D1 <- 0
            A2 <- init_a2
            B2 <- init_b2
            C2_coef <- theta[nb1 + 1:nb2]
            D2 <- 0
          }else{
            D1 <- theta[nb1 + 1]
            A2 <- init_a2
            B2 <- init_b2
            C2_coef <- theta[nb1 + 1 + 1:nb2]
            D2 <- theta[nb1 + nb2 + 2]
          }
        }else{

          SCALE_HORIZONTAL <- exp(theta[1])
          SCALE_VERTICAL <- exp(theta[2])
          A1 <- init_a1
          B1 <- init_b1
          C1_coef <- theta[2 + 1:nb1]

          if(is.null(init_d1) | is.null(init_d2)){
            D1 <- 0
            A2 <- init_a2
            B2 <- init_b2
            C2_coef <- theta[nb1 + 2 + 1:nb2]
            D2 <- 0
          }else{
            D1 <- theta[nb1 + 3]
            A2 <- init_a2
            B2 <- init_b2
            C2_coef <- theta[nb1 + 3 + 1:nb2]
            D2 <- theta[nb1 + nb2 + 4]
          }
        }
      }else{
        if(init_scale_horizontal_fix & init_scale_vertical_fix){

          if(theta[1] < 0){
            return(Inf)
          }

          SCALE_HORIZONTAL <- init_scale_horizontal
          SCALE_VERTICAL <- init_scale_vertical
          A1 <- theta[1]
          B1 <- theta[2]
          C1_coef <- theta[2 + 1:nb1]

          if(is.null(init_d1) | is.null(init_d2)){
            D1 <- 0
            A2 <- theta[nb1 + 3]
            B2 <- theta[nb1 + 4]
            C2_coef <- theta[nb1 + 4 + 1:nb2]
            D2 <- 0
          }else{
            D1 <- theta[nb1 + 3]
            A2 <- theta[nb1 + 4]
            B2 <- theta[nb1 + 5]
            C2_coef <- theta[nb1 + 5 + 1:nb2]
            D2 <- theta[nb1 + nb2 + 6]
          }
        }else{

          if(theta[3] < 0){
            return(Inf)
          }

          SCALE_HORIZONTAL <- exp(theta[1])
          SCALE_VERTICAL <- exp(theta[2])
          A1 <- theta[3]
          B1 <- theta[4]
          C1_coef <- theta[4 + 1:nb1]

          if(is.null(init_d1) | is.null(init_d2)){
            D1 <- 0
            A2 <- theta[nb1 + 5]
            B2 <- theta[nb1 + 6]
            C2_coef <- theta[nb1 + 6 + 1:nb2]
            D2 <- 0
          }else{
            D1 <- theta[nb1 + 5]
            A2 <- theta[nb1 + 6]
            B2 <- theta[nb1 + 7]
            C2_coef <- theta[nb1 + 7 + 1:nb2]
            D2 <- theta[nb1 + nb2 + 8]
          }
        }
      }
    }

    C1 <- basis1 %*% matrix(C1_coef, ncol = 1)
    C2 <- basis2 %*% matrix(C2_coef, ncol = 1)

    cov_mat <- cov_bi_differential(location = location, beta = BETA,
                        scale_horizontal = SCALE_HORIZONTAL, scale_vertical = SCALE_VERTICAL,
                        a1 = A1, b1 = B1, c1 = C1, d1 = D1, a2 = A2, b2 = B2, c2 = C2, d2 = D2,
                        radius = radius)
    #cor_mat <- cov2cor(cov_mat)

    #print(cor_mat[1, 1:55])

    variance1 <- diag(cov_mat[1:nrow(location), 1:nrow(location)])
    variance2 <- diag(cov_mat[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)])

    emp_variance1 <- diag(empirical_values[1:nrow(location), 1:nrow(location)])
    emp_variance2 <- diag(empirical_values[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)])

    covariance1 <- cov_mat[1:nrow(location), 1:nrow(location)]
    covariance2 <- cov_mat[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)]
    covariance12 <- diag(cov_mat[1:nrow(location), nrow(location) + 1:nrow(location)])
    correlation12 <- covariance12 / outer(sqrt(variance1), sqrt(variance2), '*')

    emp_covariance1 <- empirical_values[1:nrow(location), 1:nrow(location)]
    emp_covariance2 <- empirical_values[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)]
    emp_covariance12 <- empirical_values[1:nrow(location), nrow(location) + 1:nrow(location)]
    emp_correlation12 <- emp_covariance12 / outer(sqrt(emp_variance1), sqrt(emp_variance2), '*')

    #cholmat <- tryCatch(chol(cov_mat), error = function(a) numeric(0))

    #if(cor_mat[1, 5] < 0.9 | cor_mat[nrow(location) + 1, nrow(location) + 5] < 0.9 | D1 < 0 | D2 < 0){
    #if(cor_mat[1, 2] < 0.95 | cor_mat[1, 5] > 0.6 | cor_mat[1, 51] > 1e-15 | cor_mat[nrow(location) + 1, nrow(location) + 2] < 0.95 | cor_mat[nrow(location) + 1, nrow(location) + 5] > 0.6 | cor_mat[nrow(location) + 1, nrow(location) + 51] > 1e-15 | D1 < 0 | D2 < 0){
    #if(cor_mat[1, 2] < 0.95 | cor_mat[1, 3] > 0.75 | cor_mat[1, 5] > 0.6 | cor_mat[1, 6] < 0.4 | cor_mat[1, 51] > 1e-15 | D1 < 0 | D2 < 0){
    #if(D1 < 0 | D2 < 0){
    #if(length(cholmat) == 0 | D1 < 0 | D2 < 0){
      #return(Inf)
    #}else{

      plot_bi_differential(location = location[1:50, ], est_beta = BETA,
                           est_scale_horizontal = SCALE_HORIZONTAL, est_scale_vertical = SCALE_VERTICAL,
                           est_a1 = A1, est_b1 = B1, est_c1 = C1_coef, est_d1 = D1, est_a2 = A2, est_b2 = B2, est_c2 = C2_coef, est_d2 = D2,
                           basis1 = basis1[1:50, ], basis2 = basis2[1:50, ], radius = radius)

      #out <- sum((variance1 - empirical_values[, 1])^2 * 100 + (variance2 - empirical_values[, 2])^2 * 1000 + (correlation12 - empirical_values[, 3])^2 * 50)


      #emp_covariance12 <- empirical_values[1:nrow(location), nrow(location) + 1:nrow(location)]
      #emp_correlation <- emp_covariance12 / outer(sqrt(emp_variance1), sqrt(emp_variance2), '*')


      #theo_covariance12 <- cov_mat[1:nrow(location), nrow(location) + 1:nrow(location)]
      #theo_correlation <- theo_covariance12 / outer(sqrt(theo_variance1), sqrt(theo_variance2), '*')

      #out <- sum((emp_variance1 - theo_variance1)^2 * 100 +
      #             (emp_variance2 - theo_variance2)^2 * 1000 +
      #             (emp_correlation - theo_correlation)^2 * 50)

      out <- sum((emp_covariance1 - covariance1)^2 * 10 +
                   (emp_covariance2 - covariance2)^2 * 10000 +
                   (emp_correlation12 - correlation12)^2 * 0)

      return(out)
    #}
  }

  if(!init_beta_fix){
    if(init_a1_fix){
      if(init_scale_horizontal_fix & init_scale_vertical_fix){
        if(is.null(init_d1) | is.null(init_d2)){
          fit <- optim(par = c(init_beta, init_c1, init_c2), fn = NEGLOGLIK, empirical_values = empirical_values, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = 500))
        }else{
          fit <- optim(par = c(init_beta, init_c1, init_d1, init_c2, init_d2), fn = NEGLOGLIK, empirical_values = empirical_values, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = 500))
        }
      }else{
        if(is.null(init_d1) | is.null(init_d2)){
          fit <- optim(par = c(init_beta, init_scale_horizontal, init_scale_vertical, init_c1, init_c2), fn = NEGLOGLIK, empirical_values = empirical_values, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = 500))
        }else{
          fit <- optim(par = c(init_beta, init_scale_horizontal, init_scale_vertical, init_c1, init_d1, init_c2, init_d2), fn = NEGLOGLIK, empirical_values = empirical_values, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = 500))
        }
      }
    }else{
      if(init_scale_horizontal_fix & init_scale_vertical_fix){
        if(is.null(init_d1) | is.null(init_d2)){
          fit <- optim(par = c(init_beta, init_a1, init_b1, init_c1, init_a2, init_b2, init_c2), fn = NEGLOGLIK, empirical_values = empirical_values, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = 500))
        }else{
          fit <- optim(par = c(init_beta, init_a1, init_b1, init_c1, init_d1, init_a2, init_b2, init_c2, init_d2), fn = NEGLOGLIK, empirical_values = empirical_values, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = 500))
        }
      }else{
        if(is.null(init_d1) | is.null(init_d2)){
          fit <- optim(par = c(init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1, init_a2, init_b2, init_c2), fn = NEGLOGLIK, empirical_values = empirical_values, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = 500))
        }else{
          fit <- optim(par = c(init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1, init_d1, init_a2, init_b2, init_c2, init_d2), fn = NEGLOGLIK, empirical_values = empirical_values, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = 500))
        }
      }
    }
  }else{
    if(init_a1_fix){
      if(init_scale_horizontal_fix & init_scale_vertical_fix){
        if(is.null(init_d1) | is.null(init_d2)){
          fit <- optim(par = c(init_c1, init_c2), fn = NEGLOGLIK, empirical_values = empirical_values, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = 500))
        }else{
          fit <- optim(par = c(init_c1, init_d1, init_c2, init_d2), fn = NEGLOGLIK, empirical_values = empirical_values, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = 500))
        }
      }else{
        if(is.null(init_d1) | is.null(init_d2)){
          fit <- optim(par = c(init_scale_horizontal, init_scale_vertical, init_c1, init_c2), fn = NEGLOGLIK, empirical_values = empirical_values, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = 500))
        }else{
          fit <- optim(par = c(init_scale_horizontal, init_scale_vertical, init_c1, init_d1, init_c2, init_d2), fn = NEGLOGLIK, empirical_values = empirical_values, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = 500))
        }
      }
    }else{
      if(init_scale_horizontal_fix & init_scale_vertical_fix){
        if(is.null(init_d1) | is.null(init_d2)){
          fit <- optim(par = c(init_a1, init_b1, init_c1, init_a2, init_b2, init_c2), fn = NEGLOGLIK, empirical_values = empirical_values, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = 500))
        }else{
          fit <- optim(par = c(init_a1, init_b1, init_c1, init_d1, init_a2, init_b2, init_c2, init_d2), fn = NEGLOGLIK, empirical_values = empirical_values, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = 500))
        }
      }else{
        if(is.null(init_d1) | is.null(init_d2)){
          fit <- optim(par = c(init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1, init_a2, init_b2, init_c2), fn = NEGLOGLIK, empirical_values = empirical_values, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = 500))
        }else{
          fit <- optim(par = c(init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1, init_d1, init_a2, init_b2, init_c2, init_d2), fn = NEGLOGLIK, empirical_values = empirical_values, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = 500))
        }
      }
    }
  }

  for(aa in 1:20){
    fit <- optim(par = fit$par, fn = NEGLOGLIK, empirical_values = empirical_values, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = 500))
  }

  est_theta <- fit$par

  return(est_theta)
}

#' Maximum Likelihood Estimation of the parameters of the bivariate differential operator cross-covariance function
#'
#' @description
#' \code{est_bi_differential} performs maximum likelihood estimation using a
#' Newton-type algorithm to carry out a minimization of the negative of the
#' Gaussian log-likelihood function:
#' \deqn{l(\boldsymbol{\theta}) = -\frac{q n}{2} \log (2 \pi) - \frac{1}{2} \log |\boldsymbol{\Sigma}(\boldsymbol{\theta})| - \frac{1}{2} \mathbf{Z}^{\top} \boldsymbol{\Sigma}(\boldsymbol{\theta})^{-1} \mathbf{Z},}
#' Here \eqn{q} is the number of parameters and \eqn{\boldsymbol{\Sigma}}.
#'
#' @usage est_bi_differential(residuals, location, init_beta, init_scale_horizontal, 
#' init_scale_vertical, init_scale_horizontal_fix, init_scale_vertical_fix,
#' init_a1, init_b1, init_c1, init_d1, init_a1_fix, init_b1_fix, 
#' init_a2, init_b2, init_c2, init_d2, init_a2_fix, init_b2_fix, 
#' init_beta_fix, radius, basis1, nb1, basis2, nb2, splines_degree,
#' knots1, knots2, MAXIT, RERUNS, STEPMAX)
#'
#' @param residuals A vector of residuals
#' @param location An nx3 matrix of coordinates.
#' @param init_beta A number for colocated correlation parameter.
#' @param init_beta_fix An indicator whether colocated correlation parameter should be estimated.
#' @param init_scale_horizontal A number for the horizontal scale parameter.
#' @param init_scale_vertical A number for the vertical scale parameter.
#' @param init_scale_horizontal_fix An indicator whether horizontal scale parameter should be estimated.
#' @param init_scale_vertical_fix An indicator whether vertical scale parameter should be estimated.
#' @param init_a1 A number for the anisotropy parameter in Latitude associated with variable 1.
#' @param init_b1 A number for the anisotropy parameter in longitude associated with variable 1.
#' @param init_c1 A number for vector for the nonstationarity parameter in depth associated with variable 1.
#' @param init_d1 A number for the variance parameter of the fully isotropic associated with variable 1.
#' @param init_a1_fix An indicator whether a1 parameter should be estimated.
#' @param init_b1_fix An indicator whether b1 parameter should be estimated.
#' @param init_a2 A number for the anisotropy parameter in Latitude associated with variable 2.
#' @param init_b2 A number for the anisotropy parameter in longitude associated with variable 2.
#' @param init_c2 A number for vector for the nonstationarity parameter in depth associated with variable 2.
#' @param init_d2 A number for the variance parameter of the fully isotropic associated with variable 2.
#' @param init_a2_fix An indicator whether a1 parameter should be estimated.
#' @param init_b2_fix An indicator whether b1 parameter should be estimated.
#' @param radius A number for the radius of the sphere.
#' @param basis1 A matrix of basis function values for variable 1.
#' @param nb1 A number indicating the number of bases for variable 1.
#' @param basis2 A matrix of basis function values for variable 2.
#' @param nb2 A number indicating the number of bases for variable 2.
#' @param splines_degree A number indicating the degree of the splines.
#' @param knots1 A vector of knot locations for variable 1.
#' @param knots2 A vector of knot locations for variable 2.
#' @param MAXIT A number indicating the maximum number of iterations for optim.
#' @param RERUNS A number indicating the number of times optim is re-run from previous MLE.
#' @param STEPMAX A number indicating the stepmax of nlm.
#' @import stats
#' @import mvtnorm
#' @return A vector of estimated parameter values.
#' @export
est_bi_differential <- function(residuals, location, init_beta, init_scale_horizontal, init_scale_vertical, init_scale_horizontal_fix = F, init_scale_vertical_fix = F, init_a1, init_b1, init_c1, init_d1 = NULL, init_a1_fix = F, init_b1_fix = F, init_a2, init_b2, init_c2, init_d2 = NULL, init_a2_fix = F, init_b2_fix = F, init_beta_fix = F, radius, basis1, nb1 = ncol(basis1), basis2, nb2 = ncol(basis2), splines_degree = 4, knots1, knots2, MAXIT = 2000, RERUNS = 20, STEPMAX = 1){

  NEGLOGLIK <- function(theta, residuals, basis1, basis2){

    print(paste("theta = c(", paste(round(theta, 8), collapse=","), ")", sep = ''))

    k <- length(theta)

    if(!init_beta_fix){

      if(theta[1] < 0 | theta[1] > 1){
        return(Inf)
      }

      BETA <- theta[1]

      if(init_a1_fix){

        A1 <- init_a1
        B1 <- init_b1
        A2 <- init_a2
        B2 <- init_b2

        if(init_scale_horizontal_fix & init_scale_vertical_fix){

          SCALE_HORIZONTAL <- init_scale_horizontal
          SCALE_VERTICAL <- init_scale_vertical

          if(splines_degree == 0){
            C1_coef <- theta[2]
            if(is.null(init_d1) | is.null(init_d2)){
              D1 <- 0
              C2_coef <- theta[3]
              D2 <- 0
            }else{
              D1 <- theta[3]
              C2_coef <- theta[4]
              D2 <- theta[5]
            }
          }else if(splines_degree > 0){
            C1_coef <- theta[1 + 1:nb1]

            if(is.null(init_d1) | is.null(init_d2)){
              D1 <- 0
              C2_coef <- theta[nb1 + 1 + 1:nb2]
              D2 <- 0
            }else{
              D1 <- theta[nb1 + 2]
              C2_coef <- theta[nb1 + 2 + 1:nb2]
              D2 <- theta[nb1 + nb2 + 3]
            }
          }
        }else if(!init_scale_horizontal_fix & !init_scale_vertical_fix){

          SCALE_HORIZONTAL <- exp(theta[2])
          SCALE_VERTICAL <- exp(theta[3])

          if(splines_degree == 0){

            C1_coef <- theta[4]

            if(is.null(init_d1) | is.null(init_d2)){
              D1 <- 0
              C2_coef <- theta[5]
              D2 <- 0
            }else{
              D1 <- theta[5]
              C2_coef <- theta[6]
              D2 <- theta[7]
            }
          }else if(splines_degree > 0){

            C1_coef <- theta[3 + 1:nb1]

            if(is.null(init_d1) | is.null(init_d2)){
              D1 <- 0
              C2_coef <- theta[nb1 + 3 + 1:nb2]
              D2 <- 0
            }else{
              D1 <- theta[nb1 + 4]
              C2_coef <- theta[nb1 + 4 + 1:nb2]
              D2 <- theta[nb1 + nb2 + 5]
            }
          }
        }
      }else if(!init_a1_fix){
        if(init_scale_horizontal_fix & init_scale_vertical_fix){
          if(theta[2] < 0){
            return(Inf)
          }

          SCALE_HORIZONTAL <- init_scale_horizontal
          SCALE_VERTICAL <- init_scale_vertical
          A1 <- theta[2]
          B1 <- theta[3]

          if(splines_degree == 0){
            C1_coef <- theta[4]

            if(is.null(init_d1) | is.null(init_d2)){
              D1 <- 0
              A2 <- theta[5]
              B2 <- theta[6]
              C2_coef <- theta[7]
              D2 <- 0
            }else{
              D1 <- theta[5]
              A2 <- theta[6]
              B2 <- theta[7]
              C2_coef <- theta[8]
              D2 <- theta[9]
            }
          }else if(splines_degree > 0){
            C1_coef <- theta[3 + 1:nb1]

            if(is.null(init_d1) | is.null(init_d2)){
              D1 <- 0
              A2 <- theta[nb1 + 4]
              B2 <- theta[nb1 + 5]
              C2_coef <- theta[nb1 + 5 + 1:nb2]
              D2 <- 0
            }else{
              D1 <- theta[nb1 + 4]
              A2 <- theta[nb1 + 5]
              B2 <- theta[nb1 + 6]
              C2_coef <- theta[nb1 + 6 + 1:nb2]
              D2 <- theta[nb1 + nb2 + 7]
            }
          }
        }else if(!init_scale_horizontal_fix & !init_scale_vertical_fix){

          if(theta[4] < 0){
            return(Inf)
          }

          SCALE_HORIZONTAL <- exp(theta[2])
          SCALE_VERTICAL <- exp(theta[3])
          A1 <- theta[4]
          B1 <- theta[5]

          if(splines_degree == 0){
            C1_coef <- theta[6]

            if(is.null(init_d1) | is.null(init_d2)){
              D1 <- 0
              A2 <- theta[7]
              B2 <- theta[8]
              C2_coef <- theta[9]
              D2 <- 0
            }else{
              D1 <- theta[7]
              A2 <- theta[8]
              B2 <- theta[9]
              C2_coef <- theta[10]
              D2 <- theta[11]
            }
          }else if(splines_degree > 0){
            C1_coef <- theta[5 + 1:nb1]

            if(is.null(init_d1) | is.null(init_d2)){
              D1 <- 0
              A2 <- theta[nb1 + 6]
              B2 <- theta[nb1 + 7]
              C2_coef <- theta[nb1 + 7 + 1:nb2]
              D2 <- 0
            }else{
              D1 <- theta[nb1 + 6]
              A2 <- theta[nb1 + 7]
              B2 <- theta[nb1 + 8]
              C2_coef <- theta[nb1 + 8 + 1:nb2]
              D2 <- theta[nb1 + nb2 + 9]
            }
          }
        }
      }
    }else if(init_beta_fix){

      BETA <- init_beta

      if(init_a1_fix){

        A1 <- init_a1
        B1 <- init_b1
        A2 <- init_a2
        B2 <- init_b2

        if(init_scale_horizontal_fix & init_scale_vertical_fix){

          SCALE_HORIZONTAL <- init_scale_horizontal
          SCALE_VERTICAL <- init_scale_vertical

          if(splines_degree == 0){
            C1_coef <- theta[1]
            if(is.null(init_d1) | is.null(init_d2)){
              D1 <- 0
              C2_coef <- theta[2]
              D2 <- 0
            }else{
              D1 <- theta[2]
              C2_coef <- theta[3]
              D2 <- theta[4]
            }
          }else if(splines_degree > 0){
            C1_coef <- theta[1:nb1]

            if(is.null(init_d1) | is.null(init_d2)){
              D1 <- 0
              C2_coef <- theta[nb1 + 1:nb2]
              D2 <- 0
            }else{
              D1 <- theta[nb1 + 1]
              C2_coef <- theta[nb1 + 1 + 1:nb2]
              D2 <- theta[nb1 + nb2 + 2]
            }
          }
        }else if(!init_scale_horizontal_fix & !init_scale_vertical_fix){

          SCALE_HORIZONTAL <- exp(theta[1])
          SCALE_VERTICAL <- exp(theta[2])

          if(splines_degree == 0){

            C1_coef <- theta[3]

            if(is.null(init_d1) | is.null(init_d2)){
              D1 <- 0
              C2_coef <- theta[4]
              D2 <- 0
            }else{
              D1 <- theta[4]
              C2_coef <- theta[5]
              D2 <- theta[6]
            }
          }else if(splines_degree > 0){

            C1_coef <- theta[2 + 1:nb1]

            if(is.null(init_d1) | is.null(init_d2)){
              D1 <- 0
              C2_coef <- theta[nb1 + 2 + 1:nb2]
              D2 <- 0
            }else{
              D1 <- theta[nb1 + 3]
              C2_coef <- theta[nb1 + 3 + 1:nb2]
              D2 <- theta[nb1 + nb2 + 4]
            }
          }
        }
      }else if(!init_a1_fix){
        if(init_scale_horizontal_fix & init_scale_vertical_fix){
          if(theta[1] < 0){
            return(Inf)
          }

          SCALE_HORIZONTAL <- init_scale_horizontal
          SCALE_VERTICAL <- init_scale_vertical
          A1 <- theta[1]
          B1 <- theta[2]

          if(splines_degree == 0){
            C1_coef <- theta[3]

            if(is.null(init_d1) | is.null(init_d2)){
              D1 <- 0
              A2 <- theta[4]
              B2 <- theta[5]
              C2_coef <- theta[6]
              D2 <- 0
            }else{
              D1 <- theta[4]
              A2 <- theta[5]
              B2 <- theta[6]
              C2_coef <- theta[7]
              D2 <- theta[8]
            }
          }else if(splines_degree > 0){
            C1_coef <- theta[2 + 1:nb1]

            if(is.null(init_d1) | is.null(init_d2)){
              D1 <- 0
              A2 <- theta[nb1 + 3]
              B2 <- theta[nb1 + 4]
              C2_coef <- theta[nb1 + 4 + 1:nb2]
              D2 <- 0
            }else{
              D1 <- theta[nb1 + 3]
              A2 <- theta[nb1 + 4]
              B2 <- theta[nb1 + 5]
              C2_coef <- theta[nb1 + 5 + 1:nb2]
              D2 <- theta[nb1 + nb2 + 6]
            }
          }
        }else if(!init_scale_horizontal_fix & !init_scale_vertical_fix){

          if(theta[3] < 0){
            return(Inf)
          }

          SCALE_HORIZONTAL <- exp(theta[1])
          SCALE_VERTICAL <- exp(theta[2])
          A1 <- theta[3]
          B1 <- theta[4]

          if(splines_degree == 0){
            C1_coef <- theta[5]

            if(is.null(init_d1) | is.null(init_d2)){
              D1 <- 0
              A2 <- theta[6]
              B2 <- theta[7]
              C2_coef <- theta[8]
              D2 <- 0
            }else{
              D1 <- theta[6]
              A2 <- theta[7]
              B2 <- theta[8]
              C2_coef <- theta[9]
              D2 <- theta[10]
            }
          }else if(splines_degree > 0){
            C1_coef <- theta[4 + 1:nb1]

            if(is.null(init_d1) | is.null(init_d2)){
              D1 <- 0
              A2 <- theta[nb1 + 5]
              B2 <- theta[nb1 + 6]
              C2_coef <- theta[nb1 + 6 + 1:nb2]
              D2 <- 0
            }else{
              D1 <- theta[nb1 + 5]
              A2 <- theta[nb1 + 6]
              B2 <- theta[nb1 + 7]
              C2_coef <- theta[nb1 + 7 + 1:nb2]
              D2 <- theta[nb1 + nb2 + 8]
            }
          }
        }
      }
    }

    if(splines_degree == 0){
      C1 <- C1_coef
      C2 <- C2_coef
    }else if(splines_degree > 0){
      C1 <- basis1 %*% matrix(C1_coef, ncol = 1)
      C2 <- basis2 %*% matrix(C2_coef, ncol = 1)
    }

    cov_mat <- cov_bi_differential(location = location, beta = BETA,
                                   scale_horizontal = SCALE_HORIZONTAL, scale_vertical = SCALE_VERTICAL,
                                   a1 = A1, b1 = B1, c1 = C1, d1 = D1, a2 = A2, b2 = B2, c2 = C2, d2 = D2,
                                   radius = radius)

    variance1 <- diag(cov_mat[1:nrow(location), 1:nrow(location)])
    variance2 <- diag(cov_mat[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)])

    print(range(cov_mat))

    cat(c("variance1: ", round(max(variance1), 4), " -- variance2: ", round(max(variance2), 4)), '\n')

    cholmat <- tryCatch(chol(cov_mat), error = function(a) numeric(0))
    if(length(cholmat) == 0){
      print("CHOL PROBLEM")
      #Beig <- eigen(cov_mat)
      #g <- Beig$vectors
      #l <- Beig$values
      #ind2 <- which(l <= 1e-12)
      #l[ind2] <- 1e-12

      #l_dd_full <- diag(l)
      #cov_mat <- g %*% l_dd_full %*% t(g)
      #cholmat <- chol(cov_mat)
      return(Inf)
    }

    if(max(variance1) > 2 | max(variance2) > 0.1){
      return(Inf)
    }else{

      ref_lat <- c(40, 0, -40, 40, 0, -40)
      ref_long <- c(-175, -175, -175, -30, -30, -30)

      loc3d_eval <- cbind(ref_long[1], ref_lat[1], seq(0, max(location[, 3]), length.out = 100))

      new_basis1 <- bsplineBasis(loc3d_eval[, 3], splines_degree, knots1)
      new_basis2 <- bsplineBasis(loc3d_eval[, 3], splines_degree, knots2)

      plot_bi_differential(location = loc3d_eval, est_beta = BETA,
                           est_scale_horizontal = SCALE_HORIZONTAL, est_scale_vertical = SCALE_VERTICAL,
                           est_a1 = A1, est_b1 = B1, est_c1 = C1_coef, est_d1 = D1, est_a2 = A2, est_b2 = B2, est_c2 = C2_coef, est_d2 = D2,
                           basis1 = new_basis1, basis2 = new_basis2, radius = radius, splines_degree = splines_degree)

      Sigma_inv <- solve(cov_mat)
      Sigma_log_det <- 2 * sum(log(diag(t(cholmat))))

      out1 <- as.numeric(as.matrix(t(residuals) %*% Sigma_inv %*% residuals))
      out2 <- .5 * nrow(location) * k * log(2 * pi) + .5 * nrow(location) * Sigma_log_det
      out <- out1 + out2

      #out <- -sum(dmvnorm(residuals, mean = rep(0, ncol(cov_mat)), sigma = cov_mat, log = T))

      cat(c("negloglik: ", round(out, 4)), '\n')

      return(out)
    }
  }

  if(!init_beta_fix){
    if(init_a1_fix){
      if(init_scale_horizontal_fix & init_scale_vertical_fix){
        if(is.null(init_d1) | is.null(init_d2)){
          fit <- optim(par = c(init_beta, init_c1, init_c2), fn = NEGLOGLIK, residuals = residuals, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = MAXIT), hessian = T)
        }else{
          fit <- optim(par = c(init_beta, init_c1, init_d1, init_c2, init_d2), fn = NEGLOGLIK, residuals = residuals, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = MAXIT), hessian = T)
        }
      }else{
        if(is.null(init_d1) | is.null(init_d2)){
          fit <- optim(par = c(init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1, init_a2, init_b2, init_c2), fn = NEGLOGLIK, residuals = residuals, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = MAXIT), hessian = T)
        }else{
          fit <- optim(par = c(init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1, init_d1, init_a2, init_b2, init_c2, init_d2), fn = NEGLOGLIK, residuals = residuals, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = MAXIT), hessian = T)
        }
      }
    }else{
      if(init_scale_horizontal_fix & init_scale_vertical_fix){
        if(is.null(init_d1) | is.null(init_d2)){
          fit <- optim(par = c(init_beta, init_a1, init_b1, init_c1, init_a2, init_b2, init_c2), fn = NEGLOGLIK, residuals = residuals, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = MAXIT), hessian = T)
        }else{
          fit <- optim(par = c(init_beta, init_a1, init_b1, init_c1, init_d1, init_a2, init_b2, init_c2, init_d2), fn = NEGLOGLIK, residuals = residuals, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = MAXIT), hessian = T)
        }
      }else{
        if(is.null(init_d1) | is.null(init_d2)){
          print("YOU ARE HERE")
          #fit <- optim(par = c(init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1, init_a2, init_b2, init_c2), fn = NEGLOGLIK, residuals = residuals, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = MAXIT), hessian = T)
          fit <- nlm(NEGLOGLIK, c(init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1, init_a2, init_b2, init_c2), residuals = residuals, basis1 = basis1, basis2 = basis2, hessian = T,
              print.level = 2, iterlim = MAXIT, stepmax = STEPMAX)
        }else{
          fit <- optim(par = c(init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1, init_d1, init_a2, init_b2, init_c2, init_d2), fn = NEGLOGLIK, residuals = residuals, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = MAXIT), hessian = T)
        }
      }
    }
  }else{
    if(init_a1_fix){
      if(init_scale_horizontal_fix & init_scale_vertical_fix){
        if(is.null(init_d1) | is.null(init_d2)){
          fit <- optim(par = c(init_c1, init_c2), fn = NEGLOGLIK, residuals = residuals, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = MAXIT), hessian = T)
        }else{
          fit <- optim(par = c(init_c1, init_d1, init_c2, init_d2), fn = NEGLOGLIK, residuals = residuals, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = MAXIT), hessian = T)
        }
      }else{
        if(is.null(init_d1) | is.null(init_d2)){
          fit <- optim(par = c(init_scale_horizontal, init_scale_vertical, init_c1, init_c2), fn = NEGLOGLIK, residuals = residuals, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = MAXIT), hessian = T)
        }else{
          fit <- optim(par = c(init_scale_horizontal, init_scale_vertical, init_c1, init_d1, init_c2, init_d2), fn = NEGLOGLIK, residuals = residuals, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = MAXIT), hessian = T)
        }
      }
    }else{
      if(init_scale_horizontal_fix & init_scale_vertical_fix){
        if(is.null(init_d1) | is.null(init_d2)){
          fit <- optim(par = c(init_a1, init_b1, init_c1, init_a2, init_b2, init_c2), fn = NEGLOGLIK, residuals = residuals, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = MAXIT), hessian = T)
        }else{
          fit <- optim(par = c(init_a1, init_b1, init_c1, init_d1, init_a2, init_b2, init_c2, init_d2), fn = NEGLOGLIK, residuals = residuals, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = MAXIT), hessian = T)
        }
      }else{
        if(is.null(init_d1) | is.null(init_d2)){
          fit <- optim(par = c(init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1, init_a2, init_b2, init_c2), fn = NEGLOGLIK, residuals = residuals, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = MAXIT), hessian = T)
        }else{
          fit <- optim(par = c(init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1, init_d1, init_a2, init_b2, init_c2, init_d2), fn = NEGLOGLIK, residuals = residuals, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = MAXIT), hessian = T)
        }
      }
    }
  }

  #if(RERUNS > 0){
    #for(aa in 1:RERUNS){
      #fit <- optim(par = fit$par, fn = NEGLOGLIK, residuals = residuals, basis1 = basis1, basis2 = basis2, control = list(trace = 5, maxit = MAXIT), hessian = T)
    #}
  #}

  #est_theta <- fit$par
  #est_hessian <- fit$hessian

  return(fit)
}

#' Compute the value of the bspline bases
#' @param x A matrix of basis function values.
#' @param degree A number indicating the degrees.
#' @param innerknots A vector of knot locations.
#' @param lowknot The first knot location.
#' @param highknot The last knot location.
#' @return A matrix of basis function values.
#' @useDynLib DiffOp, .registration=TRUE
#' @export
bsplineBasis <- function (x, degree, innerknots, lowknot = min(x,innerknots), highknot = max(x,innerknots)) {
  innerknots <- unique (sort (innerknots))
  knots <-
    c(rep(lowknot, degree + 1), innerknots, rep(highknot, degree + 1))
  n <- length (x)
  m <- length (innerknots) + 2 * (degree + 1)
  nf <- length (innerknots) + degree + 1
  basis <- rep (0,  n * nf)
  res <- .C(
    "splinebasis", d = as.integer(degree),
    n = as.integer(n), m = as.integer (m), x = as.double (x), knots = as.double (knots), basis = as.double(basis)
  )
  basis <- matrix (res$basis, n, nf)
  #basis <- basis[,which(colSums(basis) > 0)]
  return (basis)
}

#' Prediction function
#' @param residuals A vector of residuals
#' @param location An nx3 matrix of coordinates.
#' @param location_new An nx3 matrix of coordinates where prediction is required.
#' @param masked_residuals A vector of residuals removed to test prediction performance.
#' @param est_beta A number for colocated correlation parameter.
#' @param est_scale_horizontal A number for the horizontal scale parameter.
#' @param est_scale_vertical A number for the vertical scale parameter.
#' @param est_a1 A number for the anisotropy parameter in Latitude associated with variable 1.
#' @param est_b1 A number for the anisotropy parameter in longitude associated with variable 1.
#' @param est_c1 A number for vector for the nonstationarity parameter in depth associated with variable 1.
#' @param est_d1 A number for the variance parameter of the fully isotropic associated with variable 1.
#' @param est_a2 A number for the anisotropy parameter in Latitude associated with variable 2.
#' @param est_b2 A number for the anisotropy parameter in longitude associated with variable 2.
#' @param est_c2 A number for vector for the nonstationarity parameter in depth associated with variable 2.
#' @param est_d2 A number for the variance parameter of the fully isotropic associated with variable 2.
#' @param radius A number for the radius of the sphere.
#' @param splines_degree A number indicating the degree of the splines.
#' @param knots1 A vector of knot locations for variable 1.
#' @param knots2 A vector of knot locations for variable 2.
#' @import stats
#' @import mvtnorm
#' @return A vector of prediction values.
#' @export
predict_bi_differential <- function(residuals, location, location_new, masked_residuals = NULL, est_beta, est_scale_horizontal, est_scale_vertical, est_a1, est_b1, est_c1, est_d1 = NULL, est_a2, est_b2, est_c2, est_d2 = NULL, radius, splines_degree = 4, knots1, knots2){

  BETA <- est_beta
  SCALE_HORIZONTAL <- est_scale_horizontal
  SCALE_VERTICAL <- est_scale_vertical
  A1 <- est_a1
  B1 <- est_b1
  A2 <- est_a2
  B2 <- est_b2
  C1_coef <- est_c1
  C2_coef <- est_c2

  if(is.null(est_d1) | is.null(est_d2)){
    D1 <- 0
    D2 <- 0
  }else{
    D1 <- est_d1
    D2 <- est_d2
  }

  location_full <- rbind(location, location_new)

  basis1 <- bsplineBasis(location_full[, 3], splines_degree, knots1)
  nb1 <- ncol(basis1)
  basis2 <- bsplineBasis(location_full[, 3], splines_degree, knots2)
  nb2 <- ncol(basis2)

  C1 <- basis1 %*% matrix(C1_coef, ncol = 1)
  C2 <- basis2 %*% matrix(C2_coef, ncol = 1)

  cov_mat <- cov_bi_differential(location = location_full, beta = BETA,
                                 scale_horizontal = SCALE_HORIZONTAL, scale_vertical = SCALE_VERTICAL,
                                 a1 = A1, b1 = B1, c1 = C1, d1 = D1, a2 = A2, b2 = B2, c2 = C2, d2 = D2,
                                 radius = radius)

  if(!is.null(masked_residuals)){

    INDEX_OUT = c(nrow(location) + 1:nrow(location_new), nrow(location) + nrow(location_new) + nrow(location) + 1:nrow(location_new))

    Sigma_insample_dd <- cov_mat[-INDEX_OUT, -INDEX_OUT]
    Sigma_outsample_dd <- cov_mat[INDEX_OUT, INDEX_OUT]
    Sigma_inoutsample_dd <- cov_mat[-INDEX_OUT, INDEX_OUT]

    Sigma_inv <- solve(Sigma_insample_dd)

    print(dim(Sigma_inoutsample_dd))
    print(dim(Sigma_inv))
    print(length(residuals))

    pred_dd <- t(Sigma_inoutsample_dd) %*% Sigma_inv %*% matrix(residuals, ncol = 1)

    PREDICTIONS <- cbind(pred_dd, masked_residuals)
    error <- pred_dd - masked_residuals
    MSE <- as.numeric(as.matrix(t(error) %*% error)) / length(INDEX_OUT)

    print(MSE)
  }
}


#' Plotting the bivariate differential operator cross-covariance function
#'
#' @param location An nx3 matrix of coordinates with latitude and longitude of the reference location and the pressure coordinate is ordered from surface to bottom.
#' @param est_beta A number for colocated correlation parameter.
#' @param est_scale_horizontal A number for the horizontal scale parameter.
#' @param est_scale_vertical A number for the vertical scale parameter.
#' @param est_a1 A number for the anisotropy parameter in Latitude associated with variable 1.
#' @param est_b1 A number for the anisotropy parameter in longitude associated with variable 1.
#' @param est_c1 A number for vector for the nonstationarity parameter in depth associated with variable 1.
#' @param est_d1 A number for the variance parameter of the fully isotropic associated with variable 1.
#' @param est_a2 A number for the anisotropy parameter in Latitude associated with variable 2.
#' @param est_b2 A number for the anisotropy parameter in longitude associated with variable 2.
#' @param est_c2 A number for vector for the nonstationarity parameter in depth associated with variable 2.
#' @param est_d2 A number for the variance parameter of the fully isotropic associated with variable 2.
#' @param radius A number for the radius of the sphere.
#' @param basis1 A matrix of basis function values for variable 1.
#' @param nb1 A number indicating the number of bases for variable 1.
#' @param basis2 A matrix of basis function values for variable 2.
#' @param nb2 A number indicating the number of bases for variable 2.
#' @param splines_degree A number indicating the degree of the splines.
#' @return Figure/plots
#' @export
#' @importFrom graphics par
plot_bi_differential <- function(location, est_beta, est_scale_horizontal, est_scale_vertical, est_a1, est_b1, est_c1, est_d1 = NULL, est_a2, est_b2, est_c2, est_d2 = NULL, radius, basis1, nb1 = ncol(basis1), basis2, nb2 = ncol(basis2), splines_degree = 4){

  BETA <- est_beta
  SCALE_HORIZONTAL <- est_scale_horizontal
  SCALE_VERTICAL <- est_scale_vertical
  A1 <- est_a1
  B1 <- est_b1
  A2 <- est_a2
  B2 <- est_b2

  if(splines_degree == 0){
    C1 <- est_c1
    C2 <- est_c2
  }else if(splines_degree > 0){
    C1 <- basis1 %*% matrix(est_c1, ncol = 1)
    C2 <- basis2 %*% matrix(est_c2, ncol = 1)
  }

  if(is.null(est_d1) | is.null(est_d2)){
    D1 <- 0
    D2 <- 0
  }else{
    D1 <- est_d1
    D2 <- est_d2
  }

  cov_mat <- cov_bi_differential(location = location, beta = BETA,
                                 scale_horizontal = SCALE_HORIZONTAL, scale_vertical = SCALE_VERTICAL,
                                 a1 = A1, b1 = B1, c1 = C1, d1 = D1, a2 = A2, b2 = B2, c2 = C2, d2 = D2,
                                 radius = radius)

  par(mfrow = c(1, 3))

  variance1 <- diag(cov_mat[1:nrow(location), 1:nrow(location)])
  variance2 <- diag(cov_mat[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)])
  covariance12 <- diag(cov_mat[1:nrow(location), nrow(location) + 1:nrow(location)])

  plot(variance1, location[, 3], type = 'l', xlab = 'Variance', ylab = 'Depth', ylim = c(max(location[, 3]), 0))
  plot(variance2, location[, 3], type = 'l', xlab = 'Variance', ylab = 'Depth', ylim = c(max(location[, 3]), 0))
  plot(covariance12 / sqrt(variance1 * variance2), location[, 3], type = 'l', xlab = 'Colocated Correlation', ylab = 'Depth', ylim = c(max(location[, 3]), 0))

}
