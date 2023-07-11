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

#' Compute the kernel smooth empirical cross-covariance values
#'
#' @description
#' \code{compute_emp_cov} computes the kernel smoothed empirical cross-covariance,
#' \eqn{\hat{C}_{ij}(\mathbf{s}_1, \mathbf{s}_2)} between any two locations
#' \eqn{\mathbf{s}_1} and \eqn{\mathbf{s}_2}, which is required by the
#' weighted least squares estimation routine.
#' The formula used to compute \eqn{\hat{C}_{ij}(\mathbf{s}_1, \mathbf{s}_2)}
#' was proposed by Kleiber, W., & Nychka, D. (2012) and has the following form:
#' \deqn{\hat{C}_{ij}(\mathbf{s}_1, \mathbf{s}_2) = \frac{\sum_{\tilde{\mathbf{s}} \in \tilde{\mathcal{D}}} K_{\lambda_{h}, \lambda_{v}}(\mathbf{s}_1, \tilde{\mathbf{s}})^{1/2} K_{\lambda_{h}, \lambda_{v}}(\mathbf{s}_2, \tilde{\mathbf{s}})^{1/2} Z_i (\tilde{\mathbf{s}}) Z_j (\tilde{\mathbf{s}})}{ \{\sum_{\tilde{\mathbf{s}} \in \tilde{\mathcal{D}}} K_{\lambda_{h}, \lambda_{v}}(\mathbf{s}_1, \tilde{\mathbf{s}})\}^{1/2} \{\sum_{\tilde{\mathbf{s}} \in \tilde{\mathcal{D}}} K_{\lambda_{h}, \lambda_{v}}(\mathbf{s}_2, \tilde{\mathbf{s}})\}^{1/2}},}
#' where \eqn{K_{\lambda_{h}, \lambda_{v}}} is a nonnegative kernel function
#' with horizontal and vertical bandwidths, \eqn{\lambda_{h}} and
#' \eqn{\lambda_{v}}, respectively, which control the degree of smoothing.
#'
#' @usage compute_emp_cov(location, variable1_residuals, variable2_residuals,
#' bandwidth_horizontal, bandwidth_vertical, radius)
#'
#' @param location An \eqn{n \times 3} matrix of coordinates.
#' @param variable1_residuals A numeric vector of variable 1 residuals.
#' @param variable2_residuals A numeric vector of variable 2 residuals.
#' @param bandwidth_horizontal A numeric constant parameter for the bandwidth
#' in the horizontal direction.
#' @param bandwidth_vertical A numeric constant parameter for the bandwidth
#' in the vertical direction.
#' @param radius A numeric constant indicating the radius of the sphere.
#'
#' @return A matrix of dimension \eqn{2 n \times 2 n}.
#'
#' @author Mary Lai Salvana \email{yourlainess@gmail.com}
#'
#' @references Kleiber, W., & Nychka, D. (2012). Nonstationary modeling for multivariate spatial processes. \emph{Journal of Multivariate Analysis, 112}, 76-91.
#'
#' @examples
#'
#' data(argo_ref_loc1)
#'
#' loc3d <- cbind(argo_ref_loc1$Longitude, argo_ref_loc1$Latitude, argo_ref_loc1$Pressure)
#'
#' earthRadiusKm = 6371
#'
#' \dontrun{emp_cov <- compute_emp_cov(location = loc3d,
#' variable1_residuals = argo_ref_loc1$TemperatureResiduals,
#' variable2_residuals = argo_ref_loc1$SalinityResiduals,
#' bandwidth_horizontal = 0.009, bandwidth_vertical = 0.03,
#' radius = earthRadiusKm)}
#'
#'
#' @export
compute_emp_cov <- function(location, variable1_residuals, variable2_residuals, bandwidth_horizontal = 0.009, bandwidth_vertical = 0.03, radius){
  LAT1D <- matrix(location[, 2], nrow(location), nrow(location), byrow = F)
  LON1D <- matrix(location[, 1], nrow(location), nrow(location), byrow = F)
  PRES1 <- matrix(location[, 3], nrow(location), nrow(location), byrow = F)
  LAT2D <- matrix(location[, 2], nrow(location), nrow(location), byrow = T)
  LON2D <- matrix(location[, 1], nrow(location), nrow(location), byrow = T)
  PRES2 <- matrix(location[, 3], nrow(location), nrow(location), byrow = T)

  dist0 <- sqrt(h_new(bandwidth_horizontal, bandwidth_vertical, LAT1D, LON1D, PRES1, LAT2D, LON2D, PRES2, radius))
  kernel <- exp(-dist0)
  kernel_sum <- rowSums(kernel)

  denominator <- outer(sqrt(kernel_sum), sqrt(kernel_sum), '*')

  variable1_residuals_matrix <- matrix(variable1_residuals, length(variable1_residuals), length(variable1_residuals), byrow = T)

  emp_covariance1_temp <- sqrt(kernel) * variable1_residuals_matrix
  emp_covariance1 <- emp_covariance1_temp %*% t(emp_covariance1_temp) / denominator

  variable2_residuals_matrix <- matrix(variable2_residuals, length(variable2_residuals), length(variable2_residuals), byrow = T)

  emp_covariance2_temp <- sqrt(kernel) * variable2_residuals_matrix
  emp_covariance2 <- emp_covariance2_temp %*% t(emp_covariance2_temp) / denominator

  emp_crosscovariance <- emp_covariance1_temp %*% t(emp_covariance2_temp) / denominator

  emp_covariance <- rbind(cbind(emp_covariance1, emp_crosscovariance), cbind(t(emp_crosscovariance), emp_covariance2))

  return(emp_covariance)
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
#' a1, b1, c1, d1, a2, b2, c2, d2, radius, splines_degree, knots1, knots2,
#' c1_coef, c2_coef)
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
#' @param splines_degree A number indicating the degree of the splines when
#' using splines to characterize the nonstationary parameters c1 and c2.
#' @param knots1 A vector of knot locations for variable 1 when using splines to
#' characterize c1.
#' @param knots2 A vector of knot locations for variable 2 when using splines to
#' characterize c2.
#' @param c1_coef A numeric vector indicating the splines coefficients for the
#' nonstationary with depth parameter c1 associated with variable 1.
#' @param c2_coef A numeric vector indicating the splines coefficients for the
#' nonstationary with depth parameter c2 associated with variable 2.
#'
#' @useDynLib DiffOp, .registration=TRUE
#'
#' @return A matrix of dimension \eqn{2 n \times 2 n}.
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
#' C1 = sin((loc3d[, 3] + 0.1) * pi / 0.5)
#' C2 = cos((loc3d[, 3] + 0.1) * pi / 0.5)
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
cov_bi_differential <- function(location, beta, scale_horizontal, scale_vertical, a1, b1, c1 = NULL, d1, a2, b2, c2 = NULL, d2, radius, splines_degree = NULL, knots1 = NULL, knots2 = NULL, c1_coef = NULL, c2_coef = NULL){

  LAT1D <- matrix(location[, 2], nrow(location), nrow(location), byrow = F)
  LON1D <- matrix(location[, 1], nrow(location), nrow(location), byrow = F)
  PRES1 <- matrix(location[, 3], nrow(location), nrow(location), byrow = F)
  LAT2D <- matrix(location[, 2], nrow(location), nrow(location), byrow = T)
  LON2D <- matrix(location[, 1], nrow(location), nrow(location), byrow = T)
  PRES2 <- matrix(location[, 3], nrow(location), nrow(location), byrow = T)

  if(!is.null(splines_degree)){
    basis1 <- bsplineBasis(location[, 3], splines_degree, knots1)
    nb1 <- ncol(basis1)
    basis2 <- bsplineBasis(location[, 3], splines_degree, knots2)
    nb2 <- ncol(basis2)

    if(splines_degree == 0){
      c1 <- c1_coef
      c2 <- c2_coef
    }else if(splines_degree > 0){
      c1 <- basis1 %*% matrix(c1_coef, ncol = 1)
      c2 <- basis2 %*% matrix(c2_coef, ncol = 1)
    }
  }

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

#' Weighted Least Squares Estimation of the parameters of the bivariate differential operator cross-covariance function
#'
#' @description
#' \code{est_bi_differential_wls} performs weighted least squares estimation using a
#' Newton-type algorithm to carry out a minimization of the objective function:
#' \deqn{\mathcal{Q}(\boldsymbol{\theta}) = \sum_{i, j = 1}^{q} \sum_{a = 1}^{n} \sum_{b = 1}^{n} w_{ij} \left\{ \hat{C}_{ij}(\mathbf{s}_a, \mathbf{s}_b) - C_{ij}(\mathbf{s}_a, \mathbf{s}_b | \boldsymbol{\theta})\right\}^2,}
#' where \eqn{\boldsymbol{\theta}} is a vector parameters,
#' \eqn{q} is the number of parameters,
#' \eqn{\hat{C}_{ij}(\mathbf{s}_a, \mathbf{s}_b)} is the kernel smoothed
#' empirical cross-covariance computed from the data, and \eqn{w_{ij}} is the
#' weight assigned to component \eqn{(i,j)} of this least squares problem.
#'
#' @usage est_bi_differential_wls(empirical_values, location, init_beta, init_scale_horizontal,
#' init_scale_vertical, init_a1, init_b1, init_c1_coef, init_d1,
#' init_a2, init_b2, init_c2_coef, init_d2,
#' beta_fix, scale_horizontal_fix, scale_vertical_fix,
#' a1_fix, b1_fix, c1_fix, d1_fix, a2_fix, b2_fix, c2_fix, d2_fix,
#' radius, splines_degree, knots1, knots2,
#' w1, w2, w12, iterlim, stepmax, hessian)
#'
#' @param empirical_values A matrix of dimension \eqn{2 n \times 2 n} containing
#' the empirical marginal and cross-covariance values.
#' @param location An \eqn{n \times 3} matrix of coordinates.
#' @param init_beta A numeric constant indicating the initial value for the
#' colocated correlation parameter.
#' @param init_scale_horizontal A numeric constant indicating the initial value for the
#' horizontal scale parameter.
#' @param init_scale_vertical A numeric constant indicating the initial value for the
#' vertical scale parameter.
#' @param init_a1 A numeric constant indicating the initial value for the
#' anisotropy in latitude parameter associated with variable 1.
#' @param init_b1 A numeric constant indicating the initial value for the
#' anisotropy in longitude parameter associated with variable 1.
#' @param init_c1_coef A numeric vector indicating the initial value for the
#' splines coefficients for the nonstationary with depth parameter c1
#' associated with variable 1.
#' @param init_d1 A numeric constant indicating the initial value for the
#' variance parameter from the fully isotropic component associated with variable 1.
#' @param init_a2 A numeric constant indicating the initial value for the
#' anisotropy in latitude parameter associated with variable 2.
#' @param init_b2 A numeric constant indicating the initial value for the
#' anisotropy in longitude parameter associated with variable 2.
#' @param init_c2_coef A numeric vector indicating the initial value for the
#' splines coefficients for the nonstationary with depth parameter c2
#' associated with variable 12.
#' @param init_d2 A numeric constant indicating the initial value for the
#' variance parameter from the fully isotropic component associated with variable 2.
#' @param beta_fix An indicator that takes in the values \code{TRUE} or \code{FALSE} whether the
#' colocated correlation parameter should not be estimated.
#' @param scale_horizontal_fix An indicator that takes in the values \code{TRUE} or \code{FALSE} whether the
#' horizontal scale parameter should not be estimated.
#' @param scale_vertical_fix An indicator that takes in the values \code{TRUE} or \code{FALSE} whether the
#' vertical scale parameter should not be estimated.
#' @param a1_fix If \code{TRUE}, the a1 parameter is not be estimated.
#' @param b1_fix If \code{TRUE}, the b1 parameter is not be estimated.
#' @param c1_fix If \code{TRUE}, the c1 parameters are not be estimated.
#' @param d1_fix If \code{TRUE}, the d1 parameter is not be estimated.
#' @param a2_fix If \code{TRUE}, the a2 parameter is not be estimated.
#' @param b2_fix If \code{TRUE}, the b2 parameter is not be estimated.
#' @param c2_fix If \code{TRUE}, the c2 parameters are not be estimated.
#' @param d2_fix If \code{TRUE}, the d2 parameter is not be estimated.
#' @param radius A numeric constant indicating the radius of the sphere.
#' @param splines_degree A number indicating the degree of the splines.
#' @param knots1 A vector of knot locations for variable 1.
#' @param knots2 A vector of knot locations for variable 2.
#' @param w1 A number indicating the weight for the marginal covariance
#' for variable 1 in the weighted least squares objective function.
#' @param w2 A number indicating the weight for the marginal covariance
#' for variable 2 in the weighted least squares objective function.
#' @param w12 A number indicating the weight for the cross-covariance
#' in the weighted least squares objective function.
#' @param iterlim A number indicating the maximum number of iterations for \code{nlm}.
#' @param stepmax A number indicating the stepmax of nlm.
#' @param hessian If \code{TRUE}, the hessian is not returned.
#'
#' @import stats
#'
#' @useDynLib DiffOp, .registration=TRUE
#'
#' @return A list containing the values of the convergence code,
#' minimum of the objective function,
#' estimated parameters and their corresponding standard deviations,
#' splines degree, and pre-specified knots.
#'
#' @author Mary Lai Salvana \email{yourlainess@gmail.com}
#'
#' @examples
#'
#' data(argo_ref_loc1)
#'
#' loc3d <- cbind(argo_ref_loc1$Longitude, argo_ref_loc1$Latitude, argo_ref_loc1$Pressure)
#'
#' earthRadiusKm = 6371
#'
#' \dontrun{emp_cov <- compute_emp_cov(location = loc3d,
#' variable1_residuals = argo_ref_loc1$TemperatureResiduals,
#' variable2_residuals = argo_ref_loc1$SalinityResiduals,
#' bandwidth_horizontal = 0.009, bandwidth_vertical = 0.03,
#' radius = earthRadiusKm)}
#'
#' KNOTS1 <- c(0.1, 0.5, 0.9)
#' KNOTS2 <- c(0.1, 0.5, 0.9)
#'
#' INIT_BETA = 0.6
#' INIT_SCALE_HORIZONTAL = log(0.02)
#' INIT_SCALE_VERTICAL = log(0.4)
#' INIT_A1 = INIT_A2 = 0.005
#' INIT_B1 = INIT_B2 = 0.005
#' INIT_D1 = INIT_D2 = 0
#'
#' SPLINES_DEGREE = 2
#'
#' set.seed(1235)
#' INIT_C1 <- runif(length(KNOTS1) + SPLINES_DEGREE + 1, -5, 5)
#'
#' set.seed(1236)
#' INIT_C2 <- runif(length(KNOTS2) + SPLINES_DEGREE + 1, -5, 5)
#'
#' \dontrun{est_params_wls <- est_bi_differential_wls(empirical_values = emp_cov, location = loc3d,
#'                                  init_beta = INIT_BETA,
#'                                  init_scale_horizontal = INIT_SCALE_HORIZONTAL,
#'                                  init_scale_vertical = INIT_SCALE_VERTICAL,
#'                                  init_a1 = INIT_A1, init_b1 = INIT_B1,
#'                                  init_c1_coef = INIT_C1, init_d1 = INIT_D1,
#'                                  init_a2 = INIT_A2, init_b2 = INIT_B2,
#'                                  init_c2_coef = INIT_C2, init_d2 = INIT_D2,
#'                                  d1_fix = TRUE, d2_fix = TRUE,
#'                                  radius = earthRadiusKm,
#'                                  splines_degree = SPLINES_DEGREE,
#'                                  knots1 = KNOTS1, knots2 = KNOTS2,
#'                                  w1 = 100, w2 = 50000, w12 = 1000,
#'                                  iterlim = 1, hessian = FALSE)}
#'
#' @export
est_bi_differential_wls <- function(empirical_values, location, init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1_coef, init_d1, init_a2, init_b2, init_c2_coef, init_d2, beta_fix = F, scale_horizontal_fix = F, scale_vertical_fix = F, a1_fix = F, b1_fix = F, c1_fix = F, d1_fix = F, a2_fix = F, b2_fix = F, c2_fix = F, d2_fix = F, radius, splines_degree = 2, knots1, knots2, w1 = 100, w2 = 50000, w12 = 1000, iterlim = 2000, stepmax = 1, hessian = TRUE){

  basis1 <- bsplineBasis(location[, 3], splines_degree, knots1)
  nb1 <- ncol(basis1)
  basis2 <- bsplineBasis(location[, 3], splines_degree, knots2)
  nb2 <- ncol(basis2)

  NEGLOGLIK <- function(theta){

    print(paste("theta = c(", paste(round(theta, 8), collapse=","), ")", sep = ''))

    if(c1_fix & c2_fix){

      BETA <- theta[1]
      SCALE_HORIZONTAL <- exp(theta[2])
      SCALE_VERTICAL <- exp(theta[3])
      A1 <- theta[4] * 1e-3
      B1 <- theta[5] * 1e-3
      C1_coef <- init_c1_coef
      D1 <- init_d1
      A2 <- theta[6] * 1e-3
      B2 <- theta[7] * 1e-3
      C2_coef <- init_c2_coef
      D2 <- init_d2

    }else{
      if(!beta_fix){

        BETA <- 2 * (1 / (1 + exp(-theta[1]))) - 1

        if(a1_fix & b1_fix & a2_fix & b2_fix){

          A1 <- init_a1
          B1 <- init_b1
          A2 <- init_a2
          B2 <- init_b2

          if(scale_horizontal_fix & scale_vertical_fix){

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

              if(d1_fix & d2_fix){
                D1 <- init_d1
                C2_coef <- theta[1 + nb1 + 1:nb2]
                D2 <- init_d2
              }else{
                D1 <- theta[nb1 + 2]
                C2_coef <- theta[nb1 + 2 + 1:nb2]
                D2 <- theta[nb1 + nb2 + 3]
              }
            }
          }else if(!scale_horizontal_fix & !scale_vertical_fix){

            SCALE_HORIZONTAL <- exp(theta[2])
            SCALE_VERTICAL <- exp(theta[3])

            if(splines_degree == 0){

              C1_coef <- theta[4]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                C2_coef <- theta[5]
                D2 <- init_d2
              }else{
                D1 <- theta[5]
                C2_coef <- theta[6]
                D2 <- theta[7]
              }
            }else if(splines_degree > 0){

              C1_coef <- theta[3 + 1:nb1]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                C2_coef <- theta[nb1 + 3 + 1:nb2]
                D2 <- init_d2
              }else{
                D1 <- theta[nb1 + 4]
                C2_coef <- theta[nb1 + 4 + 1:nb2]
                D2 <- theta[nb1 + nb2 + 5]
              }
            }
          }else if(!scale_horizontal_fix & !scale_vertical_fix){

            SCALE_HORIZONTAL <- exp(theta[2])
            SCALE_VERTICAL <- exp(theta[3])

            if(splines_degree == 0){

              C1_coef <- theta[4]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                C2_coef <- theta[5]
                D2 <- init_d2
              }else{
                D1 <- theta[5]
                C2_coef <- theta[6]
                D2 <- theta[7]
              }
            }else if(splines_degree > 0){

              C1_coef <- theta[3 + 1:nb1]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                C2_coef <- theta[nb1 + 3 + 1:nb2]
                D2 <- init_d2
              }else{
                D1 <- theta[nb1 + 4]
                C2_coef <- theta[nb1 + 4 + 1:nb2]
                D2 <- theta[nb1 + nb2 + 5]
              }
            }
          }
        }else if(!a1_fix){
          if(scale_horizontal_fix & scale_vertical_fix){

            SCALE_HORIZONTAL <- init_scale_horizontal
            SCALE_VERTICAL <- init_scale_vertical
            A1 <- theta[2] * 1e-3
            B1 <- theta[3] * 1e-3

            if(splines_degree == 0){
              C1_coef <- theta[4]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                A2 <- theta[5] * 1e-3
                B2 <- theta[6] * 1e-3
                C2_coef <- theta[7]
                D2 <- init_d2
              }else{
                D1 <- theta[5]
                A2 <- theta[6] * 1e-3
                B2 <- theta[7] * 1e-3
                C2_coef <- theta[8]
                D2 <- theta[9]
              }
            }else if(splines_degree > 0){
              C1_coef <- theta[3 + 1:nb1]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                A2 <- theta[3 + nb1 + 1] * 1e-3
                B2 <- theta[3 + nb1 + 2] * 1e-3
                C2_coef <- theta[3 + nb1 + 2 + 1:nb2]
                D2 <- init_d2
              }else{
                D1 <- theta[nb1 + 4]
                A2 <- theta[nb1 + 5] * 1e-3
                B2 <- theta[nb1 + 6] * 1e-3
                C2_coef <- theta[nb1 + 6 + 1:nb2]
                D2 <- theta[nb1 + nb2 + 7]
              }
            }
          }else if(!scale_horizontal_fix & !scale_vertical_fix){

            SCALE_HORIZONTAL <- exp(theta[2])
            SCALE_VERTICAL <- exp(theta[3])
            A1 <- theta[4] * 1e-3
            B1 <- theta[5] * 1e-3

            if(splines_degree == 0){
              C1_coef <- theta[6]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                A2 <- theta[7] * 1e-3
                B2 <- theta[8] * 1e-3
                C2_coef <- theta[9]
                D2 <- init_d2
              }else{
                D1 <- theta[7]
                A2 <- theta[8] * 1e-3
                B2 <- theta[9] * 1e-3
                C2_coef <- theta[10]
                D2 <- theta[11]
              }
            }else if(splines_degree > 0){
              C1_coef <- theta[5 + 1:nb1]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                A2 <- theta[5 + nb1 + 1] * 1e-3
                B2 <- theta[5 + nb1 + 2] * 1e-3
                C2_coef <- theta[5 + nb1 + 2 + 1:nb2]
                D2 <- init_d2
              }else{
                D1 <- theta[nb1 + 6]
                A2 <- theta[nb1 + 7] * 1e-3
                B2 <- theta[nb1 + 8] * 1e-3
                C2_coef <- theta[nb1 + 8 + 1:nb2]
                D2 <- theta[nb1 + nb2 + 9]
              }
            }
          }else if(!scale_horizontal_fix & scale_vertical_fix){

            SCALE_HORIZONTAL <- exp(theta[2])
            SCALE_VERTICAL <- init_scale_vertical
            A1 <- theta[3] * 1e-3
            B1 <- theta[4] * 1e-3

            if(splines_degree == 0){
              C1_coef <- theta[5]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                A2 <- theta[6] * 1e-3
                B2 <- theta[7] * 1e-3
                C2_coef <- theta[8]
                D2 <- init_d2
              }else{
                D1 <- theta[6]
                A2 <- theta[7] * 1e-3
                B2 <- theta[8] * 1e-3
                C2_coef <- theta[9]
                D2 <- theta[10]
              }
            }else if(splines_degree > 0){
              C1_coef <- theta[4 + 1:nb1]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                A2 <- theta[4 + nb1 + 1] * 1e-3
                B2 <- theta[4 + nb1 + 2] * 1e-3
                C2_coef <- theta[4 + nb1 + 2 + 1:nb2]
                D2 <- init_d2
              }else{
                D1 <- theta[4 + nb1 + 1]
                A2 <- theta[4 + nb1 + 2] * 1e-3
                B2 <- theta[4 + nb1 + 3] * 1e-3
                C2_coef <- theta[4 + nb1 + 3 + 1:nb2]
                D2 <- theta[4 + nb1 + 3 + nb2 + 1]
              }
            }
          }
        }
      }else if(beta_fix){

        BETA <- init_beta

        if(a1_fix & b1_fix & a2_fix & b2_fix){

          A1 <- init_a1
          B1 <- init_b1
          A2 <- init_a2
          B2 <- init_b2

          if(scale_horizontal_fix & scale_vertical_fix){

            SCALE_HORIZONTAL <- init_scale_horizontal
            SCALE_VERTICAL <- init_scale_vertical

            if(splines_degree == 0){
              C1_coef <- theta[1]
              if(d1_fix & d2_fix){
                D1 <- init_d1
                C2_coef <- theta[2]
                D2 <- init_d2
              }else{
                D1 <- theta[2]
                C2_coef <- theta[3]
                D2 <- theta[4]
              }
            }else if(splines_degree > 0){
              C1_coef <- theta[1:nb1]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                C2_coef <- theta[nb1 + 1:nb2]
                D2 <- init_d2
              }else{
                D1 <- theta[nb1 + 1]
                C2_coef <- theta[nb1 + 1 + 1:nb2]
                D2 <- theta[nb1 + nb2 + 2]
              }
            }
          }else if(!scale_horizontal_fix & !scale_vertical_fix){

            SCALE_HORIZONTAL <- exp(theta[1])
            SCALE_VERTICAL <- exp(theta[2])

            if(splines_degree == 0){

              C1_coef <- theta[3]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                C2_coef <- theta[4]
                D2 <- init_d2
              }else{
                D1 <- theta[4]
                C2_coef <- theta[5]
                D2 <- theta[6]
              }
            }else if(splines_degree > 0){

              C1_coef <- theta[2 + 1:nb1]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                C2_coef <- theta[nb1 + 2 + 1:nb2]
                D2 <- init_d2
              }else{
                D1 <- theta[nb1 + 3]
                C2_coef <- theta[nb1 + 3 + 1:nb2]
                D2 <- theta[nb1 + nb2 + 4]
              }
            }
          }
        }else{
          if(scale_horizontal_fix & scale_vertical_fix){

            SCALE_HORIZONTAL <- init_scale_horizontal
            SCALE_VERTICAL <- init_scale_vertical
            A1 <- theta[1] * 1e-3
            B1 <- theta[2] * 1e-3

            if(splines_degree == 0){
              C1_coef <- theta[3]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                A2 <- theta[4] * 1e-3
                B2 <- theta[5] * 1e-3
                C2_coef <- theta[6]
                D2 <- init_d2
              }else{
                D1 <- theta[4]
                A2 <- theta[5] * 1e-3
                B2 <- theta[6] * 1e-3
                C2_coef <- theta[7]
                D2 <- theta[8]
              }
            }else if(splines_degree > 0){
              C1_coef <- theta[2 + 1:nb1]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                A2 <- theta[2 + nb1 + 1] * 1e-3
                B2 <- theta[2 + nb1 + 2] * 1e-3
                C2_coef <- theta[2 + nb1 + 2 + 1:nb2]
                D2 <- init_d2
              }else{
                D1 <- theta[nb1 + 3]
                A2 <- theta[nb1 + 4] * 1e-3
                B2 <- theta[nb1 + 5] * 1e-3
                C2_coef <- theta[nb1 + 5 + 1:nb2]
                D2 <- theta[nb1 + nb2 + 6]
              }
            }
          }else if(!scale_horizontal_fix & !scale_vertical_fix){

            SCALE_HORIZONTAL <- exp(theta[1])
            SCALE_VERTICAL <- exp(theta[2])
            A1 <- theta[3] * 1e-3
            B1 <- theta[4] * 1e-3

            if(splines_degree == 0){
              C1_coef <- theta[5]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                A2 <- theta[6] * 1e-3
                B2 <- theta[7] * 1e-3
                C2_coef <- theta[8]
                D2 <- init_d2
              }else{
                D1 <- theta[6]
                A2 <- theta[7] * 1e-3
                B2 <- theta[8] * 1e-3
                C2_coef <- theta[9]
                D2 <- theta[10]
              }
            }else if(splines_degree > 0){
              C1_coef <- theta[4 + 1:nb1]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                A2 <- theta[nb1 + 5] * 1e-3
                B2 <- theta[nb1 + 6] * 1e-3
                C2_coef <- theta[nb1 + 6 + 1:nb2]
                D2 <- init_d2
              }else{
                D1 <- theta[nb1 + 5]
                A2 <- theta[nb1 + 6] * 1e-3
                B2 <- theta[nb1 + 7] * 1e-3
                C2_coef <- theta[nb1 + 7 + 1:nb2]
                D2 <- theta[nb1 + nb2 + 8]
              }
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

    #param_list <- list(beta = BETA, scale_horizontal = SCALE_HORIZONTAL, scale_vertical = SCALE_VERTICAL, a1 = A1, b1 = B1, c1 = C1_coef, d1 = D1, a2 = A2, b2 = B2, c2 = C2_coef, d2 = D2)
    #print(param_list)

    cov_mat <- cov_bi_differential(location = location, beta = BETA,
                                   scale_horizontal = SCALE_HORIZONTAL, scale_vertical = SCALE_VERTICAL,
                                   a1 = A1, b1 = B1, c1 = C1, d1 = D1, a2 = A2, b2 = B2, c2 = C2, d2 = D2,
                                   radius = radius)

    covariance1 <- cov_mat[1:nrow(location), 1:nrow(location)]
    covariance2 <- cov_mat[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)]
    covariance12 <- cov_mat[1:nrow(location), nrow(location) + 1:nrow(location)]
    #correlation12 <- diag(covariance12) / sqrt(diag(covariance1) * diag(covariance2))

    emp_covariance1 <- empirical_values[1:nrow(location), 1:nrow(location)]
    emp_covariance2 <- empirical_values[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)]
    emp_covariance12 <- empirical_values[1:nrow(location), nrow(location) + 1:nrow(location)]
    #emp_correlation12 <- diag(emp_covariance12) / sqrt(diag(emp_covariance1) * diag(emp_covariance2))

    #out <- sum((diag(emp_covariance1) - diag(covariance1))^2 * 100 + (diag(emp_covariance2) - diag(covariance2))^2 * 50000 + (emp_correlation12 - correlation12)^2 * 100)
    out <- sum((emp_covariance1 - covariance1)^2 * w1 + (emp_covariance2 - covariance2)^2 * w2 + (emp_covariance12 - covariance12)^2 * w12)

    return(out)
  }

  if(c1_fix & c2_fix){
    fit <- nlm(NEGLOGLIK, c(init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_a2, init_b2), hessian = T,
               print.level = 2, iterlim = iterlim, stepmax = stepmax)

  }else{
    if(!beta_fix){
      if(a1_fix & b1_fix & a2_fix & b2_fix){
        if(scale_horizontal_fix & scale_vertical_fix){
          if(d1_fix & d2_fix){
            fit <- nlm(NEGLOGLIK, c(init_beta, init_c1_coef, init_c2_coef), hessian = hessian,
                       print.level = 2, iterlim = iterlim, stepmax = stepmax)

            BETA <- 2 * (1 / (1 + exp(-fit$estimate[1]))) - 1
            SCALE_HORIZONTAL <- init_scale_horizontal
            SCALE_VERTICAL <- init_scale_vertical
            A1 <- init_a1
            B1 <- init_b1
            C1_coef <- fit$estimate[1 + 1:nb1]
            D1 <- init_d1

            A2 <- init_a2
            B2 <- init_b2
            C2_coef <- fit$estimate[1 + nb1 + 1:nb2]
            D2 <- init_d2
          }else{
            fit <- optim(par = c(init_beta, init_c1_coef, init_d1, init_c2_coef, init_d2), fn = NEGLOGLIK, control = list(trace = 5, iterlim = iterlim), hessian = hessian)
          }
        }else{
          if(d1_fix & d2_fix){
            fit <- optim(par = c(init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1_coef, init_a2, init_b2, init_c2_coef), fn = NEGLOGLIK, control = list(trace = 5, iterlim = iterlim), hessian = hessian)
          }else{
            fit <- optim(par = c(init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1_coef, init_d1, init_a2, init_b2, init_c2_coef, init_d2), fn = NEGLOGLIK, control = list(trace = 5, iterlim = iterlim), hessian = hessian)
          }
        }
      }else{
        if(scale_horizontal_fix & scale_vertical_fix){
          if(d1_fix & d2_fix){
            fit <- nlm(NEGLOGLIK, c(init_beta, init_a1, init_b1, init_c1_coef, init_a2, init_b2, init_c2_coef), hessian = hessian,
                       print.level = 2, iterlim = iterlim, stepmax = stepmax)

            BETA <- 2 * (1 / (1 + exp(-fit$estimate[1]))) - 1
            SCALE_HORIZONTAL <- init_scale_horizontal
            SCALE_VERTICAL <- init_scale_vertical
            A1 <- fit$estimate[2] * 1e-3
            B1 <- fit$estimate[3] * 1e-3
            C1_coef <- fit$estimate[3 + 1:nb1]
            D1 <- init_d1

            A2 <- fit$estimate[3 + nb1 + 1] * 1e-3
            B2 <- fit$estimate[3 + nb1 + 2] * 1e-3
            C2_coef <- fit$estimate[3 + nb1 + 2 + 1:nb2]
            D2 <- init_d2
          }else{
            fit <- optim(par = c(init_beta, init_a1, init_b1, init_c1_coef, init_d1, init_a2, init_b2, init_c2_coef, init_d2), fn = NEGLOGLIK, control = list(trace = 5, iterlim = iterlim), hessian = hessian)
          }
        }else if(!scale_horizontal_fix & !scale_vertical_fix){
          if(d1_fix & d2_fix){
            fit <- nlm(NEGLOGLIK, c(init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1_coef, init_a2, init_b2, init_c2_coef), hessian = hessian,
                       print.level = 2, iterlim = iterlim, stepmax = stepmax)

            BETA <- 2 * (1 / (1 + exp(-fit$estimate[1]))) - 1
            SCALE_HORIZONTAL <- exp(fit$estimate[2])
            SCALE_VERTICAL <- exp(fit$estimate[3])
            A1 <- fit$estimate[4] * 1e-3
            B1 <- fit$estimate[5] * 1e-3
            C1_coef <- fit$estimate[5 + 1:nb1]
            D1 <- init_d1

            A2 <- fit$estimate[5 + nb1 + 1] * 1e-3
            B2 <- fit$estimate[5 + nb1 + 2] * 1e-3
            C2_coef <- fit$estimate[5 + nb1 + 2 + 1:nb2]
            D2 <- init_d2
          }else{
            fit <- optim(par = c(init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1_coef, init_d1, init_a2, init_b2, init_c2_coef, init_d2), fn = NEGLOGLIK, control = list(trace = 5, iterlim = iterlim), hessian = hessian)
          }
        }
      }
    }else{
      if(a1_fix & b1_fix & a2_fix & b2_fix){
        if(scale_horizontal_fix & scale_vertical_fix){
          if(d1_fix & d2_fix){
            fit <- nlm(NEGLOGLIK, c(init_c1_coef, init_c2_coef), hessian = hessian,
                       print.level = 2, iterlim = iterlim, stepmax = stepmax)

            BETA <- init_beta
            SCALE_HORIZONTAL <- init_scale_horizontal
            SCALE_VERTICAL <- init_scale_vertical
            A1 <- init_a1
            B1 <- init_b1
            C1_coef <- fit$estimate[1:nb1]
            D1 <- init_d1

            A2 <- init_a2
            B2 <- init_b2
            C2_coef <- fit$estimate[nb1 + 1:nb2]
            D2 <- init_d2
          }else{
            fit <- optim(par = c(init_c1_coef, init_d1, init_c2_coef, init_d2), fn = NEGLOGLIK, control = list(trace = 5, iterlim = iterlim), hessian = hessian)
          }
        }else{
          if(d1_fix & d2_fix){
            fit <- optim(par = c(init_scale_horizontal, init_scale_vertical, init_c1_coef, init_c2_coef), fn = NEGLOGLIK, control = list(trace = 5, iterlim = iterlim), hessian = hessian)
          }else{
            fit <- optim(par = c(init_scale_horizontal, init_scale_vertical, init_c1_coef, init_d1, init_c2_coef, init_d2), fn = NEGLOGLIK, control = list(trace = 5, iterlim = iterlim), hessian = hessian)
          }
        }
      }else{
        if(scale_horizontal_fix & scale_vertical_fix){
          if(d1_fix & d2_fix){
            fit <- nlm(NEGLOGLIK, c(init_a1, init_b1, init_c1_coef, init_a2, init_b2, init_c2_coef), hessian = hessian,
                       print.level = 2, iterlim = iterlim, stepmax = stepmax)

            BETA <- init_beta
            SCALE_HORIZONTAL <- init_scale_horizontal
            SCALE_VERTICAL <- init_scale_vertical
            A1 <- fit$estimate[1] * 1e-3
            B1 <- fit$estimate[2] * 1e-3
            C1_coef <- fit$estimate[2 + 1:nb1]
            D1 <- init_d1

            A2 <- fit$estimate[2 + nb1 + 1] * 1e-3
            B2 <- fit$estimate[2 + nb1 + 2] * 1e-3
            C2_coef <- fit$estimate[2 + nb1 + 2 + 1:nb2]
            D2 <- init_d2
          }else{
            fit <- optim(par = c(init_a1, init_b1, init_c1_coef, init_d1, init_a2, init_b2, init_c2_coef, init_d2), fn = NEGLOGLIK, control = list(trace = 5, iterlim = iterlim), hessian = hessian)
          }
        }else{
          if(d1_fix & d2_fix){
            fit <- optim(par = c(init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1_coef, init_a2, init_b2, init_c2_coef), fn = NEGLOGLIK, control = list(trace = 5, iterlim = iterlim), hessian = hessian)
          }else{
            fit <- optim(par = c(init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1_coef, init_d1, init_a2, init_b2, init_c2_coef, init_d2), fn = NEGLOGLIK, control = list(trace = 5, iterlim = iterlim), hessian = hessian)
          }
        }
      }
    }
  }
  results <- list(convergence_code = fit$code, iterations = fit$iterations, minimum = fit$minimum, theta = fit$estimate, est_beta = BETA, est_scale_horizontal = SCALE_HORIZONTAL, est_scale_vertical = SCALE_VERTICAL, est_a1 = A1, est_b1 = B1, est_c1_coef = C1_coef, est_d1 = D1, est_a2 = A2, est_b2 = B2, est_c2_coef = C2_coef, est_d2 = D2, splines_degree = splines_degree, knots1 = knots1, knots2 = knots2)
  return(results)
}


#' Maximum Likelihood Estimation of the parameters of the bivariate differential operator cross-covariance function
#'
#' @description
#' \code{est_bi_differential_mle} performs maximum likelihood estimation using a
#' Newton-type algorithm to carry out a minimization of the negative of the
#' Gaussian log-likelihood function:
#' \deqn{l(\boldsymbol{\theta}) = -\frac{q n}{2} \log (2 \pi) - \frac{1}{2} \log |\boldsymbol{\Sigma}(\boldsymbol{\theta})| - \frac{1}{2} \mathbf{Z}^{\top} \boldsymbol{\Sigma}(\boldsymbol{\theta})^{-1} \mathbf{Z},}
#' Here \eqn{\boldsymbol{\theta}} is a vector parameters,
#' \eqn{q} is the number of parameters,
#' \eqn{\boldsymbol{\Sigma}(\boldsymbol{\theta})}
#' is the \eqn{q n \times q n} cross-covariance matrix for \eqn{\mathbf{Z}}
#' with determinant \eqn{|\boldsymbol{\Sigma}(\boldsymbol{\theta})|}, and
#' \eqn{\mathbf{Z} = \{\mathbf{Z}(\mathbf{s}_1)^{\top}, \ldots, \mathbf{Z}(\mathbf{s}_n)^{\top} \}^{\top}}
#' is the vector of residuals such that \eqn{\mathbf{Z}(\mathbf{s}) = \{ Z_1 (\mathbf{s}), Z_2 (\mathbf{s})\}^{\top}}.
#'
#' @usage est_bi_differential_mle(residuals, location, init_beta, init_scale_horizontal,
#' init_scale_vertical, init_a1, init_b1, init_c1_coef, init_d1,
#' init_a2, init_b2, init_c2_coef, init_d2,
#' beta_fix, scale_horizontal_fix, scale_vertical_fix,
#' a1_fix, b1_fix, c1_fix, d1_fix, a2_fix, b2_fix, c2_fix, d2_fix,
#' radius, splines_degree, knots1, knots2, iterlim, stepmax, hessian)
#'
#' @param residuals A \eqn{2 n} vector of residuals with
#' variable 1 and 2 residuals appended next to each other.
#' @param location An \eqn{n \times 3} matrix of coordinates.
#' @param init_beta A numeric constant indicating the initial value for the
#' colocated correlation parameter.
#' @param init_scale_horizontal A numeric constant indicating the initial value for the
#' horizontal scale parameter.
#' @param init_scale_vertical A numeric constant indicating the initial value for the
#' vertical scale parameter.
#' @param init_a1 A numeric constant indicating the initial value for the
#' anisotropy in latitude parameter associated with variable 1.
#' @param init_b1 A numeric constant indicating the initial value for the
#' anisotropy in longitude parameter associated with variable 1.
#' @param init_c1_coef A numeric vector indicating the initial value for the
#' splines coefficients for the nonstationary with depth parameter c1
#' associated with variable 1.
#' @param init_d1 A numeric constant indicating the initial value for the
#' variance parameter from the fully isotropic component associated with variable 1.
#' @param init_a2 A numeric constant indicating the initial value for the
#' anisotropy in latitude parameter associated with variable 2.
#' @param init_b2 A numeric constant indicating the initial value for the
#' anisotropy in longitude parameter associated with variable 2.
#' @param init_c2_coef A numeric vector indicating the initial value for the
#' splines coefficients for the nonstationary with depth parameter c2
#' associated with variable 2.
#' @param init_d2 A numeric constant indicating the initial value for the
#' variance parameter from the fully isotropic component associated with variable 2.
#' @param beta_fix An indicator that takes in the values \code{TRUE} or \code{FALSE} whether the
#' colocated correlation parameter should not be estimated.
#' @param scale_horizontal_fix An indicator that takes in the values \code{TRUE} or \code{FALSE} whether the
#' horizontal scale parameter should not be estimated.
#' @param scale_vertical_fix An indicator that takes in the values \code{TRUE} or \code{FALSE} whether the
#' vertical scale parameter should not be estimated.
#' @param a1_fix If \code{TRUE}, the a1 parameter is not be estimated.
#' @param b1_fix If \code{TRUE}, the b1 parameter is not be estimated.
#' @param c1_fix If \code{TRUE}, the c1 parameters are not be estimated.
#' @param d1_fix If \code{TRUE}, the d1 parameter is not be estimated.
#' @param a2_fix If \code{TRUE}, the a2 parameter is not be estimated.
#' @param b2_fix If \code{TRUE}, the b2 parameter is not be estimated.
#' @param c2_fix If \code{TRUE}, the c2 parameters are not be estimated.
#' @param d2_fix If \code{TRUE}, the d2 parameter is not be estimated.
#' @param radius A numeric constant indicating the radius of the sphere.
#' @param splines_degree A number indicating the degree of the splines.
#' @param knots1 A vector of knot locations for variable 1.
#' @param knots2 A vector of knot locations for variable 2.
#' @param iterlim A number indicating the maximum number of iterations for \code{nlm}.
#' @param stepmax A number indicating the stepmax of nlm.
#' @param hessian If \code{TRUE}, the hessian is not returned.
#'
#' @import stats
#'
#' @useDynLib DiffOp, .registration=TRUE
#'
#' @return A list containing the values of the convergence code,
#' log likelihood function,
#' estimated parameters and their corresponding standard deviations,
#' splines degree, and pre-specified knots.
#'
#' @author Mary Lai Salvana \email{yourlainess@gmail.com}
#'
#' @examples
#'
#' library(dplyr)
#'
#' x <- seq(0, 1, length.out = 5)
#' y <- seq(0, 1, length.out = 5)
#'
#' loc2d <- expand.grid(x, y) %>% as.matrix()
#' depth <- seq(0, 1, length.out = 5)
#' loc3d <- cbind(rep(loc2d[, 1], each = length(depth)), rep(loc2d[, 2], each = length(depth)), depth)
#'
#' earthRadiusKm = 6371
#'
#' BETA = 0.5
#' SCALE_HORIZONTAL = 0.03
#' SCALE_VERTICAL = 0.3
#' A1 = A2 = 0.00001
#' B1 = B2 = 0.00001
#' C1 = sin((loc3d[, 3] + 0.1) * pi / 0.5)
#' C2 = cos((loc3d[, 3] + 0.1) * pi / 0.5)
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
#' library(MASS)
#'
#' set.seed(1235)
#' Z <- mvrnorm(1, mu = rep(0, ncol(cov_mat)), Sigma = cov_mat)
#' Z1 <- Z[1:nrow(loc3d)]
#' Z2 <- Z[nrow(loc3d) + 1:nrow(loc3d)]
#'
#' KNOTS1 <- c(0.1, 0.5, 0.9)
#' KNOTS2 <- c(0.1, 0.5, 0.9)
#'
#' INIT_BETA = 0.6
#' INIT_SCALE_HORIZONTAL = log(0.02)
#' INIT_SCALE_VERTICAL = log(0.4)
#' INIT_A1 = INIT_A2 = 0.005
#' INIT_B1 = INIT_B2 = 0.005
#' INIT_D1 = INIT_D2 = 0
#'
#' SPLINES_DEGREE = 2
#'
#' set.seed(1235)
#' INIT_C1 <- runif(length(KNOTS1) + SPLINES_DEGREE + 1, -5, 5)
#'
#' set.seed(1236)
#' INIT_C2 <- runif(length(KNOTS2) + SPLINES_DEGREE + 1, -5, 5)
#'
#' \dontrun{est_params_mle <- est_bi_differential_mle(residuals = Z, location = loc3d, init_beta = INIT_BETA,
#'                                  init_scale_horizontal = INIT_SCALE_HORIZONTAL,
#'                                  init_scale_vertical = INIT_SCALE_VERTICAL,
#'                                  init_a1 = INIT_A1, init_b1 = INIT_B1,
#'                                  init_c1_coef = INIT_C1, init_d1 = INIT_D1,
#'                                  init_a2 = INIT_A2, init_b2 = INIT_B2,
#'                                  init_c2_coef = INIT_C2, init_d2 = INIT_D2,
#'                                  d1_fix = TRUE, d2_fix = TRUE,
#'                                  radius = earthRadiusKm,
#'                                  splines_degree = SPLINES_DEGREE,
#'                                  knots1 = KNOTS1, knots2 = KNOTS2,
#'                                  iterlim = 1, hessian = T)}
#'
#'
#' @export
est_bi_differential_mle <- function(residuals, location, init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1_coef, init_d1, init_a2, init_b2, init_c2_coef, init_d2, beta_fix = F, scale_horizontal_fix = F, scale_vertical_fix = F, a1_fix = F, b1_fix = F, c1_fix = F, d1_fix = F, a2_fix = F, b2_fix = F, c2_fix = F, d2_fix = F, radius, splines_degree = 2, knots1, knots2, iterlim = 2000, stepmax = 1, hessian = TRUE){

  basis1 <- bsplineBasis(location[, 3], splines_degree, knots1)
  nb1 <- ncol(basis1)
  basis2 <- bsplineBasis(location[, 3], splines_degree, knots2)
  nb2 <- ncol(basis2)

  NEGLOGLIK <- function(theta){

    print(paste("theta = c(", paste(round(theta, 8), collapse=","), ")", sep = ''))

    if(c1_fix & c2_fix){

      BETA <- 1 / (1 + exp(-theta[1]))
      SCALE_HORIZONTAL <- exp(theta[2])
      SCALE_VERTICAL <- exp(theta[3])
      A1 <- exp(theta[4]) * 1e-3
      B1 <- theta[5] * 1e-3
      C1_coef <- init_c1_coef
      D1 <- init_d1
      A2 <- theta[6] * 1e-3
      B2 <- theta[7] * 1e-3
      C2_coef <- init_c2_coef
      D2 <- init_d2

    }else{
      if(!beta_fix){

        BETA <- 1 / (1 + exp(-theta[1]))

        if(a1_fix & b1_fix & a2_fix & b2_fix){

          A1 <- init_a1
          B1 <- init_b1
          A2 <- init_a2
          B2 <- init_b2

          if(scale_horizontal_fix & scale_vertical_fix){

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

              if(d1_fix & d2_fix){
                D1 <- init_d1
                C2_coef <- theta[nb1 + 1 + 1:nb2]
                D2 <- init_d2
              }else{
                D1 <- theta[nb1 + 2]
                C2_coef <- theta[nb1 + 2 + 1:nb2]
                D2 <- theta[nb1 + nb2 + 3]
              }
            }
          }else if(!scale_horizontal_fix & !scale_vertical_fix){

            SCALE_HORIZONTAL <- exp(theta[2])
            SCALE_VERTICAL <- exp(theta[3])

            if(splines_degree == 0){

              C1_coef <- theta[4]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                C2_coef <- theta[5]
                D2 <- init_d2
              }else{
                D1 <- theta[5]
                C2_coef <- theta[6]
                D2 <- theta[7]
              }
            }else if(splines_degree > 0){

              C1_coef <- theta[3 + 1:nb1]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                C2_coef <- theta[nb1 + 3 + 1:nb2]
                D2 <- init_d2
              }else{
                D1 <- theta[nb1 + 4]
                C2_coef <- theta[nb1 + 4 + 1:nb2]
                D2 <- theta[nb1 + nb2 + 5]
              }
            }
          }
        }else{
          if(scale_horizontal_fix & scale_vertical_fix){

            SCALE_HORIZONTAL <- init_scale_horizontal
            SCALE_VERTICAL <- init_scale_vertical
            A1 <- exp(theta[2]) * 1e-3
            B1 <- theta[3] * 1e-3

            if(splines_degree == 0){
              C1_coef <- theta[4]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                A2 <- theta[5] * 1e-3
                B2 <- theta[6] * 1e-3
                C2_coef <- theta[7]
                D2 <- init_d2
              }else{
                D1 <- theta[5]
                A2 <- theta[6] * 1e-3
                B2 <- theta[7] * 1e-3
                C2_coef <- theta[8]
                D2 <- theta[9]
              }
            }else if(splines_degree > 0){
              C1_coef <- theta[3 + 1:nb1]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                A2 <- theta[nb1 + 4] * 1e-3
                B2 <- theta[nb1 + 5] * 1e-3
                C2_coef <- theta[nb1 + 5 + 1:nb2]
                D2 <- init_d2
              }else{
                D1 <- theta[nb1 + 4]
                A2 <- theta[nb1 + 5] * 1e-3
                B2 <- theta[nb1 + 6] * 1e-3
                C2_coef <- theta[nb1 + 6 + 1:nb2]
                D2 <- theta[nb1 + nb2 + 7]
              }
            }
          }else if(!scale_horizontal_fix & !scale_vertical_fix){

            SCALE_HORIZONTAL <- exp(theta[2])
            SCALE_VERTICAL <- exp(theta[3])
            A1 <- exp(theta[4]) * 1e-3
            B1 <- theta[5] * 1e-3

            if(splines_degree == 0){
              C1_coef <- theta[6]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                A2 <- theta[7] * 1e-3
                B2 <- theta[8] * 1e-3
                C2_coef <- theta[9]
                D2 <- init_d2
              }else{
                D1 <- theta[7]
                A2 <- theta[8] * 1e-3
                B2 <- theta[9] * 1e-3
                C2_coef <- theta[10]
                D2 <- theta[11]
              }
            }else if(splines_degree > 0){
              C1_coef <- theta[5 + 1:nb1]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                A2 <- theta[5 + nb1 + 1] * 1e-3
                B2 <- theta[5 + nb1 + 2] * 1e-3
                C2_coef <- theta[5 + nb1 + 2 + 1:nb2]
                D2 <- init_d2
              }else{
                D1 <- theta[nb1 + 6]
                A2 <- theta[nb1 + 7] * 1e-3
                B2 <- theta[nb1 + 8] * 1e-3
                C2_coef <- theta[nb1 + 8 + 1:nb2]
                D2 <- theta[nb1 + nb2 + 9]
              }
            }
          }
        }
      }else if(beta_fix){

        BETA <- init_beta

        if(a1_fix & b1_fix & a2_fix & b2_fix){

          A1 <- init_a1
          B1 <- init_b1
          A2 <- init_a2
          B2 <- init_b2

          if(scale_horizontal_fix & scale_vertical_fix){

            SCALE_HORIZONTAL <- init_scale_horizontal
            SCALE_VERTICAL <- init_scale_vertical

            if(splines_degree == 0){
              C1_coef <- theta[1]
              if(d1_fix & d2_fix){
                D1 <- init_d1
                C2_coef <- theta[2]
                D2 <- init_d2
              }else{
                D1 <- theta[2]
                C2_coef <- theta[3]
                D2 <- theta[4]
              }
            }else if(splines_degree > 0){
              C1_coef <- theta[1:nb1]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                C2_coef <- theta[nb1 + 1:nb2]
                D2 <- init_d2
              }else{
                D1 <- theta[nb1 + 1]
                C2_coef <- theta[nb1 + 1 + 1:nb2]
                D2 <- theta[nb1 + nb2 + 2]
              }
            }
          }else if(!scale_horizontal_fix & !scale_vertical_fix){

            SCALE_HORIZONTAL <- exp(theta[1])
            SCALE_VERTICAL <- exp(theta[2])

            if(splines_degree == 0){

              C1_coef <- theta[3]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                C2_coef <- theta[4]
                D2 <- init_d2
              }else{
                D1 <- theta[4]
                C2_coef <- theta[5]
                D2 <- theta[6]
              }
            }else if(splines_degree > 0){

              C1_coef <- theta[2 + 1:nb1]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                C2_coef <- theta[nb1 + 2 + 1:nb2]
                D2 <- init_d2
              }else{
                D1 <- theta[nb1 + 3]
                C2_coef <- theta[nb1 + 3 + 1:nb2]
                D2 <- theta[nb1 + nb2 + 4]
              }
            }
          }
        }else{
          if(scale_horizontal_fix & scale_vertical_fix){

            SCALE_HORIZONTAL <- init_scale_horizontal
            SCALE_VERTICAL <- init_scale_vertical
            A1 <- exp(theta[1]) * 1e-3
            B1 <- theta[2] * 1e-3

            if(splines_degree == 0){
              C1_coef <- theta[3]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                A2 <- theta[4] * 1e-3
                B2 <- theta[5] * 1e-3
                C2_coef <- theta[6]
                D2 <- init_d2
              }else{
                D1 <- theta[4]
                A2 <- theta[5] * 1e-3
                B2 <- theta[6] * 1e-3
                C2_coef <- theta[7]
                D2 <- theta[8]
              }
            }else if(splines_degree > 0){
              C1_coef <- theta[2 + 1:nb1]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                A2 <- theta[nb1 + 3] * 1e-3
                B2 <- theta[nb1 + 4] * 1e-3
                C2_coef <- theta[nb1 + 4 + 1:nb2]
                D2 <- init_d2
              }else{
                D1 <- theta[nb1 + 3]
                A2 <- theta[nb1 + 4] * 1e-3
                B2 <- theta[nb1 + 5] * 1e-3
                C2_coef <- theta[nb1 + 5 + 1:nb2]
                D2 <- theta[nb1 + nb2 + 6]
              }
            }
          }else if(!scale_horizontal_fix & !scale_vertical_fix){

            SCALE_HORIZONTAL <- exp(theta[1])
            SCALE_VERTICAL <- exp(theta[2])
            A1 <- exp(theta[3]) * 1e-3
            B1 <- theta[4] * 1e-3

            if(splines_degree == 0){
              C1_coef <- theta[5]

              if(d1_fix & d2_fix){
                D1 <- init_d1
                A2 <- theta[6] * 1e-3
                B2 <- theta[7] * 1e-3
                C2_coef <- theta[8]
                D2 <- init_d2
              }else{
                D1 <- theta[6]
                A2 <- theta[7] * 1e-3
                B2 <- theta[8] * 1e-3
                C2_coef <- theta[9]
                D2 <- theta[10]
              }
            }else if(splines_degree > 0){
              C1_coef <- theta[4 + 1:nb1]

              if(d1_fix & d2_fix){

                D1 <- init_d1
                A2 <- theta[nb1 + 5] * 1e-3
                B2 <- theta[nb1 + 6] * 1e-3
                C2_coef <- theta[nb1 + 6 + 1:nb2]
                D2 <- init_d2
              }else{
                D1 <- theta[nb1 + 5]
                A2 <- theta[nb1 + 6] * 1e-3
                B2 <- theta[nb1 + 7] * 1e-3
                C2_coef <- theta[nb1 + 7 + 1:nb2]
                D2 <- theta[nb1 + nb2 + 8]
              }
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

    cholmat <- tryCatch(chol(cov_mat), error = function(a) numeric(0))
    if(length(cholmat) == 0){
      return(Inf)
    }

    Sigma_inv <- solve(cov_mat)
    Sigma_log_det <- 2 * sum(log(diag(t(cholmat))))

    out1 <- .5 * as.numeric(as.matrix(t(residuals) %*% Sigma_inv %*% residuals))
    out2 <- .5 * nrow(location) * 2 * log(2 * pi) + .5 * Sigma_log_det
    out <- out1 + out2

    cat(c("negloglik: ", round(out, 4)), '\n')

    return(out)
  }

  if(c1_fix & c2_fix){
    if(d1_fix & d2_fix){

      theta0 <- c(init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_a2, init_b2)

      fit <- nlm(NEGLOGLIK, theta0, hessian = T,
                 print.level = 2, iterlim = iterlim, stepmax = stepmax)

      BETA <- 1 / (1 + exp(-fit$estimate[1]))
      SCALE_HORIZONTAL <- exp(fit$estimate[2])
      SCALE_VERTICAL <- exp(fit$estimate[3])
      A1 <- exp(fit$estimate[4]) * 1e-3
      B1 <- fit$estimate[5] * 1e-3
      C1_coef <- init_c1_coef
      D1 <- init_d1

      A2 <- fit$estimate[6] * 1e-3
      B2 <- fit$estimate[7] * 1e-3
      C2_coef <- init_c2_coef
      D2 <- init_d2

      est_gradient <- fit$gradient
      est_hessian <- fit$hessian
      j <- diag(c(exp(-fit$estimate[1]) / (1 + exp(-fit$estimate[1]))^2, exp(fit$estimate[2:3]), exp(fit$estimate[4]) * 1e-3, rep(1e-3, 3)), nrow = length(theta0), ncol = length(theta0))
      fisher_info<-solve(est_hessian)
      variance <- j %*% fisher_info %*% t(j)
      est_sd <- sqrt(diag(variance))

      results <- list(convergence_code = fit$code, iterations = fit$iterations, loglikelihood_value = -fit$minimum, theta = fit$estimate,
                      est_beta = BETA, est_beta_sd = est_sd[1], est_scale_horizontal = SCALE_HORIZONTAL, est_scale_horizontal_sd = est_sd[2],
                      est_scale_vertical = SCALE_VERTICAL, est_scale_vertical_sd = est_sd[3], est_a1 = A1, est_a1_sd = est_sd[4], est_b1 = B1, est_b1_sd = est_sd[5],
                      est_c1_coef = C1_coef, est_d1 = D1, est_a2 = A2, est_a2_sd = est_sd[6], est_b2 = B2, est_b2_sd = est_sd[7], est_c2_coef = C2_coef, est_d2 = D2,
                      splines_degree = splines_degree, knots1 = knots1, knots2 = knots2)
    }
  }else{
    if(!beta_fix){
      if(a1_fix & b1_fix & a2_fix & b2_fix){
        if(scale_horizontal_fix & scale_vertical_fix){
          if(d1_fix & d2_fix){
            fit <- optim(par = c(init_beta, init_c1_coef, init_c2_coef), fn = NEGLOGLIK, control = list(trace = 5, iterlim = iterlim), hessian = T)
          }else{
            fit <- optim(par = c(init_beta, init_c1_coef, init_d1, init_c2_coef, init_d2), fn = NEGLOGLIK, control = list(trace = 5, iterlim = iterlim), hessian = T)
          }
        }else{
          if(d1_fix & d2_fix){

            theta0 <- c(init_beta, init_scale_horizontal, init_scale_vertical, init_c1_coef, init_c2_coef)

            fit <- nlm(NEGLOGLIK, theta0, hessian = hessian,
                       print.level = 2, iterlim = iterlim, stepmax = stepmax)

            BETA <- 1 / (1 + exp(-fit$estimate[1]))
            SCALE_HORIZONTAL <- exp(fit$estimate[2])
            SCALE_VERTICAL <- exp(fit$estimate[3])
            A1 <- init_a1
            B1 <- init_b1
            C1_coef <- fit$estimate[3 + 1:length(init_c1_coef)]
            D1 <- init_d1

            A2 <- init_a2
            B2 <- init_b2
            C2_coef <- fit$estimate[3 + length(init_c1_coef) + 1:length(init_c2_coef)]
            D2 <- init_d2

            est_gradient <- fit$gradient
            est_hessian <- fit$hessian
            j <- diag(c(exp(-fit$estimate[1]) / (1 + exp(-fit$estimate[1]))^2, exp(fit$estimate[2:3]), rep(1, length(init_c1_coef) + length(init_c2_coef))), nrow = length(theta0), ncol = length(theta0))
            fisher_info<-solve(est_hessian)
            variance <- j %*% fisher_info %*% t(j)
            est_sd <- sqrt(diag(variance))

            results <- list(convergence_code = fit$code, iterations = fit$iterations, loglikelihood_value = -fit$minimum, theta = fit$estimate,
                            est_beta = BETA, est_beta_sd = est_sd[1], est_scale_horizontal = SCALE_HORIZONTAL, est_scale_horizontal_sd = est_sd[2],
                            est_scale_vertical = SCALE_VERTICAL, est_scale_vertical_sd = est_sd[3], est_a1 = A1, est_b1 = B1,
                            est_c1_coef = C1_coef, est_c1_coef_sd = est_sd[3 + 1:length(init_c1_coef)], est_d1 = D1, est_a2 = A2, est_b2 = B2,
                            est_c2_coef = C2_coef, est_c2_coef_sd = est_sd[3 + length(init_c1_coef) + 1:length(init_c2_coef)], est_d2 = D2,
                            splines_degree = splines_degree, knots1 = knots1, knots2 = knots2)

          }else{
            fit <- optim(par = c(init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1_coef, init_d1, init_a2, init_b2, init_c2_coef, init_d2), fn = NEGLOGLIK, control = list(trace = 5, iterlim = iterlim), hessian = T)
          }
        }
      }else{
        if(scale_horizontal_fix & scale_vertical_fix){
          if(d1_fix & d2_fix){

            theta0 <- c(init_beta, init_a1, init_b1, init_c1_coef, init_a2, init_b2, init_c2_coef)

            fit <- nlm(NEGLOGLIK, theta0, hessian = hessian,
                       print.level = 2, iterlim = iterlim, stepmax = stepmax)

            BETA <- 1 / (1 + exp(-fit$estimate[1]))
            SCALE_HORIZONTAL <- init_scale_horizontal
            SCALE_VERTICAL <- init_scale_vertical
            A1 <- exp(fit$estimate[2]) * 1e-3
            B1 <- fit$estimate[3] * 1e-3

            C1_coef <- fit$estimate[3 + 1:length(init_c1_coef)]
            D1 <- init_d1

            A2 <- fit$estimate[3 + length(init_c1_coef) + 1] * 1e-3
            B2 <- fit$estimate[3 + length(init_c1_coef) + 2] * 1e-3
            C2_coef <- fit$estimate[3 + length(init_c1_coef) + 2 + 1:length(init_c2_coef)]
            D2 <- init_d2

            est_gradient <- fit$gradient
            est_hessian <- fit$hessian
            j <- diag(c(exp(-fit$estimate[1]) / (1 + exp(-fit$estimate[1]))^2, exp(fit$estimate[2]) * 1e-3, 1e-3, rep(1, length(init_c1_coef)), rep(1e-3, 2), rep(1, length(init_c2_coef))), nrow = length(theta0), ncol = length(theta0))
            fisher_info<-solve(est_hessian)
            variance <- j %*% fisher_info %*% t(j)
            est_sd <- sqrt(diag(variance))

            results <- list(convergence_code = fit$code, iterations = fit$iterations, loglikelihood_value = -fit$minimum, theta = fit$estimate,
                            est_beta = BETA, est_beta_sd = est_sd[1], est_scale_horizontal = SCALE_HORIZONTAL, est_scale_vertical = SCALE_VERTICAL,
                            est_a1 = A1, est_a1_sd = est_sd[2], est_b1 = B1, est_b1_sd = est_sd[3],
                            est_c1_coef = C1_coef, est_c1_coef_sd = est_sd[3 + 1:length(init_c1_coef)], est_d1 = D1, est_a2 = A2, est_a2_sd = est_sd[3 + length(init_c1_coef) + 1],
                            est_b2 = B2, est_b2_sd = est_sd[3 + length(init_c1_coef) + 2], est_c2_coef = C2_coef, est_c1_coef_sd = est_sd[3 + length(init_c1_coef) + 2 + 1:length(init_c2_coef)], est_d2 = D2,
                            splines_degree = splines_degree, knots1 = knots1, knots2 = knots2)

          }else{
            fit <- optim(par = c(init_beta, init_a1, init_b1, init_c1_coef, init_d1, init_a2, init_b2, init_c2_coef, init_d2), fn = NEGLOGLIK, control = list(trace = 5, iterlim = iterlim), hessian = T)
          }
        }else if(!scale_horizontal_fix & !scale_vertical_fix){
          if(d1_fix & d2_fix){

            theta0 <- c(init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1_coef, init_a2, init_b2, init_c2_coef)

            fit <- nlm(NEGLOGLIK, theta0, hessian = hessian,
                       print.level = 2, iterlim = iterlim, stepmax = stepmax)

            BETA <- 1 / (1 + exp(-fit$estimate[1]))
            SCALE_HORIZONTAL <- exp(fit$estimate[2])
            SCALE_VERTICAL <- exp(fit$estimate[3])
            A1 <- exp(fit$estimate[4]) * 1e-3
            B1 <- fit$estimate[5] * 1e-3

            C1_coef <- fit$estimate[5 + 1:length(init_c1_coef)]
            D1 <- init_d1

            A2 <- fit$estimate[5 + length(init_c1_coef) + 1] * 1e-3
            B2 <- fit$estimate[5 + length(init_c1_coef) + 2] * 1e-3
            C2_coef <- fit$estimate[5 + length(init_c1_coef) + 2 + 1:length(init_c2_coef)]
            D2 <- init_d2

            est_gradient <- fit$gradient
            est_hessian <- fit$hessian
            j <- diag(c(exp(-fit$estimate[1]) / (1 + exp(-fit$estimate[1]))^2, exp(fit$estimate[2:3]), exp(fit$estimate[4]) * 1e-3, 1e-3, rep(1, length(init_c1_coef)), rep(1e-3, 2), rep(1, length(init_c2_coef))), nrow = length(theta0), ncol = length(theta0))
            fisher_info<-solve(est_hessian)
            variance <- j %*% fisher_info %*% t(j)
            est_sd <- sqrt(diag(variance))

            results <- list(convergence_code = fit$code, iterations = fit$iterations, loglikelihood_value = -fit$minimum, theta = fit$estimate,
                            est_beta = BETA, est_beta_sd = est_sd[1], est_scale_horizontal = SCALE_HORIZONTAL, est_scale_horizontal_sd = est_sd[2],
                            est_scale_vertical = SCALE_VERTICAL, est_scale_vertical_sd = est_sd[3],
                            est_a1 = A1, est_a1_sd = est_sd[4], est_b1 = B1, est_b1_sd = est_sd[5],
                            est_c1_coef = C1_coef, est_c1_coef_sd = est_sd[5 + 1:length(init_c1_coef)], est_d1 = D1, est_a2 = A2, est_a2_sd = est_sd[5 + length(init_c1_coef) + 1],
                            est_b2 = B2, est_b2_sd = est_sd[5 + length(init_c1_coef) + 2], est_c2_coef = C2_coef, est_c1_coef_sd = est_sd[5 + length(init_c1_coef) + 2 + 1:length(init_c2_coef)], est_d2 = D2,
                            splines_degree = splines_degree, knots1 = knots1, knots2 = knots2)

          }else{
            fit <- optim(par = c(init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1_coef, init_d1, init_a2, init_b2, init_c2_coef, init_d2), fn = NEGLOGLIK, control = list(trace = 5, iterlim = iterlim), hessian = hessian)
          }
        }
      }
    }else{
      if(a1_fix & b1_fix & a2_fix & b2_fix){
        if(scale_horizontal_fix & scale_vertical_fix){
          if(d1_fix & d2_fix){

            theta0 <- c(init_c1_coef, init_c2_coef)

            fit <- nlm(NEGLOGLIK, theta0, hessian = hessian,
                       print.level = 2, iterlim = iterlim, stepmax = stepmax)

            BETA <- init_beta
            SCALE_HORIZONTAL <- init_scale_horizontal
            SCALE_VERTICAL <- init_scale_vertical
            A1 <- init_a1
            B1 <- init_b1
            C1_coef <- fit$estimate[1:length(init_c1_coef)]
            D1 <- init_d1

            A2 <- init_a2
            B2 <- init_b2
            C2_coef <- fit$estimate[length(init_c1_coef) + 1:length(init_c2_coef)]
            D2 <- init_d2

            est_gradient <- fit$gradient
            est_hessian <- fit$hessian
            j <- diag(c(rep(1, length(init_c1_coef) + length(init_c2_coef))), nrow = length(theta0), ncol = length(theta0))
            fisher_info<-solve(est_hessian)
            variance <- j %*% fisher_info %*% t(j)
            est_sd <- sqrt(diag(variance))

            results <- list(convergence_code = fit$code, iterations = fit$iterations, loglikelihood_value = -fit$minimum, theta = fit$estimate,
                            est_beta = BETA, est_scale_horizontal = SCALE_HORIZONTAL, est_scale_vertical = SCALE_VERTICAL, est_a1 = A1, est_b1 = B1,
                            est_c1_coef = C1_coef, est_c1_coef_sd = est_sd[1:length(init_c1_coef)], est_d1 = D1, est_a2 = A2, est_b2 = B2,
                            est_c2_coef = C2_coef, est_c2_coef_sd = est_sd[length(init_c1_coef) + 1:length(init_c2_coef)], est_d2 = D2,
                            splines_degree = splines_degree, knots1 = knots1, knots2 = knots2)

          }else{
            fit <- optim(par = c(init_c1_coef, init_d1, init_c2_coef, init_d2), fn = NEGLOGLIK, control = list(trace = 5, iterlim = iterlim), hessian = T)
          }
        }else{
          if(d1_fix & d2_fix){
            fit <- optim(par = c(init_scale_horizontal, init_scale_vertical, init_c1_coef, init_c2_coef), fn = NEGLOGLIK, control = list(trace = 5, iterlim = iterlim), hessian = T)
          }else{
            fit <- optim(par = c(init_scale_horizontal, init_scale_vertical, init_c1_coef, init_d1, init_c2_coef, init_d2), fn = NEGLOGLIK, control = list(trace = 5, iterlim = iterlim), hessian = T)
          }
        }
      }else{
        if(scale_horizontal_fix & scale_vertical_fix){
          if(d1_fix & d2_fix){
            fit <- optim(par = c(init_a1, init_b1, init_c1_coef, init_a2, init_b2, init_c2_coef), fn = NEGLOGLIK, control = list(trace = 5, iterlim = iterlim), hessian = T)
          }else{
            fit <- optim(par = c(init_a1, init_b1, init_c1_coef, init_d1, init_a2, init_b2, init_c2_coef, init_d2), fn = NEGLOGLIK, control = list(trace = 5, iterlim = iterlim), hessian = T)
          }
        }else{
          if(d1_fix & d2_fix){
            fit <- optim(par = c(init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1_coef, init_a2, init_b2, init_c2_coef), fn = NEGLOGLIK, control = list(trace = 5, iterlim = iterlim), hessian = T)
          }else{
            fit <- optim(par = c(init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1_coef, init_d1, init_a2, init_b2, init_c2_coef, init_d2), fn = NEGLOGLIK, control = list(trace = 5, iterlim = iterlim), hessian = T)
          }
        }
      }
    }
  }
  return(results)
}

#' Predict values at unsampled locations
#'
#' @usage predict_bi_differential(residuals, location, location_new, est_beta,
#' est_scale_horizontal, est_scale_vertical, est_a1, est_b1, est_c1, est_d1,
#' est_a2, est_b2, est_c2, est_d2, radius, splines_degree, knots1, knots2)
#'
#' @param residuals A \eqn{2 n} vector of residuals with
#' variable 1 and 2 residuals appended next to each other.
#' @param location An \eqn{n \times 3} matrix of coordinates.
#' @param location_new An \eqn{n \times 3} matrix of coordinates where prediction is required.
#' @param est_beta A numeric constant indicating the value for the
#' colocated correlation parameter.
#' @param est_scale_horizontal A numeric constant indicating the value for the
#' horizontal scale parameter.
#' @param est_scale_vertical A numeric constant indicating the value for the
#' vertical scale parameter.
#' @param est_a1 A numeric constant indicating the value for the
#' anisotropy in latitude parameter associated with variable 1.
#' @param est_b1 A numeric constant indicating the value for the
#' anisotropy in longitude parameter associated with variable 1.
#' @param est_c1 A numeric vector indicating the values for the
#' nonstationary parameters with depth associated with variable 1.
#' @param est_d1 A numeric constant indicating the value for the
#' variance parameter from the fully isotropic component associated with variable 1.
#' @param est_a2 A numeric constant indicating the value for the
#' anisotropy in latitude parameter associated with variable 2.
#' @param est_b2 A numeric constant indicating the value for the
#' anisotropy in longitude parameter associated with variable 2.
#' @param est_c2 A numeric vector indicating the values for the
#' nonstationary parameters with depth associated with variable 2.
#' @param est_d2 A numeric constant indicating the value for the
#' variance parameter from the fully isotropic component associated with variable 2.
#' @param radius A numeric constant indicating the radius of the sphere.
#' @param splines_degree A number indicating the degree of the splines.
#' @param knots1 A vector of knot locations for variable 1.
#' @param knots2 A vector of knot locations for variable 2.
#'
#' @return A vector of predictions.
#'
#' @useDynLib DiffOp, .registration=TRUE
#'
#' @author Mary Lai Salvana \email{yourlainess@gmail.com}
#'
#' @export
predict_bi_differential <- function(residuals, location, location_new, est_beta, est_scale_horizontal, est_scale_vertical, est_a1, est_b1, est_c1, est_d1, est_a2, est_b2, est_c2, est_d2, radius, splines_degree = 2, knots1, knots2){

  BETA <- est_beta
  SCALE_HORIZONTAL <- est_scale_horizontal
  SCALE_VERTICAL <- est_scale_vertical
  A1 <- est_a1
  B1 <- est_b1
  A2 <- est_a2
  B2 <- est_b2
  C1_coef <- est_c1
  C2_coef <- est_c2
  D1 <- est_d1
  D2 <- est_d2

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

  INDEX_OUT = c(nrow(location) + 1:nrow(location_new), nrow(location) + nrow(location_new) + nrow(location) + 1:nrow(location_new))

  Sigma_insample_dd <- cov_mat[-INDEX_OUT, -INDEX_OUT]
  Sigma_outsample_dd <- cov_mat[INDEX_OUT, INDEX_OUT]
  Sigma_inoutsample_dd <- cov_mat[-INDEX_OUT, INDEX_OUT]

  Sigma_inv <- solve(Sigma_insample_dd)

  predictions <- t(Sigma_inoutsample_dd) %*% Sigma_inv %*% matrix(residuals, ncol = 1)

  return(predictions)
}


