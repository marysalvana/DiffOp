

deg2rad <- function(deg){
  return(deg * pi / 180)
}

rad2deg <- function(rad){
  return(rad * 180 / pi)
}

distanceEarth <- function(lat1d, lon1d, lat2d, lon2d, radius) {
  lat1r = deg2rad(lat1d)
  lon1r = deg2rad(lon1d)
  lat2r = deg2rad(lat2d)
  lon2r = deg2rad(lon2d)
  u = sin((lat2r - lat1r)/2);
  v = sin((lon2r - lon1r)/2);
  return(2.0 * radius * sqrt(u^2 + cos(lat1r) * cos(lat2r) * v^2))
}

calculateDistance <- function(x1, y1, x2, y2, radius) {
  if(is.null(radius)){
    return(sqrt((x2 - x1)^2 + (y2 - y1)^2))
  }else{
    return (distanceEarth(x1, y1, x2, y2, radius))
  }
}

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
  alphaij=stationary_param[10]
  
  c1 = nonstationary_param1
  c2 = nonstationary_param2
  H = h_new(stationary_param[1], stationary_param[2], lat1d, lon1d, pres1, lat2d, lon2d, pres2, radius);
  H1 = h1(stationary_param[1], lat1d, lon1d, lat2d, lon2d, radius);
  H2 = h1(stationary_param[1], lat2d, lon1d, lat1d, lon2d, radius);
  H3 = h3(stationary_param[1], lat1d, lon1d, lat2d, lon2d, radius);
  H4 = h4(stationary_param[2], pres1, pres2);
  return(alphaij/4 * (a1 * a2 * H1 * H2 - b1 * b2 * H3^2 - c1 * c2 * H4^2 - a1 * b2 * H1 * H3
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
  alphaij=stationary_param[10]
  
  c1 = nonstationary_param1
  c2 = nonstationary_param2
  H12 = h12(stationary_param[1], lat1d, lon1d, lat2d, lon2d, radius)
  H13 = h13(stationary_param[1], lat1d, lon1d, lat2d, lon2d, radius)
  H23 = h13(stationary_param[1], lat2d, lon1d, lat1d, lon2d, radius)
  H33 = h33(stationary_param[1], lat1d, lon1d, lat2d, lon2d, radius)
  H44 = h44(stationary_param[2])
  return(-alphaij/2 * (a1 * a2 * H12 - b1 * b2 * H33 - c1 * c2 * H44 - a1 * b2 * H13 + a2 * b1 * H23) + 2 * alphaij*nu * d1 * d2) ## THIS LAST TERM LOOKS SUSPICIOUS
}

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

Legendre=function(x, degree){
  
  # temp= legendre(degree,2*x-1)
  # return(t(temp))
  temp = sapply(0 : degree, function(k) {
    if( k == 0 ) {
      return(legendre(k, x * 2 - 1))
    } else {
      return(legendre(k, x * 2 - 1)[1, ])
    }
  })
  return(temp)
}




uni_differential <- function(PARAM, fd_eval_mat_loc1, fd_eval_mat_loc2, LAT1D, LON1D, PRES1, LAT2D, LON2D, PRES2, radius){
  sigma1 = PARAM[1]
  sigma2=PARAM[2]
  SCALE_HORIZONTAL_SPACE = PARAM[3]
  SCALE_VERTICAL_SPACE = PARAM[4]
  smoothness = PARAM[5]
  beta=PARAM[6]
  
  a1 <- PARAM[7]
  b1 <- PARAM[8]
  d1 <- PARAM[9]
  a2 <- PARAM[10]
  b2 <- PARAM[11]
  d2 <- PARAM[12]
  c1 <- fd_eval_mat_loc1
  
  c2 <- fd_eval_mat_loc2
  con = 2^(smoothness - 1) * gamma(smoothness) ## COMMON SMOOTHNESS FOR NOW
  con = 1.0 / con
  con = beta*sigma1*sigma2* con ## BETA IS 1 IF UNIVARIATE
  
  alphaij=con
  
  expr <- sqrt(h_new(SCALE_HORIZONTAL_SPACE, SCALE_VERTICAL_SPACE, LAT1D, LON1D, PRES1, LAT2D, LON2D, PRES2, radius))
  
  
  STATIONARY_PARAM <- c(SCALE_HORIZONTAL_SPACE, SCALE_VERTICAL_SPACE, a1, b1, d1, a2, b2, d2, smoothness,alphaij)
  
  f <- expr^(smoothness) * besselK(expr, smoothness)
  
  diag(f)= 2^(smoothness-1)*gamma(smoothness) ##
  
  f_prime <- expr^(smoothness - 1) * besselK(expr, smoothness - 1)
  
  diag(f_prime)= 2^(smoothness-2)*gamma(smoothness-1) ##
  
  C1_val <- C1(stationary_param = STATIONARY_PARAM, nonstationary_param1 = c1, nonstationary_param2 = c2,
               lat1d = LAT1D, lon1d = LON1D, pres1 = PRES1, lat2d = LAT2D, lon2d = LON2D, pres2 = PRES2, radius)
  C2_val <- C2(stationary_param = STATIONARY_PARAM, nonstationary_param1 = c1, nonstationary_param2 = c2,
               lat1d = LAT1D, lon1d = LON1D, pres1 = PRES1, lat2d = LAT2D, lon2d = LON2D, pres2 = PRES2, radius)
  
  val <-  (C1_val * f_prime + C2_val * f + d1 * d2 * (expr^2 * f_prime + 2 * (smoothness) * f))
  ## diag(val) <- con * (diag(C1_val) + diag(C2_val)) + sigma_square * d1 * d2
  
  
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
#' a1, b1, c1, d1, a2, b2, c2, d2, radius_of_sphere, splines_degree,
#' inner_knots1, inner_knots2)
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
#' if splines_degree > 0, this is A numeric vector indicating the splines coefficients for the
#' nonstationary with depth parameter c1 associated with variable 1.
#' @param d2 A numeric constant indicating the variance parameter from the fully isotropic component associated with variable 2.
#' @param radius_of_sphere A numeric constant indicating the radius of the sphere.
#' @param splines_degree A number indicating the degree of the splines when
#' using splines to characterize the nonstationary parameters c1 and c2.
#' @param inner_knots1 A vector of knot locations for variable 1 when using splines to
#' characterize c1.
#' @param inner_knots2 A vector of knot locations for variable 2 when using splines to
#' characterize c2.
#'
#' @importFrom pracma legendre
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
#'                                radius_of_sphere = earthRadiusKm)
#'
#'
#' @export
cov_bi_differential=function (location, sigma1, sigma2, beta, scale_horizontal, scale_vertical, a1, 
                              b1, c1 = NULL, d1, a2, b2, c2 = NULL, d2, radius, Legen_degree = NULL, 
                              inner_knots1 = NULL, inner_knots2 = NULL, c1_coef = NULL, 
                              c2_coef = NULL) 
{
  LAT1D <- matrix(location[, 2], nrow(location), nrow(location), 
                  byrow = F)
  LON1D <- matrix(location[, 1], nrow(location), nrow(location), 
                  byrow = F)
  PRES1 <- matrix(location[, 3], nrow(location), nrow(location), 
                  byrow = F)
  LAT2D <- matrix(location[, 2], nrow(location), nrow(location), 
                  byrow = T)
  LON2D <- matrix(location[, 1], nrow(location), nrow(location), 
                  byrow = T)
  PRES2 <- matrix(location[, 3], nrow(location), nrow(location), 
                  byrow = T)
  if (!is.null(Legen_degree)) {
    if (Legen_degree == 0) {
      nb1 <- nb2 <- 1
    }
    else {
      ## basis1 <- bsplineBasis(location[, 3], splines_degree, 
      ##        inner_knots1)
      
      basis1=Legendre(location[,3], Legen_degree)
      
      nb1 <- ncol(basis1)
      ##  basis2 <- bsplineBasis(location[, 3], splines_degree, 
      ##                       inner_knots2)
      
      basis2=Legendre(location[,3], Legen_degree)
      nb2 <- ncol(basis2)
    }
    if (Legen_degree == 0) {
      c1 <- c1_coef
      c2 <- c2_coef
    }
    else if (Legen_degree > 0) {
      c1 <- basis1 %*% matrix(c1_coef, ncol = 1)
      c2 <- basis2 %*% matrix(c2_coef, ncol = 1)
      
    }
  }
  
  fd_eval_mat_loc1 <- matrix(c1, nrow(location), nrow(location), 
                             byrow = F)
  fd_eval_mat_loc2 <- matrix(c1, nrow(location), nrow(location), 
                             byrow = T)
  fd_eval2_mat_loc1 <- matrix(c2, nrow(location), nrow(location), 
                              byrow = F)
  fd_eval2_mat_loc2 <- matrix(c2, nrow(location), nrow(location), 
                              byrow = T)
  
  PARAM <- c(sigma1, sigma1,scale_horizontal, scale_vertical, 2, 1,a1, b1, 
             d1, a1, b1, d1)
  cov_val <- uni_differential(PARAM, fd_eval_mat_loc1, fd_eval_mat_loc2, 
                              LAT1D, LON1D, PRES1, LAT2D, LON2D, PRES2, radius)
  PARAM <- c(sigma2,sigma2, scale_horizontal, scale_vertical, 2,1, a2, b2, 
             d2, a2, b2, d2)
  
  cov_val2 <- uni_differential(PARAM, fd_eval2_mat_loc1, fd_eval2_mat_loc2, 
                               LAT1D, LON1D, PRES1, LAT2D, LON2D, PRES2, radius)
  PARAM <- c(sigma1,sigma2, scale_horizontal, scale_vertical, 2,beta, a1, 
             b1, d1, a2, b2, d2)
  cov_val3 <- uni_differential(PARAM, fd_eval_mat_loc1, fd_eval2_mat_loc2, 
                               LAT1D, LON1D, PRES1, LAT2D, LON2D, PRES2, radius)
  
  Sigma <- rbind(cbind(cov_val, cov_val3), cbind(t(cov_val3), 
                                                 cov_val2))
  
  
  return(Sigma)
}

