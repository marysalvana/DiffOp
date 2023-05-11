
library(broman)

source("/Users/laisalvana/Documents/documentation_pbdR/codes/R/yarger/R/get_profile_data_subset.R")

############ PREPARING RDATA FOR EACH LOCAL REGION ############

load("/Users/laisalvana/Documents/documentation_pbdR/MERRA2/data/jan_march_residuals.RData")

ref_lat <- c(40, 0, -40, 40, 0, -40)
ref_long <- c(-175, -175, -175, -30, -30, -30)
REFERENCE_LOCATIONS_MATRIX <- cbind(ref_long, ref_lat)

for(REF_LOC_IND in 1:1){
#for(REF_LOC_IND in 1:nrow(REFERENCE_LOCATIONS_MATRIX)){

  REFERENCE_LOCATION = REFERENCE_LOCATIONS_MATRIX[REF_LOC_IND, ]

  ReferenceLongitude <- REFERENCE_LOCATION[1]
  ReferenceLatitude <- REFERENCE_LOCATION[2]

  RADIUS = 900
  df <- get_profile_data_subset(long = ReferenceLongitude, lat = ReferenceLatitude, day = 45.25, h_time = 45.25, RG_def = FALSE,
                                h_space = RADIUS, mode = 'all',  min_prof = 5, exclusion = TRUE)

  sorted_df = df[order(df$total_dist_km), ]

  profile_num_unique <- unique(sorted_df$profile_num)

  num_sub_floats = 50
  num_vertical_obs = 50

  set.seed(1234)
  profile_num_unique <- c(profile_num_unique[1], sample(profile_num_unique, num_sub_floats - 1))

  keep_in <- NULL
  for(prof in profile_num_unique){
    ind1 <- which(sorted_df$profile_num == prof)
    #ind1 <- which(sorted_df$profile_num == prof & sorted_df$pressure <= 500)
    set.seed(5678)
    keep_in <- c(keep_in, sort(sample(ind1, min(length(ind1), num_vertical_obs))))
  }

  residuals <- sorted_df[keep_in, ]

  n <- nrow(residuals)

  cat('REF_LOC_IND: ', REF_LOC_IND, '   n: ', n, '\n')

  if(!is.null(n)){
    if(n > 500){
      file_name <- 'argo_ref_loc_new'

      TemperatureResiduals = residuals$temperature
      SalinityResiduals = residuals$salinity
      Pressure = residuals$pressure
      Latitude = residuals$latitude
      Longitude = residuals$longitude
      ProfileNumber = residuals$profile_num
      DistanceKmFromRefLoc = residuals$total_dist_km

      save(TemperatureResiduals, SalinityResiduals, Pressure, Latitude, Longitude, ReferenceLatitude, ReferenceLongitude, ProfileNumber, DistanceKmFromRefLoc, file = paste('/Users/laisalvana/Documents/documentation_pbdR/MERRA2/data/', file_name, REF_LOC_IND, '.RData', sep = ""))

    }
  }
}


for(REF_LOC_IND in 1:1){
  #for(REF_LOC_IND in 1:nrow(REFERENCE_LOCATIONS_MATRIX)){

  REFERENCE_LOCATION = REFERENCE_LOCATIONS_MATRIX[REF_LOC_IND, ]

  long <- REFERENCE_LOCATION[1]
  lat <- REFERENCE_LOCATION[2]

  RADIUS = 900
  df <- get_profile_data_subset(long = long, lat = lat, day = 45.25, h_time = 45.25, RG_def = FALSE,
                                h_space = RADIUS, mode = 'all',  min_prof = 5, exclusion = TRUE)

  sorted_df = df[order(df$total_dist_km), ]

  residuals <- sorted_df

  n <- nrow(residuals)

  cat('REF_LOC_IND: ', REF_LOC_IND, '   n: ', n, '\n')

  if(!is.null(n)){
    if(n > 500){

      Y <- as.matrix(residuals[, c('temperature', 'salinity')])

      #4 column location matrix: x, y, z, time
      locs <- cbind(residuals$longitude, residuals$latitude, residuals$pressure, rep(0, nrow(residuals)))

      setClass("RealData", representation(measurements = "matrix", locs = "matrix", num_locs = "numeric", num_time = "numeric", ref_lon = "numeric", ref_lat = "numeric"))

      DAT <- new("RealData", measurements = matrix(c(Y), ncol = 1), locs = locs, num_locs = nrow(residuals), num_time = 1, ref_lon = long, ref_lat = lat)

      file_name <- 'argo_ref_loc_full'
      save(DAT, file = paste('/Users/laisalvana/Documents/documentation_pbdR/MERRA2/data/', file_name, REF_LOC_IND, '.RData', sep = ""))
    }
  }
}


############ COMPUTING EMPIRICAL VARIANCE AND CORRELATION ALONG THE VERTICAL ############

REF_LOC_IND = 1

file_name <- 'argo_ref_loc_full'
load(paste('/Users/laisalvana/Documents/documentation_pbdR/MERRA2/data/', file_name, REF_LOC_IND, '.RData', sep = ""))

loc3d <- DAT@locs
Z <- DAT@measurements

emp_vals <- matrix(0, ncol = 3, nrow = 20)
loc3d_for_empirical <- cbind(ref_long[1], ref_lat[1], seq(0, 1900, by = 100))

for(PRES in 1:20){
  ind <- which(loc3d[, 3] > (PRES - 1) * 100 & loc3d[, 3] <= PRES * 100)
  emp_vals[PRES, 1] <- var(Z[ind])
  emp_vals[PRES, 2] <- var(Z[nrow(loc3d) + ind])
  emp_vals[PRES, 3] <- cor(Z[ind], Z[nrow(loc3d) + ind])
}

ind <- which(!is.na(rowSums(emp_vals)))
EMPIRICAL_VALUES <- emp_vals[ind, ]
LOCATION <- loc3d_for_empirical[ind, ]


plot(LOCATION[, 3], EMPIRICAL_VALUES[, 1])
plot(LOCATION[, 3], EMPIRICAL_VALUES[, 2])
plot(LOCATION[, 3], EMPIRICAL_VALUES[, 3])

############ VISUALIZATION ############

file_name <- 'argo_ref_loc_new'
load(paste('/Users/laisalvana/Documents/documentation_pbdR/MERRA2/data/', file_name, REF_LOC_IND, '.RData', sep = ""))

profile_no_unique <- unique(ProfileNumber)
profile_no_closest <- which(DistanceKmFromRefLoc == min(DistanceKmFromRefLoc))

png(paste('/Users/laisalvana/Library/CloudStorage/GoogleDrive-yourlainess@gmail.com/My Drive/Work/Research/UH/bivariate_argo/figures/residuals_profile_plots.png', sep = ''), width = 8, height = 4, units='in', res = 300)

split.screen( rbind(c(0.05,0.99,0.14,0.99), c(0.98,0.99,0.14,0.99)))
split.screen( figs = c(1, 2), screen = 1)

screen(3)
par(pty = 's')
par(mai = c(0.2, 0.2, 0.2, 0.2))

#plot(0, 0, xlim = range(TemperatureResiduals), ylim = c(2000, 0), xlab = 'Temperature Residuals', ylab = 'Depth (meters)', type = 'n')
plot(0, 0, xlim = range(TemperatureResiduals), ylim = c(2000, 0), xlab = '', ylab = '', type = 'n')
mtext('Temperature Residuals', side = 1, line = 2)
mtext('Depth (meters)', side = 2, line = 2)

for(prof in profile_no_unique){
  subset <- which(ProfileNumber == prof)
  lines(TemperatureResiduals[subset], Pressure[subset], lwd = 0.5, col = 'gray')
}

lines(TemperatureResiduals[profile_no_closest], Pressure[profile_no_closest], lwd = 1.5, col = 'black')

screen(4)
par(pty = 's')
par(mai = c(0.2, 0.2, 0.2, 0.2))

#plot(0, 0, xlim = range(SalinityResiduals), ylim = c(2000, 0), xlab = 'Salinity Residuals', ylab = 'Depth (meters)', type = 'n')
plot(0, 0, xlim = range(SalinityResiduals), ylim = c(2000, 0), xlab = '', ylab = '', type = 'n')
mtext('Salinity Residuals', side = 1, line = 2)
mtext('Depth (meters)', side = 2, line = 2)

for(prof in profile_no_unique){
  subset <- which(ProfileNumber == prof)
  lines(SalinityResiduals[subset], Pressure[subset], lwd = 0.5, col = 'gray')
}

lines(SalinityResiduals[profile_no_closest], Pressure[profile_no_closest], lwd = 1.5, col = 'black')

close.screen( all=TRUE)
dev.off()

############ STEP 0: MLE FOR THE CLASSICAL MODEL ############

earthRadiusKm = 6371

INIT_BETA = INIT_A1 = INIT_B1 = INIT_A2 = INIT_B2 = 0
SCALE_HORIZONTAL = 0.002998395
SCALE_VERTICAL = 0.01644395

nn <- length(TemperatureResiduals)

ind <- 51:nn
ind_pred <- profile_no_closest

loc3d <- cbind(Longitude, Latitude, Pressure)
locs_insample <- loc3d[-ind_pred, ]
locs_outsample <- loc3d[ind_pred, ]

Z_insample <- cbind(DAT@measurements[ind, ], DAT@measurements[nn + ind, ])
Z_outsample <- cbind(DAT@measurements[ind_pred, ], DAT@measurements[nn + ind_pred, ])
Z <- c(rbind(Z_insample, Z_outsample))

KNOTS1 <- seq(0, 2100, length.out = 2)
KNOTS2 <- seq(0, 2100, length.out = 2)

basis1 <- bsplineBasis(locs_insample[, 3], 0, KNOTS1)
nb1 <- ncol(basis1)
basis2 <- bsplineBasis(locs_insample[, 3], 0, KNOTS2)
nb2 <- ncol(basis2)

set.seed(1235)
INIT_C1 <- runif(nb1, -0.5, 0.5)

set.seed(1236)
INIT_C2 <- runif(nb2, -0.5, 0.5)

est_params <- est_bi_differential(residuals = c(Z_insample), location = locs_insample, init_beta = INIT_BETA,
                                  init_scale_horizontal = log(SCALE_HORIZONTAL), init_scale_vertical = log(SCALE_VERTICAL),
                                  init_a1 = INIT_A1, init_b1 = INIT_B1, init_c1 = INIT_C1, init_a2 = INIT_A2, init_b2 = INIT_B2, init_c2 = INIT_C2,
                                  basis1 = basis1, basis2 = basis2, radius = earthRadiusKm, splines_degree = 0, knots1 = KNOTS1, knots2 = KNOTS2)

theta = c(0.722525,-1.751571,-4.502191,0.000488,3.9e-05,2.687413,0.129758,2.687051,5.3e-05,-1.6e-05,-1.206159,2.317378,-0.100077)

#RERUN FROM STOPPING POINT THETA

est_params <- est_bi_differential(residuals = c(Z_insample), location = locs_insample, init_beta = theta[1],
                                  init_scale_horizontal = theta[2], init_scale_vertical = theta[3],
                                  init_a1 = theta[4], init_b1 = theta[5], init_c1 = theta[5 + 1:nb1], init_a2 = theta[nb1 + 6], init_b2 = theta[nb1 + 7], init_c2 = theta[nb1 + 7 + 1:nb2],
                                  basis1 = basis1, basis2 = basis2, radius = earthRadiusKm, splines_degree = 0, knots1 = KNOTS1, knots2 = KNOTS2)

#MLE
theta = c(0.72250994,-1.75153982,-4.05269659,0.00031493,0.0003175,2.68750921,0.12975717,2.68716172,7.864e-05,-3.818e-05,-1.2060537,2.317533,-0.09990902)
#neg loglik = -5818.5248

############ STEP 1: WLS ############

INIT_BETA = theta[1]
INIT_SCALE_HORIZONTAL = exp(theta[2])
INIT_SCALE_VERTICAL = exp(theta[3])
INIT_A1 = theta[4]
INIT_A2 = theta[9]
INIT_B1 = theta[5]
INIT_B2 = theta[10]
INIT_D1 = INIT_D2 = 0

basis1 <- bsplineBasis(LOCATION[, 3], 4, KNOTS1)
nb1 <- ncol(basis1)
basis2 <- bsplineBasis(LOCATION[, 3], 4, KNOTS2)
nb2 <- ncol(basis2)

set.seed(1235)
INIT_C1 <- runif(nb1, -100, 100)

set.seed(1236)
INIT_C2 <- runif(nb2, -10, 10)

est_params <- est_bi_differential_wls(empirical_values = EMPIRICAL_VALUES, location = LOCATION, init_beta = INIT_BETA,
                                      init_scale_horizontal = INIT_SCALE_HORIZONTAL, init_scale_vertical = INIT_SCALE_VERTICAL, init_scale_horizontal_fix = T, init_scale_vertical_fix = T,
                                      init_a1 = INIT_A1, init_b1 = INIT_B1, init_c1 = INIT_C1, init_a2 = INIT_A2, init_b2 = INIT_B2, init_c2 = INIT_C2, init_beta_fix = T,
                                      basis1 = basis1, basis2 = basis2, radius = earthRadiusKm, splines_degree = 4, knots1 = KNOTS1, knots2 = KNOTS2)

#RERUN FOR WLS
theta = c(0,0.000137,-43.478744,106.037114,-33.149147,75.450859,58.867069,-11.714446,8.646101,2e-06,-1.9e-05,-5.850154,14.24039,-6.484243,-5.210463,1.811206,-0.15922,5.527517)
theta = c(2.2e-07,6.804e-05,-51.52184895,71.74782542,-66.36078398,81.95555724,69.22938586,-10.68999432,6.89582243,2.46e-06,-3.811e-05,-5.82206381,13.2819233,-3.61022406,-9.577206,-3.84409763,1.44576138,8.03768149)

INIT_A1 = theta[1]
INIT_A2 = theta[nb1 + 3]
INIT_B1 = theta[2]
INIT_B2 = theta[nb1 + 4]
INIT_C1 <- theta[2 + 1:nb1]
INIT_C2 <- theta[nb1 + 4 + 1:nb2]

############ STEP 2: MLE ############

basis1 <- bsplineBasis(locs_insample[, 3], 4, KNOTS1)
nb1 <- ncol(basis1)
basis2 <- bsplineBasis(locs_insample[, 3], 4, KNOTS2)
nb2 <- ncol(basis2)

#initial parameters: n = 2500, including the prediction locations
theta = c(0.7, log(INIT_SCALE_HORIZONTAL), log(INIT_SCALE_VERTICAL), 0,0.000137,-43.478744,106.037114,-33.149147,75.450859,58.867069,-11.714446,8.646101,2e-06,-1.9e-05,-5.850154,14.24039,-6.484243,-5.210463,1.811206,-0.15922,5.527517)

#initial parameters: n = 2500, excluding the prediction locations
theta = c(INIT_BETA, log(INIT_SCALE_HORIZONTAL), log(INIT_SCALE_VERTICAL), 2.2e-07,6.804e-05,-51.52184895,71.74782542,-66.36078398,81.95555724,69.22938586,-10.68999432,6.89582243,2.46e-06,-3.811e-05,-5.82206381,13.2819233,-3.61022406,-9.577206,-3.84409763,1.44576138,8.03768149)

#initial parameters: n = 2500, excluding the prediction locations with additional reruns for WLS
theta = c(INIT_BETA, log(INIT_SCALE_HORIZONTAL), log(INIT_SCALE_VERTICAL), 5.13e-06,2.42e-05,-52.73111221,71.88003595,-62.47084003,83.34465397,67.75742297,-12.37265794,3.18738643,-1.54e-06,-1.543e-05,-4.87730718,12.34853351,-7.01972981,-6.86201992,7.72033671,-5.78460723,8.77574025)

est_params <- est_bi_differential(residuals = c(Z_insample), location = locs_insample, init_beta = theta[1],
                                  init_scale_horizontal = theta[2], init_scale_vertical = theta[3],
                                  init_a1 = theta[4], init_b1 = theta[5], init_c1 = theta[5 + 1:nb1], init_a2 = theta[nb1 + 6], init_b2 = theta[nb1 + 7], init_c2 = theta[nb1 + 7 + 1:nb2],
                                  basis1 = basis1, basis2 = basis2, radius = earthRadiusKm, splines_degree = 4, knots1 = KNOTS1, knots2 = KNOTS2)

#MLE: n = 2500, including the prediction locations
theta = c(0.700363,-1.751355,-4.502849,0,0.000167,-43.478301,106.036813,-33.149263,75.451049,58.867345,-11.714336,8.64584,9e-06,-8e-06,-5.849888,22.464088,-6.484308,-3.223079,1.811588,-0.159048,5.527556)

#MLE: n = 2500, excluding the prediction locations
theta = c(0.73583128,-1.75206027,-4.05319112,1.52e-06,0.00024333,-51.52124876,72.38970266,-66.36080013,81.95269361,69.22969203,-10.68909221,6.89646304,1.392e-05,-4.99e-06,-5.8216641,13.28222786,-3.6103994,-4.45451418,-1.42743919,1.44666584,8.0391838)
#neg loglik = -6858.6145

#MLE: n = 2500, excluding the prediction locations with additional reruns for WLS
theta = c(0.72545012,-1.64953152,-4.05563362,2.53e-06,0.00021368,-52.73222083,76.33190581,-62.47251572,83.34264221,67.75549558,-12.3718472,3.18806276,-1.366e-05,-1.92e-06,-4.87616132,12.34728918,-6.4992468,-3.60197678,7.72209309,-5.78397767,8.7763426)
#neg loglik = -7017.1351

#TODO: FIX codes for init_a1_fix = T, NLM
#############################################################

EST_BETA <- theta[1]
EST_SCALE_HORIZONTAL <- exp(theta[2])
EST_SCALE_VERTICAL <- exp(theta[3])
EST_A1 <- theta[4]
EST_B1 <- theta[5]
EST_C1 <- theta[5 + 1:nb1]
EST_A2 <- theta[nb1 + 6]
EST_B2 <- theta[nb1 + 7]
EST_C2 <- theta[nb1 + 7 + 1:nb2]

Z_pred <- predict_bi_differential(residuals = c(Z_insample), location = locs_insample, location_new = locs_outsample, masked_residuals = c(Z_outsample),
                                  est_beta = EST_BETA, est_scale_horizontal = EST_SCALE_HORIZONTAL, est_scale_vertical = EST_SCALE_VERTICAL,
                                  est_a1 = EST_A1, est_b1 = EST_B1, est_c1 = EST_C1, est_a2 = EST_A2, est_b2 = EST_B2, est_c2 = EST_C2,
                                  radius = earthRadiusKm, splines_degree = 4, knots1 = KNOTS1, knots2 = KNOTS2)

#############################################################

theta = c(0.7,0,0,-0.257408,0.015356,-0.400578,0.401536,0.338903,-0.076439,0.044499,0.245316,0.045463,-0.287664,0.455265,0.042894,0.01908,0,0,-0.318483,0.419723,-0.188277,-0.435848,-0.103578,0.331421,0.331023,0.067886,-0.451616,-0.441957,0.134718,0.352717,-0.319511)

theta = c(0.7,0.001,0.001,-0.257408,0.015356,-0.400578,0.401536,0.338903,-0.076439,0.044499,0.245316,0.045463,-0.287664,0.455265,0.042894,0.01908,-0.001,-0.001,-0.318483,0.419723,-0.188277,-0.435848,-0.103578,0.331421,0.331023,0.067886,-0.451616,-0.441957,0.134718,0.352717,-0.319511)

theta = c(0.63,0.004516,0.004516,-0.252892,0.019872,-0.396062,0.406052,0.343419,-0.071923,0.049016,0.249832,0.04998,-0.283148,0.459781,0.04741,0.023596,0.004516,0.004516,-0.313967,0.424239,-0.18376,-0.431332,-0.099062,0.335937,0.335539,0.072402,-0.4471,-0.43744,0.139234,0.357233,-0.384994)

theta = c(0.7,0,0,-2.574076,0.153559,-4.005783,4.015359,3.389029,-0.764389,0.444994,2.453162,0.454634,-2.876644,4.552647,0.42894,0.190804,0,0,0,-3.184829,4.197226,-1.882765,-4.358484,-1.035784,3.314213,3.31023,0.678864,-4.516163,-4.419566,1.347183,3.52717,-3.667473,0.227632)

theta = c(0.7,0,0,-2.574076,0.153559,-4.005783,4.015359,3.389029,-0.764389,0.444994,2.566979,0.454634,-2.876644,4.552647,0.42894,0.190804,0.113816,0,0,-3.184829,4.197226,-1.882765,-4.358484,-1.035784,3.314213,3.31023,0.678864,-4.516163,-4.419566,1.347183,3.52717,-3.895105,0.227632)

ref_lat <- c(40, 0, -40, 40, 0, -40)
ref_long <- c(-175, -175, -175, -30, -30, -30)

loc3d_eval <- cbind(ref_long[1], ref_lat[1], seq(0, 2000, length.out = 100))

KNOTS1 <- seq(0, max(loc3d[, 3]), length.out = 10)
KNOTS2 <- seq(0, max(loc3d[, 3]), length.out = 10)

new_basis1 <- bsplineBasis(loc3d_eval[, 3], 2, KNOTS1)
new_basis2 <- bsplineBasis(loc3d_eval[, 3], 2, KNOTS2)

plot_bi_differential(location = loc3d_eval, est_beta = theta[1],
                     est_scale_horizontal = 0.002998395, est_scale_vertical = 0.01644395,
                     est_a1 = theta[2], est_b1 = theta[3], est_c1 = theta[3 + 1:nb1], est_a2 = theta[nb1 + 4], est_b2 = theta[nb1 + 5], est_c2 = theta[nb1 + 5 + 1:nb2],
                     basis1 = new_basis1, basis2 = new_basis2, radius = earthRadiusKm)

