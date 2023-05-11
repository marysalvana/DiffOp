
library(broman)

source("/Users/laisalvana/Documents/documentation_pbdR/codes/R/yarger/R/get_profile_data_subset.R")

############ PREPARING RDATA FOR EACH LOCAL REGION ############

load("/Users/laisalvana/Documents/documentation_pbdR/MERRA2/data/jan_march_residuals.RData")

ref_lat <- c(40, 0, -40, 40, 0, -40)
ref_long <- c(-175, -175, -175, -30, -30, -30)
REFERENCE_LOCATIONS_MATRIX <- cbind(ref_long, ref_lat)

for(REF_LOC_IND in 5:5){
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

      file_name <- 'argo_ref_loc_full'

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

############ COMPUTING EMPIRICAL VARIANCE AND CORRELATION ALONG THE VERTICAL ############

REF_LOC_IND = 1

file_name <- 'argo_ref_loc_new'
load(paste('/Users/laisalvana/Documents/documentation_pbdR/MERRA2/data/', file_name, REF_LOC_IND, '.RData', sep = ""))

loc3d <- cbind(Longitude, Latitude, Pressure)

emp_vals <- matrix(0, ncol = 3, nrow = 200)
loc3d_for_empirical <- cbind(ref_long[1], ref_lat[1], seq(0, 2000 - 10, by = 10))

for(PRES in 1:200){
  ind <- which(loc3d[, 3] > (PRES - 1) * 10 & loc3d[, 3] <= PRES * 10)
  emp_vals[PRES, 1] <- var(TemperatureResiduals[ind])
  emp_vals[PRES, 2] <- var(SalinityResiduals[ind])
  emp_vals[PRES, 3] <- cor(TemperatureResiduals[ind], SalinityResiduals[ind])
}

ind <- which(!is.na(rowSums(emp_vals)))
EMPIRICAL_VALUES <- emp_vals[ind, ]
LOCATION <- loc3d_for_empirical[ind, ]

plot(LOCATION[, 3], EMPIRICAL_VALUES[, 1])
plot(LOCATION[, 3], EMPIRICAL_VALUES[, 2])
plot(LOCATION[, 3], EMPIRICAL_VALUES[, 3])

############ VISUALIZATION OF THE ARGO FLOAT MEASUREMENTS ############

profile_no_unique <- unique(ProfileNumber)
profile_no_closest <- which(DistanceKmFromRefLoc == min(DistanceKmFromRefLoc))

#png(paste('/Users/laisalvana/Library/CloudStorage/GoogleDrive-yourlainess@gmail.com/My Drive/Work/Research/UH/bivariate_argo/figures/residuals_profile_plots_ref_loc', REF_LOC_IND, '.png', sep = ''), width = 8, height = 4, units='in', res = 300)
png(paste('/Users/laisalvana/Library/CloudStorage/GoogleDrive-yourlainess@gmail.com/My Drive/Work/Research/UH/bivariate_argo/figures/residuals_profile_plots_ref_loc', REF_LOC_IND, '_with_predictions.png', sep = ''), width = 8, height = 4, units='in', res = 300)

split.screen( rbind(c(0.05,0.99,0.14,0.99), c(0.98,0.99,0.14,0.99)))
split.screen( figs = c(1, 2), screen = 1)

screen(3)
par(pty = 's')
par(mai = c(0.2, 0.2, 0.2, 0.4))

plot(0, 0, xlim = range(TemperatureResiduals), ylim = c(2000, 0), xlab = '', ylab = '', type = 'n')
mtext('Temperature Residuals', side = 1, line = 2)
mtext('Depth (meters)', side = 2, line = 2)

for(prof in profile_no_unique){
  subset <- which(ProfileNumber == prof)
  lines(TemperatureResiduals[subset], Pressure[subset], lwd = 0.5, col = 'gray')
}

lines(TemperatureResiduals[profile_no_closest], Pressure[profile_no_closest], lwd = 1.5, col = 'black')
lines(pred_temp, Pressure[profile_no_closest], lwd = 1.5, col = 'red')


screen(4)
par(pty = 's')
par(mai = c(0.2, 0.4, 0.2, 0.2))

plot(0, 0, xlim = range(SalinityResiduals), ylim = c(2000, 0), xlab = '', ylab = '', type = 'n')
mtext('Salinity Residuals', side = 1, line = 2)
mtext('Depth (meters)', side = 2, line = 2)

for(prof in profile_no_unique){
  subset <- which(ProfileNumber == prof)
  lines(SalinityResiduals[subset], Pressure[subset], lwd = 0.5, col = 'gray')
}

lines(SalinityResiduals[profile_no_closest], Pressure[profile_no_closest], lwd = 1.5, col = 'black')
lines(pred_sal, Pressure[profile_no_closest], lwd = 1.5, col = 'red')

close.screen( all=TRUE)
dev.off()

pred_temp = c(0.677461,0.667309,0.785966,0.937077,1.019026,1.086106,1.135802,1.171143,1.0085,0.823459,0.675698,0.644475,0.515659,0.494628,0.421473,0.390112,0.339016,0.31774,0.269555,0.191044,0.170523,0.159334,0.135797,0.069826,0.031291,-0.001154,-0.017648,-0.047624,-0.048391,-0.044989,-0.024426,-0.019837,-0.016294,-0.009158,-0.008799,-0.007058,-0.000583,0.002679,0.004968,0.00507,0.004393,0.00668,0.007878,0.008984,0.010055,0.012081,0.012578,0.011845,0.008378,0.00711)
pred_sal = c(0.060194,0.058126,0.069306,0.087689,0.098449,0.107977,0.115888,0.129617,0.118284,0.100995,0.086008,0.082605,0.067105,0.064275,0.053572,0.048551,0.039916,0.03625,0.028459,0.022193,0.020425,0.019426,0.017276,0.01123,0.007863,0.005172,0.003901,0.001848,0.001906,0.000914,-0.00027,-0.000246,-0.00011,0.000956,0.001038,0.000843,-0.000315,-0.000936,-0.001337,-0.001183,-0.000823,-0.000683,-0.000733,-0.000881,-0.001109,-0.001542,-0.001727,-0.001635,-0.000896,-0.000657)

############ COMPUTING THE EMPIRICAL COVARIANCE ############

earthRadiusKm = 6371

REF_LOC_IND = 1

file_name <- 'argo_ref_loc_new'
load(paste('/Users/laisalvana/Documents/documentation_pbdR/MERRA2/data/', file_name, REF_LOC_IND, '.RData', sep = ""))

loc3d <- cbind(Longitude, Latitude, Pressure)

LAT1D <- matrix(loc3d[, 2], nrow(loc3d), nrow(loc3d), byrow = F)
LON1D <- matrix(loc3d[, 1], nrow(loc3d), nrow(loc3d), byrow = F)
PRES1 <- matrix(loc3d[, 3], nrow(loc3d), nrow(loc3d), byrow = F)
LAT2D <- matrix(loc3d[, 2], nrow(loc3d), nrow(loc3d), byrow = T)
LON2D <- matrix(loc3d[, 1], nrow(loc3d), nrow(loc3d), byrow = T)
PRES2 <- matrix(loc3d[, 3], nrow(loc3d), nrow(loc3d), byrow = T)

#dist0 <- sqrt(h_new(0.0011, 0.03, LAT1D, LON1D, PRES1, LAT2D, LON2D, PRES2, radius = earthRadiusKm))
dist0 <- sqrt(h_new(0.009, 0.03, LAT1D, LON1D, PRES1, LAT2D, LON2D, PRES2, radius = earthRadiusKm))
kernel <- exp(-dist0)
kernel_sum <- rowSums(kernel)

denominator <- outer(sqrt(kernel_sum), sqrt(kernel_sum), '*')

TemperatureResiduals_matrix <- matrix(TemperatureResiduals, length(TemperatureResiduals), length(TemperatureResiduals), byrow = T)

emp_covariance1_temp <- sqrt(kernel) * TemperatureResiduals_matrix
emp_covariance1 <- emp_covariance1_temp %*% t(emp_covariance1_temp) / denominator

SalinityResiduals_matrix <- matrix(SalinityResiduals, length(SalinityResiduals), length(SalinityResiduals), byrow = T)

emp_covariance2_temp <- sqrt(kernel) * SalinityResiduals_matrix
emp_covariance2 <- emp_covariance2_temp %*% t(emp_covariance2_temp) / denominator

emp_covariance <- emp_covariance1_temp %*% t(emp_covariance2_temp) / denominator

plot(0, 0, type = 'n', xlim = range(diag(emp_covariance1)), ylim = c(2000, 0), xlab = '', ylab = '')
for(ll in 1:50){
  lines(diag(emp_covariance1[(ll - 1) * 50 + 1:50, (ll - 1) * 50 + 1:50]), loc3d[(ll - 1) * 50 + 1:50, 3])
}
mtext('Empirical Variance', side = 1, line = 2)
mtext('Depth (meters)', side = 2, line = 2)
mtext('Temperature Residuals', side = 3, line = 0)

plot(0, 0, type = 'n', xlim = range(diag(emp_covariance2)), ylim = c(2000, 0), xlab = '', ylab = '')
for(ll in 1:50){
  lines(diag(emp_covariance2[(ll - 1) * 50 + 1:50, (ll - 1) * 50 + 1:50]), loc3d[(ll - 1) * 50 + 1:50, 3])
}
mtext('Empirical Variance', side = 1, line = 2)
mtext('Depth (meters)', side = 2, line = 2)
mtext('Salinity Residuals', side = 3, line = 0)

plot(0, 0, type = 'n', xlim = c(-1, 1), ylim = c(2000, 0), xlab = '', ylab = '')
mtext('Empirical Colocated Correlation', side = 1, line = 2)
mtext('Depth (meters)', side = 2, line = 2)
mtext('Temperature & Salinity Residuals', side = 3, line = 0)
for(ll in 1:50){
  lines(diag(emp_covariance[(ll - 1) * 50 + 1:50, (ll - 1) * 50 + 1:50]) / sqrt(diag(emp_covariance1[(ll - 1) * 50 + 1:50, (ll - 1) * 50 + 1:50]) * diag(emp_covariance2[(ll - 1) * 50 + 1:50, (ll - 1) * 50 + 1:50])), loc3d[(ll - 1) * 50 + 1:50, 3])
}

#, ylim = c(2000, 0), xlab = "Covariance", ylab = "Lag in Depth (meters)"
############ STEP 0: WLS ############

ind_pred <- profile_no_closest
locs_insample <- loc3d[-ind_pred, ]
locs_outsample <- loc3d[ind_pred, ]

KNOTS1 <- KNOTS2 <- c(50, 100, 300, 500, 700, 1000)
#KNOTS1 <- KNOTS2 <- c(50, 100, 300, 500, 1000)

n <- 1000 #start with 100 to get initial values quickly
empirical_values <- rbind(cbind(emp_covariance1[1:n, 1:n], emp_covariance[1:n, 1:n]),
                          cbind(t(emp_covariance[1:n, 1:n]), emp_covariance2[1:n, 1:n]))

location = locs_insample[1:n, ]

basis1 <- bsplineBasis(locs_insample[, 3], 2, KNOTS1)
basis1 <- basis1[1:n, ]
nb1 <- ncol(basis1)
basis2 <- bsplineBasis(locs_insample[, 3], 2, KNOTS2)
basis2 <- basis2[1:n, ]
nb2 <- ncol(basis2)

set.seed(1235)
init1 <- runif(nb1, -10, 10)

set.seed(1236)
init2 <- runif(nb2, -1, 1)

theta = c(0.001, 0.001, init1, 0.001, 0.001, init2)
#theta = c(0.001, 0.001, 0.001, 0.001)
fit <- optim(par = theta, fn = NEGLOGLIK, empirical_values = empirical_values, basis1 = basis1, basis2 = basis2, radius = earthRadiusKm, control = list(trace = 5, maxit = 500))

for(aa in 1:100){
  fit <- optim(par = fit$par, fn = NEGLOGLIK, empirical_values = empirical_values, basis1 = basis1, basis2 = basis2, radius = earthRadiusKm, control = list(trace = 5, maxit = 500))
}

############ STEP 2: PLOTTING MLE ESTIMATES FROM SABINE ############

radius = earthRadiusKm

location = cbind(locs_insample[1, 1], locs_insample[1, 2], seq(0, 2000, length.out = 1000))
location = cbind(locs_insample[51, 1], locs_insample[51, 2], seq(0, 2000, length.out = 1000))

basis1 <- bsplineBasis(location[, 3], 2, KNOTS1)
basis2 <- bsplineBasis(location[, 3], 2, KNOTS2)

###---I3---###

BETA <- 0
SCALE_HORIZONTAL <- 0.019154
SCALE_VERTICAL <- 0.015625

A1 <- 1e-06
B1 <- 0.001001
C1_coef <- c(13.676718,25.385469,7.034342,152.88769,2.982689,-30.162745,16.766043,-18.108888,4.544768)

A2 <- -0.000234
B2 <- 2.6e-05
C2_coef <- c(5.510551,3.409904,3.606057,16.453042,-5.117168,1.281018,-1.849734,-0.29168,-3.258092)

C1_param <- basis1 %*% matrix(C1_coef, ncol = 1)
C2_param <- basis2 %*% matrix(C2_coef, ncol = 1)

cov_mat <- cov_bi_differential(location = location, beta = BETA,
                               scale_horizontal = SCALE_HORIZONTAL, scale_vertical = SCALE_VERTICAL,
                               a1 = A1, b1 = B1, c1 = C1_param, d1 = 0, a2 = A2, b2 = B2, c2 = C2_param, d2 = 0,
                               radius = radius)

variance1_I3 <- diag(cov_mat[1:nrow(location), 1:nrow(location)])
variance2_I3 <- diag(cov_mat[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)])
covariance12_I3 <- diag(cov_mat[1:nrow(location), nrow(location) + 1:nrow(location)])

###---B3---###

BETA <- 1
SCALE_HORIZONTAL <- 0.040733
SCALE_VERTICAL <- 0.012211

A1 <- 0.00152
B1 <- 0.002033
C1_coef <- c(31.35856,12.256453,38.2301,98.955665,1.066466,-30.909747,9.946264,-7.338746,1.37714)

A2 <- -0.000268
B2 <- -0.000128
C2_coef <- c(5.715104,2.190801,7.604064,15.965054,-2.48662,-0.602377,1.873185,-1.720387,-0.26431)

C1_param <- basis1 %*% matrix(C1_coef, ncol = 1)
C2_param <- basis2 %*% matrix(C2_coef, ncol = 1)

cov_mat <- cov_bi_differential(location = location, beta = BETA,
                               scale_horizontal = SCALE_HORIZONTAL, scale_vertical = SCALE_VERTICAL,
                               a1 = A1, b1 = B1, c1 = C1_param, d1 = 0, a2 = A2, b2 = B2, c2 = C2_param, d2 = 0,
                               radius = radius)

variance1_B3 <- diag(cov_mat[1:nrow(location), 1:nrow(location)])
variance2_B3 <- diag(cov_mat[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)])
covariance12_B3 <- diag(cov_mat[1:nrow(location), nrow(location) + 1:nrow(location)])

###---B4---###

BETA <- 0.566572
SCALE_HORIZONTAL <- 0.051631
SCALE_VERTICAL <- 0.007955

A1 <- 0.000222
B1 <- 2.8e-05
C1_coef <- c(69.520954,46.056877,169.712968,65.781759,60.943436,24.552186,2.421603,-0.69497,0)

A2 <- -4.5e-05
B2 <- -6e-06
C2_coef <- c(16.668018,8.699763,56.276284,7.009517,6.364778,1.74203,-1.214807,0.695601,-0.194455)

C1_param <- basis1 %*% matrix(C1_coef, ncol = 1)
C2_param <- basis2 %*% matrix(C2_coef, ncol = 1)

cov_mat <- cov_bi_differential(location = location, beta = BETA,
                               scale_horizontal = SCALE_HORIZONTAL, scale_vertical = SCALE_VERTICAL,
                               a1 = A1, b1 = B1, c1 = C1_param, d1 = 0, a2 = A2, b2 = B2, c2 = C2_param, d2 = 0,
                               radius = radius)

variance1_B4 <- diag(cov_mat[1:nrow(location), 1:nrow(location)])
variance2_B4 <- diag(cov_mat[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)])
covariance12_B4 <- diag(cov_mat[1:nrow(location), nrow(location) + 1:nrow(location)])


###---I1---###

variance1_I1 <- rep(0.0738046, nrow(location))
variance2_I1 <- rep(0.00350332, nrow(location))
covariance12_I1 <- rep(0, nrow(location))

###---B1---###

variance1_B1 <- rep(0.0713464, nrow(location))
variance2_B1 <- rep(0.003275922, nrow(location))
correlation12_B1 <- rep(0.254685, nrow(location))

###

MODEL_NAMES <- c("I1", "I2", "I3", "B1", "B2", "B3", "B4")
col_vec <- c("#38b5c3", "#8f00ff", "#FF00FF", "#0000bd", "#f6e34d", "#f7b329", "#e65013")

VARIABLE_TYPE = c('Temperature (°C)', 'Salinity (PSU)', 'Temperature Empirical Variance', 'Salinity Empirical Variance', 'Empirical Cross-Correlation')
LINE_TYPE  = c(1, 0, 2, 3, 0, 1, 5)

VARIABLE_NAME_LINE1 = c('Variance', 'Variance', 'Colocated Correlation')
VARIABLE_NAME_LINE2 = c('Temperature', 'Salinity', 'Temperature & Salinity')


pdf(file = paste('/Users/laisalvana/Library/CloudStorage/GoogleDrive-yourlainess@gmail.com/My Drive/Work/Research/UH/bivariate_argo/figures/argo_ref_loc', REF_LOC_IND, '_fitted_vs_empirical_covariance_new.pdf', sep = ''), width = 12, height = 5)

split.screen( rbind(c(0.09,0.99,0.13,0.95), c(0.98,0.99,0.07,0.95)))
split.screen( figs = c(1, 3), screen = 1)

screen(3)
par(mai = c(0.1, 0.1, 0.1, 0.25))

plot(0, 0, type = 'n', xlim = range(diag(emp_covariance1)), ylim = c(2000, 0), xlab = '', ylab = '')
for(ll in 1:50){
  lines(diag(emp_covariance1[(ll - 1) * 50 + 1:50, (ll - 1) * 50 + 1:50]), loc3d[(ll - 1) * 50 + 1:50, 3], lwd = 0.2, col = 'gray')
}

lines(variance1_I1, location[, 3], col = col_vec[1], lty = LINE_TYPE[1], lwd = 2)
lines(variance1_I3, location[, 3], col = col_vec[3], lty = LINE_TYPE[3], lwd = 2)
lines(variance1_B1, location[, 3], col = col_vec[4], lty = LINE_TYPE[4], lwd = 2)
lines(variance1_B3, location[, 3], col = col_vec[6], lty = LINE_TYPE[6], lwd = 2)
lines(variance1_B4, location[, 3], col = col_vec[7], lty = LINE_TYPE[7], lwd = 2)

abline(h = KNOTS1, col = "#964B00", lty = 3, lwd = 0.5)

mtext(VARIABLE_NAME_LINE1[1], side = 1, line = 2.5, cex = 1)
mtext('Depth (meters)', side = 2, line = 2)
mtext(VARIABLE_NAME_LINE2[1], side = 3, line = 0.5, cex = 1.2, col = '#336699', font = 2)
mtext(paste('Reference Location ', REF_LOC_IND, sep = ''), side = 2, line = 3, cex = 1.2, col = "#336699", font = 2)

legend(1.7, 1400, legend=c("I1", "I3", "B1", "B3", "B4"), col = col_vec[c(1, 3, 4, 6, 7)], lty = LINE_TYPE[c(1, 3, 4, 6, 7)], cex = 0.8, box.lty = 0, lwd = 2)

screen(4)
par(mai = c(0.1, 0.1, 0.1, 0.25))

plot(0, 0, type = 'n', xlim = range(diag(emp_covariance2)), ylim = c(2000, 0), xlab = '', ylab = '', yaxt = 'n')
for(ll in 1:50){
  lines(diag(emp_covariance2[(ll - 1) * 50 + 1:50, (ll - 1) * 50 + 1:50]), loc3d[(ll - 1) * 50 + 1:50, 3], lwd = 0.2, col = 'gray')
}

lines(variance2_I1, location[, 3], col = col_vec[1], lty = LINE_TYPE[1], lwd = 2)
lines(variance2_I3, location[, 3], col = col_vec[3], lty = LINE_TYPE[3], lwd = 2)
lines(variance2_B1, location[, 3], col = col_vec[4], lty = LINE_TYPE[4], lwd = 2)
lines(variance2_B3, location[, 3], col = col_vec[6], lty = LINE_TYPE[6], lwd = 2)
lines(variance2_B4, location[, 3], col = col_vec[7], lty = LINE_TYPE[7], lwd = 2)

abline(h = KNOTS2, col = "#964B00", lty = 3, lwd = 0.5)

mtext(VARIABLE_NAME_LINE1[2], side = 1, line = 2.5, cex = 1)
mtext(VARIABLE_NAME_LINE2[2], side = 3, line = 0.5, cex = 1.2, col = '#336699', font = 2)

legend(0.065, 1400, legend=c("I1", "I3", "B1", "B3", "B4"), col = col_vec[c(1, 3, 4, 6, 7)], lty = LINE_TYPE[c(1, 3, 4, 6, 7)], cex = 0.8, box.lty = 0, lwd = 2)

screen(5)
par(mai = c(0.1, 0.1, 0.1, 0.25))

plot(0, 0, type = 'n', xlim = c(-1, 1), ylim = c(2000, 0), xlab = '', ylab = '', yaxt = 'n')
for(ll in 1:50){
  lines(diag(emp_covariance[(ll - 1) * 50 + 1:50, (ll - 1) * 50 + 1:50]) / sqrt(diag(emp_covariance1[(ll - 1) * 50 + 1:50, (ll - 1) * 50 + 1:50]) * diag(emp_covariance2[(ll - 1) * 50 + 1:50, (ll - 1) * 50 + 1:50])), loc3d[(ll - 1) * 50 + 1:50, 3], lwd = 0.2, col = 'gray')
}

lines(covariance12_I1 / sqrt(variance1_I1 * variance2_I1), location[, 3], col = col_vec[1], lty = LINE_TYPE[1], lwd = 2)
lines(covariance12_I3 / sqrt(variance1_I3 * variance2_I3), location[, 3], col = col_vec[3], lty = LINE_TYPE[3], lwd = 2)
lines(correlation12_B1, location[, 3], col = col_vec[4], lty = LINE_TYPE[4], lwd = 2)
lines(covariance12_B3 / sqrt(variance1_B3 * variance2_B3), location[, 3], col = col_vec[6], lty = LINE_TYPE[6], lwd = 2)
lines(covariance12_B4 / sqrt(variance1_B4 * variance2_B4), location[, 3], col = col_vec[7], lty = LINE_TYPE[7], lwd = 2)

mtext(VARIABLE_NAME_LINE1[3], side = 1, line = 2.5, cex = 1)
mtext(VARIABLE_NAME_LINE2[3], side = 3, line = 0.5, cex = 1.2, col = '#336699', font = 2)

legend(0.5, 1400, legend=c("I1", "I3", "B1", "B3", "B4"), col = col_vec[c(1, 3, 4, 6, 7)], lty = LINE_TYPE[c(1, 3, 4, 6, 7)], cex = 0.8, box.lty = 0, lwd = 2)

close.screen( all=TRUE)
dev.off()

############ STEP 00: PLOTTING EMPIRICAL VARIANCES AND COLOCATED CORRELATIONS FOR JSS PAPER ############

VARIABLE_TYPE = c('Temperature (°C)', 'Salinity (PSU)', 'Temperature Empirical Variance', 'Salinity Empirical Variance', 'Empirical Cross-Correlation')
LINE_TYPE  = c(1, 0, 2, 3, 0, 1, 5)

VARIABLE_NAME_LINE1 = c('Variance', 'Variance', 'Colocated Correlation')
VARIABLE_NAME_LINE2 = c('Temperature', 'Salinity', 'Temperature & Salinity')

pdf(file = paste('/Users/laisalvana/Library/CloudStorage/GoogleDrive-yourlainess@gmail.com/My Drive/Work/Research/UH/bivariate_argo/figures/argo_ref_loc', REF_LOC_IND, '_empirical_variances_colocated_correlation_jss.pdf', sep = ''), width = 12, height = 5)

split.screen( rbind(c(0.09,0.99,0.13,0.95), c(0.98,0.99,0.07,0.95)))
split.screen( figs = c(1, 3), screen = 1)

screen(3)
par(mai = c(0.1, 0.1, 0.1, 0.25))

plot(0, 0, type = 'n', xlim = range(diag(emp_covariance1)), ylim = c(2000, 0), xlab = '', ylab = '')
for(ll in 1:50){
  lines(diag(emp_covariance1[(ll - 1) * 50 + 1:50, (ll - 1) * 50 + 1:50]), loc3d[(ll - 1) * 50 + 1:50, 3], lwd = 0.2, col = 'gray')
}

mtext(VARIABLE_NAME_LINE1[1], side = 1, line = 2.5, cex = 1)
mtext('Depth (meters)', side = 2, line = 2)
mtext(VARIABLE_NAME_LINE2[1], side = 3, line = 0.5, cex = 1.2, col = '#336699', font = 2)
mtext(paste('Reference Location ', REF_LOC_IND, sep = ''), side = 2, line = 3, cex = 1.2, col = "#336699", font = 2)

screen(4)
par(mai = c(0.1, 0.1, 0.1, 0.25))

plot(0, 0, type = 'n', xlim = range(diag(emp_covariance2)), ylim = c(2000, 0), xlab = '', ylab = '', yaxt = 'n')
for(ll in 1:50){
  lines(diag(emp_covariance2[(ll - 1) * 50 + 1:50, (ll - 1) * 50 + 1:50]), loc3d[(ll - 1) * 50 + 1:50, 3], lwd = 0.2, col = 'gray')
}

mtext(VARIABLE_NAME_LINE1[2], side = 1, line = 2.5, cex = 1)
mtext(VARIABLE_NAME_LINE2[2], side = 3, line = 0.5, cex = 1.2, col = '#336699', font = 2)

screen(5)
par(mai = c(0.1, 0.1, 0.1, 0.25))

plot(0, 0, type = 'n', xlim = c(-1, 1), ylim = c(2000, 0), xlab = '', ylab = '', yaxt = 'n')
for(ll in 1:50){
  lines(diag(emp_covariance[(ll - 1) * 50 + 1:50, (ll - 1) * 50 + 1:50]) / sqrt(diag(emp_covariance1[(ll - 1) * 50 + 1:50, (ll - 1) * 50 + 1:50]) * diag(emp_covariance2[(ll - 1) * 50 + 1:50, (ll - 1) * 50 + 1:50])), loc3d[(ll - 1) * 50 + 1:50, 3], lwd = 0.2, col = 'gray')
}

mtext(VARIABLE_NAME_LINE1[3], side = 1, line = 2.5, cex = 1)
mtext(VARIABLE_NAME_LINE2[3], side = 3, line = 0.5, cex = 1.2, col = '#336699', font = 2)

close.screen( all=TRUE)
dev.off()

############ STEP 01: PLOTTING EMPIRICAL COVARIANCES AND CROSS-COVARIANCES FOR JSS PAPER ############

VARIABLE_TYPE = c('Temperature (°C)', 'Salinity (PSU)', 'Temperature Empirical Variance', 'Salinity Empirical Variance', 'Empirical Cross-Correlation')
LINE_TYPE  = c(1, 0, 2, 3, 0, 1, 5)

VARIABLE_NAME_LINE1 = c('Correlation', 'Correlation', 'Cross-Correlation')
VARIABLE_NAME_LINE2 = c('Temperature', 'Salinity', 'Temperature & Salinity')

pdf(file = paste('/Users/laisalvana/Library/CloudStorage/GoogleDrive-yourlainess@gmail.com/My Drive/Work/Research/UH/bivariate_argo/figures/argo_ref_loc', REF_LOC_IND, '_empirical_covariances_jss.pdf', sep = ''), width = 12, height = 5)

split.screen( rbind(c(0.09,0.99,0.13,0.95), c(0.98,0.99,0.07,0.95)))
split.screen( figs = c(1, 3), screen = 1)

screen(3)
par(mai = c(0.1, 0.1, 0.1, 0.25))

plot(0, 0, type = 'n', xlim = c(0, 1), ylim = c(2000, 0), xlab = '', ylab = '')
for(ll in 1:50){
  cov_sub <- emp_covariance1[(ll - 1) * 50 + 1, (ll - 1) * 50 + 1:50] #/ emp_covariance1[ll, ll]
  locs_sub <- loc3d[(ll - 1) * 50 + 1:50, 3] - loc3d[(ll - 1) * 50 + 1, 3]
  ind_sub <- which(locs_sub >= 0)
  lines(cov_sub[ind_sub], locs_sub[ind_sub])
}
lines(cov_mat[1, 1:1000], location[, 3], col = col_vec[1], lty = LINE_TYPE[1], lwd = 2)


mtext(VARIABLE_NAME_LINE1[1], side = 1, line = 2.5, cex = 1)
mtext('Lag in Depth (meters)', side = 2, line = 2)
mtext(VARIABLE_NAME_LINE2[1], side = 3, line = 0.5, cex = 1.2, col = '#336699', font = 2)
mtext(paste('Reference Location ', REF_LOC_IND, sep = ''), side = 2, line = 3, cex = 1.2, col = "#336699", font = 2)

screen(4)
par(mai = c(0.1, 0.1, 0.1, 0.25))

plot(0, 0, type = 'n', xlim = c(0, 1.5), ylim = c(2000, 0), xlab = '', ylab = '', yaxt = 'n')
for(ll in 1:50){
  lines(emp_covariance2[(ll - 1) * 50 + 1, (ll - 1) * 50 + 1:50] / emp_covariance2[(ll - 1) * 50 + 1, (ll - 1) * 50 + 1], loc3d[(ll - 1) * 50 + 1:50, 3] - loc3d[(ll - 1) * 50 + 1, 3], lwd = 0.2, col = 'gray')
}

mtext(VARIABLE_NAME_LINE1[2], side = 1, line = 2.5, cex = 1)
mtext(VARIABLE_NAME_LINE2[2], side = 3, line = 0.5, cex = 1.2, col = '#336699', font = 2)

screen(5)
par(mai = c(0.1, 0.1, 0.1, 0.25))

plot(0, 0, type = 'n', xlim = c(-1, 1), ylim = c(2000, 0), xlab = '', ylab = '', yaxt = 'n')
for(ll in 1:50){
  lines(emp_covariance[(ll - 1) * 50 + 1, (ll - 1) * 50 + 1:50] / sqrt(emp_covariance1[(ll - 1) * 50 + 1, (ll - 1) * 50 + 1] * emp_covariance2[(ll - 1) * 50 + 1, (ll - 1) * 50 + 1]), loc3d[(ll - 1) * 50 + 1:50, 3] - loc3d[(ll - 1) * 50 + 1, 3], lwd = 0.2, col = 'gray')
}

mtext(VARIABLE_NAME_LINE1[3], side = 1, line = 2.5, cex = 1)
mtext(VARIABLE_NAME_LINE2[3], side = 3, line = 0.5, cex = 1.2, col = '#336699', font = 2)

close.screen( all=TRUE)
dev.off()

##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################

Z_insample <- c(TemperatureResiduals_insample, SalinityResiduals_insample)
Z_outsample <- c(TemperatureResiduals_outsample, SalinityResiduals_outsample)


INIT_BETA = 0
INIT_A1 = INIT_B1 = 0
INIT_A2 = INIT_B2 = 0
SCALE_HORIZONTAL = 0.001
SCALE_VERTICAL = 0.01


KNOTS1 <- seq(0, 2100, length.out = 2)
KNOTS2 <- seq(0, 2100, length.out = 2)

basis1 <- bsplineBasis(locs_insample[, 3], 0, KNOTS1)
nb1 <- ncol(basis1)
basis2 <- bsplineBasis(locs_insample[, 3], 0, KNOTS2)
nb2 <- ncol(basis2)

INIT_C1 <- 10
INIT_C2 <- 10

est_params <- est_bi_differential(residuals = matrix(Z_insample, ncol = 1), location = locs_insample, init_beta = INIT_BETA,
                                  init_scale_horizontal = log(SCALE_HORIZONTAL), init_scale_vertical = log(SCALE_VERTICAL),
                                  init_a1 = INIT_A1, init_b1 = INIT_B1, init_c1 = INIT_C1, init_a2 = INIT_A2, init_b2 = INIT_B2, init_c2 = INIT_C2,
                                  basis1 = basis1, basis2 = basis2, radius = earthRadiusKm, splines_degree = 0, knots1 = KNOTS1, knots2 = KNOTS2, STEPMAX = 10)

#xxxxxxxxxxxxxxxxxxxxxxx TRYING OUT TO SEE WHAT GOOD INITIAL VALUES TO USE

theta = c(0.185367,log(0.03397),log(0.03397),6e-06,8e-06,0.001293,2e-06,-3e-06,0.000359)

location = locs_insample

location <- cbind(-174.983, 39.919, seq(0, 2000, length.out = 100)/2000)

BETA <- 0.185367
SCALE_HORIZONTAL <- 0.05881
SCALE_VERTICAL <- 22
A1 <- 9e-06
B1 <- 1.5e-05
A2 <- -1e-06
B2 <- 1.4e-05
C1 <- 0.0009
C2 <- 0.0001
D1 <- 0
D2 <- 0

cov_mat <- cov_bi_differential(location = location, beta = BETA,
                               scale_horizontal = SCALE_HORIZONTAL, scale_vertical = SCALE_VERTICAL,
                               a1 = A1, b1 = B1, c1 = C1, d1 = D1, a2 = A2, b2 = B2, c2 = C2, d2 = D2,
                               radius = earthRadiusKm)

plot(location[1:100, 3], cov_mat[1, 1:100], type = 'l', lwd = 2)

#xxxxxxxxxxxxxxxxxxxxxxx

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

INIT_BETA = 0
INIT_SCALE_HORIZONTAL = 0.0011
INIT_SCALE_VERTICAL = 20
INIT_A1 = 5.7e-05
INIT_A2 = 1.1e-05
INIT_B1 = -3.5e-05
INIT_B2 = 1e-06
INIT_D1 = INIT_D2 = 0

KNOTS1 <- c(0.1, 0.5)
KNOTS2 <- c(0.1, 0.5)

#basis1 <- bsplineBasis(LOCATION[, 3] / 2000, 4, KNOTS1)
basis1 <- bsplineBasis(locs_insample[, 3], 2, KNOTS1)
basis1 <- basis1[1:50, ]
nb1 <- ncol(basis1)
#basis2 <- bsplineBasis(LOCATION[, 3] / 2000, 4, KNOTS2)
basis2 <- bsplineBasis(locs_insample[, 3], 2, KNOTS2)
basis2 <- basis2[1:50, ]
nb2 <- ncol(basis2)

set.seed(1235)
INIT_C1 <- runif(nb1, -0, 0)

set.seed(1236)
INIT_C2 <- runif(nb2, -0, 0)

est_params <- est_bi_differential_wls(empirical_values = EMPIRICAL_VALUES, location = locs_insample[1:50, ], init_beta = INIT_BETA,
                                      init_scale_horizontal = INIT_SCALE_HORIZONTAL, init_scale_vertical = INIT_SCALE_VERTICAL, init_scale_horizontal_fix = T, init_scale_vertical_fix = T,
                                      init_a1 = INIT_A1, init_b1 = INIT_B1, init_c1 = INIT_C1, init_a2 = INIT_A2, init_b2 = INIT_B2, init_c2 = INIT_C2, init_beta_fix = F, init_a1_fix = T, init_b1_fix = T, init_a2_fix = T, init_b2_fix = T,
                                      basis1 = basis1, basis2 = basis2, radius = earthRadiusKm, splines_degree = 2, knots1 = KNOTS1, knots2 = KNOTS2)

#RERUN FOR WLS
theta = c(0.68450632,0.02050917,0.01830979,-0.00447512,0.00185042,0.00331421,0.00014238,-0.00039106,-5.882e-05)

INIT_A1 = theta[1]
INIT_A2 = theta[nb1 + 3]
INIT_B1 = theta[2]
INIT_B2 = theta[nb1 + 4]
INIT_C1 <- theta[2 + 1:nb1]
INIT_C2 <- theta[nb1 + 4 + 1:nb2]

############ STEP 2: MLE ############

basis1 <- bsplineBasis(locs_insample[, 3], 2, KNOTS1)
nb1 <- ncol(basis1)
basis2 <- bsplineBasis(locs_insample[, 3], 2, KNOTS2)
nb2 <- ncol(basis2)

#initial parameters: n = 2500, including the prediction locations
theta = c(0.7, log(INIT_SCALE_HORIZONTAL), log(INIT_SCALE_VERTICAL), 0,0.000137,-43.478744,106.037114,-33.149147,75.450859,58.867069,-11.714446,8.646101,2e-06,-1.9e-05,-5.850154,14.24039,-6.484243,-5.210463,1.811206,-0.15922,5.527517)

#initial parameters: n = 2500, excluding the prediction locations
theta = c(INIT_BETA, log(INIT_SCALE_HORIZONTAL), log(INIT_SCALE_VERTICAL), 2.2e-07,6.804e-05,-51.52184895,71.74782542,-66.36078398,81.95555724,69.22938586,-10.68999432,6.89582243,2.46e-06,-3.811e-05,-5.82206381,13.2819233,-3.61022406,-9.577206,-3.84409763,1.44576138,8.03768149)

#initial parameters: n = 2500, excluding the prediction locations with additional reruns for WLS
theta = c(INIT_BETA, log(INIT_SCALE_HORIZONTAL), log(INIT_SCALE_VERTICAL), 5.13e-06,2.42e-05,-52.73111221,71.88003595,-62.47084003,83.34465397,67.75742297,-12.37265794,3.18738643,-1.54e-06,-1.543e-05,-4.87730718,12.34853351,-7.01972981,-6.86201992,7.72033671,-5.78460723,8.77574025)

#initial parameters: n = 2500, excluding the prediction locations, splines_degree = 2
theta = c(INIT_BETA, log(INIT_SCALE_HORIZONTAL), log(INIT_SCALE_VERTICAL), 3.944e-05,3.202e-05,-43.20737327,71.17479335,-47.12204867,52.61168296,39.0776386,-13.49731551,-6.89e-06,-5.71e-06,3.15369315,7.23586132,-8.81128232,4.81225399,-2.5139546,13.39680324)

#initial parameters: n = 2500, WLS until the end
theta = c(0.81629834, log(INIT_SCALE_HORIZONTAL), log(INIT_SCALE_VERTICAL),0.01751006,0.05029142,114.49443241,133.77856476,-38.71210424,-4.28989463,6.92862956,16.86219532,-0.00070489,-0.00133181,24.98619127,0.80020364,0.2593444,-0.00564173,-4.51826103,-2.75454255)

est_params <- est_bi_differential(residuals = c(Z_insample), location = locs_insample, init_beta = theta[1],
                                  init_scale_horizontal = theta[2], init_scale_vertical = theta[3],
                                  init_a1 = theta[4], init_b1 = theta[5], init_c1 = theta[5 + 1:nb1], init_a2 = theta[nb1 + 6], init_b2 = theta[nb1 + 7], init_c2 = theta[nb1 + 7 + 1:nb2],
                                  basis1 = basis1, basis2 = basis2, radius = earthRadiusKm, splines_degree = 2, knots1 = KNOTS1, knots2 = KNOTS2, STEPMAX = 500, MAXIT = 2000, RERUNS = 0)

#MLE: n = 2500
theta = c(0.72545012,-1.64953152,-4.05563362,2.53e-06,0.00021368,-52.73222083,76.33190581,-62.47251572,83.34264221,67.75549558,-12.3718472,3.18806276,-1.366e-05,-1.92e-06,-4.87616132,12.34728918,-6.4992468,-3.60197678,7.72209309,-5.78397767,8.7763426)
#neg loglik = -7017.1351

fisher_info<-solve(-est_params)

#xxxxxxxxxxxxxxxxxxxxxxx VISUALIZING IMPLIED COVARIANCE STRUCTURE

theta = c(0.68450632, log(INIT_SCALE_HORIZONTAL), log(INIT_SCALE_VERTICAL),INIT_A1, INIT_B1,0.02050917,0.01830979,-0.00447512,0.1, INIT_A2, INIT_B2,0.00331421,0.00014238,-0.00039106,-5.882e-05)

C1_coef <- theta[5 + 1:nb1]
C2_coef <- theta[nb1 + 7 + 1:nb2]

C1 <- basis1 %*% matrix(C1_coef, ncol = 1)
C2 <- basis2 %*% matrix(C2_coef, ncol = 1)

cov_mat <- cov_bi_differential(location = locs_insample, beta = theta[1],
                               scale_horizontal = exp(theta[2]), scale_vertical = exp(theta[3]),
                               a1 = theta[4], b1 = theta[5], c1 = C1, d1 = 0, a2 = theta[nb1 + 6], b2 = theta[nb1 + 7], c2 = C2, d2 = 0,
                               radius = earthRadiusKm)

pdf(paste('/Users/laisalvana/Library/CloudStorage/GoogleDrive-yourlainess@gmail.com/My Drive/Work/Research/UH/bivariate_argo/figures/estimated_covariance_contours.pdf', sep = ''), width = 15, height = 5)

split.screen( rbind(c(0.03,0.99,0.06,0.99), c(0.98,0.99,0.06,0.99)))
split.screen( figs = c(1, 3), screen = 1)

screen(3)
par(pty = 's')
par(mai = c(0.4, 0.4, 0.4, 0.4))

plot(cov_mat[1, 1:50], locs_insample[1:50, 3], ylab = '', xlab = '', ylim = c(max(locs_insample[1:50, 3]), 0), type = 'l', lwd = 2)
mtext(expression(hat(C)[11]), side = 1, line = 2.5, cex = 1.2)
mtext('Vertical Distance (meters)', side = 2, line = 2)

screen(4)
par(pty = 's')
par(mai = c(0.4, 0.4, 0.4, 0.4))

plot(cov_mat[2450 + 1, 2450 + 1:50], locs_insample[1:50, 3], ylab = '', xlab = '', ylim = c(2000, 0), type = 'l', lwd = 2)
mtext(expression(hat(C)[22]), side = 1, line = 2.5, cex = 1.2)
mtext('Vertical Distance (meters)', side = 2, line = 2)

screen(5)
par(pty = 's')
par(mai = c(0.4, 0.4, 0.4, 0.4))

plot(cov_mat[1, 2450 + 1:50], locs_insample[1:50, 3], ylab = '', xlab = '', ylim = c(2000, 0), type = 'l', lwd = 2)
mtext(expression(hat(C)[12]), side = 1, line = 2.5, cex = 1.2)
mtext('Vertical Distance (meters)', side = 2, line = 2)

close.screen( all=TRUE)
dev.off()

pdf(paste('/Users/laisalvana/Library/CloudStorage/GoogleDrive-yourlainess@gmail.com/My Drive/Work/Research/UH/bivariate_argo/figures/estimated_variance_correlation_curves.pdf', sep = ''), width = 15, height = 5)

split.screen( rbind(c(0.03,0.99,0.06,0.99), c(0.98,0.99,0.06,0.99)))
split.screen( figs = c(1, 3), screen = 1)

screen(3)
par(pty = 's')
par(mai = c(0.4, 0.4, 0.4, 0.4))

plot(diag(cov_mat[1:50, 1:50]), locs_insample[1:50, 3], ylab = '', xlab = '', ylim = c(max(locs_insample[1:50, 3]), 0), type = 'l', lwd = 2)
mtext('Estimated Variance', side = 1, line = 2.5, cex = 1.2)
mtext('Depth (meters)', side = 2, line = 2)
mtext('Temperature', side = 3, line = 1, cex = 1.2)

screen(4)
par(pty = 's')
par(mai = c(0.4, 0.4, 0.4, 0.4))

plot(diag(cov_mat[2450 + 1:50, 2450 + 1:50]), locs_insample[1:50, 3], ylab = '', xlab = '', ylim = c(2000, 0), type = 'l', lwd = 2)
mtext('Estimated Variance', side = 1, line = 2.5, cex = 1.2)
mtext('Depth (meters)', side = 2, line = 2)
mtext('Salinity', side = 3, line = 1, cex = 1.2)

screen(5)
par(pty = 's')
par(mai = c(0.4, 0.4, 0.4, 0.4))

plot(diag(cov_mat[1:50, 2450 + 1:50]), locs_insample[1:50, 3], ylab = '', xlab = '', ylim = c(2000, 0), type = 'l', lwd = 2)
mtext('Estimated Colocated Correlation', side = 1, line = 2.5, cex = 1.2)
mtext('Depth (meters)', side = 2, line = 2)
mtext('Temperature & Salinity', side = 3, line = 1, cex = 1.2)

close.screen( all=TRUE)
dev.off()

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


########################## NEW WLS ##########################

#KNOTS1 <- KNOTS2 <- c(seq(50, 500, by = 50), 700, 1000, 1300, 1500)
KNOTS1 <- KNOTS2 <- c(50, 100, 300, 500, 700, 1000)

n <- 1500
empirical_values <- rbind(cbind(emp_covariance1[1:n, 1:n], emp_covariance[1:n, 1:n]),
                          cbind(t(emp_covariance[1:n, 1:n]), emp_covariance2[1:n, 1:n]))

location = locs_insample[1:n, ]

basis1 <- bsplineBasis(locs_insample[, 3], 2, KNOTS1)
basis1 <- basis1[1:n, ]
nb1 <- ncol(basis1)
basis2 <- bsplineBasis(locs_insample[, 3], 2, KNOTS2)
basis2 <- basis2[1:n, ]
nb2 <- ncol(basis2)

set.seed(1235)
init1 <- runif(nb1, -10, 10)

set.seed(1236)
init2 <- runif(nb2, -1, 1)

theta = c(0.001, 0.001, init1, 0.001, 0.001, init2)
fit <- optim(par = theta, fn = NEGLOGLIK, empirical_values = empirical_values, basis1 = basis1, basis2 = basis2, radius = earthRadiusKm, control = list(trace = 5, maxit = 500))

for(aa in 1:100){
  fit <- optim(par = fit$par, fn = NEGLOGLIK, empirical_values = empirical_values, basis1 = basis1, basis2 = basis2, radius = earthRadiusKm, control = list(trace = 5, maxit = 500))
}

#n = 100, beta = 0.9 as we decrease SCALE_VERTICAL
theta = c(0.00246273,0.00030454,8.691133,13.66296296,10.48374036,9.85908405,-3.07955378,-1.74557456,0.47835926,-1.83597167,0.42006333,-0.00022336,-0.000179,1.03262241,1.54002773,1.28793131,0.97721898,-0.74341365,0.57893729,-0.25085098,1.147632,0.03080385)
theta = c(0.00182011,0.00213812,26.84453859,39.70320758,33.51312336,26.28802037,-3.98823905,-3.61716544,-6.07273884,-4.66245635,-2.78249983,-0.00012476,-5.842e-05,3.47459305,4.63655288,3.88476418,3.23205214,-2.8582196,1.18835578,0.23935119,-0.20851424,-0.11917048)
theta = c(0.00127428,0.00204556,78.77136832,126.46469432,94.44429823,88.48827232,-18.59373985,-18.88540596,-12.65328204,-0.88318001,-7.65037147,2.25e-06,-0.00015574,9.40701963,14.75755799,11.6663534,9.717654,-9.09087822,9.24741368,0.38303944,3.22505272,2.37411392)
theta = c(0.0024502,0.00181957,112.49654321,171.98004576,139.93266558,116.93197958,-20.05766198,-15.79292226,-21.19784691,-26.01094529,-31.22777421,-0.00015747,-5.108e-05,11.72314134,21.81409625,16.02939546,15.28234516,-12.82949102,4.91863381,0.67171598,-1.2906754,-0.08738533)
theta = c(0.00214882,0.00193461,165.13073313,244.77546622,188.59696616,176.58852246,-37.53694197,-7.33643604,-28.95033678,10.70904738,-42.20482341,-0.000213,-3.928e-05,20.43317613,30.73726563,22.24734922,22.32307169,-19.71249676,7.00241988,0.70604981,-3.83749137,0.33721687)
theta = c(0.00251538,0.00177685,260.88703733,415.61611505,319.05816469,280.80096987,-61.67379775,-44.76183115,7.18980369,-92.67025642,-34.89432237,-0.00017041,-0.00013596,29.99133117,50.28258277,34.57602393,41.39969698,-34.02841533,16.719796,-5.7028914,0.94306126,-3.76563713)
theta = c(0.00254991,0.00077242,812.69761809,1252.09205045,957.74967706,855.93283776,-77.41254158,-219.84547406,-98.93929691,-131.75604803,-113.23806571,-6.108e-05,-5.612e-05,93.94272737,147.15388304,116.55006739,100.39741991,-98.85023546,103.29105938,-15.80650128,81.73391265,22.66374368)
theta = c(0.0028038,0.00181274,1141.84614839,1692.10037884,1355.04761435,1163.13246018,-107.01408024,-213.13409006,-24.46851992,78.13540981,-313.63341313,-7.089e-05,-0.00014412,120.15051145,220.40940843,162.11958937,154.57111856,-142.89150778,60.52269326,-11.07401855,26.10198093,57.03272586)
theta = c(0.00272329,0.00162622,1222.67490108,2031.68519004,1610.37906853,1361.31739128,-194.74720104,-237.61843783,-297.95100644,-140.95825854,-381.76481599,-5.828e-05,-0.00012204,159.76973339,245.00973937,180.09962557,177.13344636,-163.12696862,78.99675592,6.4866513,138.05352629,5.17238747)

#n = 100, beta = 1, THIS IS WHAT IS IN SABINE
theta = c(0.00270846,0.00134791,1407.40731682,2045.95155486,1608.95637922,1426.34770783,-298.3636312,5.57373751,-52.53823064,-55.47290071,28.49712777,-0.00011204,-0.00016021,152.20580117,253.44214037,184.14986691,181.00982675,-158.709912,44.81169714,-1.41371192,84.88382948,12.47930132)

#n = 500, beta = 0.9
theta = c(0.00377438,0.0011499,1127.66382777,1451.93160997,1172.08337279,1674.03651586,-1264.63115544,240.06657449,-7.67761641,84.77458857,43.11874137,-8.324e-05,-0.00025795,240.11759237,297.35862948,274.29788798,265.36145771,-212.96484906,64.19380948,-9.6715385,38.74047081,-7.83948635)

#n = 500, beta = 1
theta = c(0.004,0.0002319,1161.2445868,1432.45853915,1172.36542887,1641.23712635,-1207.40945716,170.65686755,244.44896655,-98.22303306,91.27300748,-0.00015615,-0.0003584,239.30700447,286.78269202,272.00270608,256.63874006,-204.45616822,65.23760888,-45.97313344,4.75203489,-0.40076906)

#MLE results for plotting:
theta = c(6.9e-05, -4.6e-05,752.761599,667.342355,1505.32545,2676.630244,-1986.208588,-241.954419,-99.565546,-23.200609,-20.132394,
          -3.8e-05,1e-06,133.801299,123.670912,437.70727,618.161216,-177.746145,5.282143,22.495273,-0.614592,-0.199912)

#for estimating only spline coefficients
theta = c(461.077506,355.016643,786.448522,2634.159983,-1941.119526,-405.988123,-61.347585,-312.278314,-61.169192,
          135.739637,124.168438,350.327986,519.558449,-322.510729,25.809393,-2.178382,63.809757,8.475811)

##########################################
est_params <- est_bi_differential_wls(empirical_values = EMPIRICAL_VALUES, , init_beta = INIT_BETA,
                                      init_scale_horizontal = INIT_SCALE_HORIZONTAL, init_scale_vertical = INIT_SCALE_VERTICAL, init_scale_horizontal_fix = F, init_scale_vertical_fix = F,
                                      init_a1 = INIT_A1, init_b1 = INIT_B1, init_c1 = INIT_C1, init_a2 = INIT_A2, init_b2 = INIT_B2, init_c2 = INIT_C2, init_beta_fix = F, init_a1_fix = F, init_b1_fix = F, init_a2_fix = F, init_b2_fix = F,
                                      basis1 = basis1, basis2 = basis2, radius = earthRadiusKm, splines_degree = 2, knots1 = KNOTS1, knots2 = KNOTS2)

EMPIRICAL_VALUES <- rbind(cbind(emp_covariance1[1:1000, 1:1000], emp_covariance[1:1000, 1:1000]),
                          cbind(t(emp_covariance[1:1000, 1:1000]), emp_covariance2[1:1000, 1:1000]))
basis1 <- bsplineBasis(locs_insample[, 3], 2, KNOTS1)
basis1 <- basis1[1:1000, ]
nb1 <- ncol(basis1)
basis2 <- bsplineBasis(locs_insample[, 3], 2, KNOTS2)
basis2 <- basis2[1:1000, ]
nb2 <- ncol(basis2)

#WLS for n = 1000
theta = c(0.00062188,0.00126067,-114.30383942,-223.29745228,-167.39299729,-125.105902,-148.49184518,-95.75156444,-92.89504895,42.58379899,-20.31261471,23.63371261,5.31404379,-7.23588144,1.06748394,5.88909492,18.52800183,5.1496904,2.39408807)
theta = c(0.00105905,0.00266644,-95.68125219,-146.91209989,-99.91469332,-122.8035243,-105.68877791,-97.0385092,-85.28740214,32.30971907,-96.80087517,-12.4946009,-30.31868798,-0.50041972,2.37400337,14.53954995,7.86831353,-11.31296512,25.54145192)

theta = c(0.001059,0.002666, -102.541309,-171.51506,-156.5047,-55.70442,-331.886245,-215.424259,-411.777048,-481.514641,-247.968563,-32.336227,117.68406,306.137219,93.028573,-0.987142,-9.014911,-24.179167,18.053546)
theta = c(0.001059,0.002666,-310.423698,-117.777033,-193.406835,-457.93796,-323.861474,-275.45146,-492.230912,-638.747933,-310.492838,-112.947322,-69.183244,217.473221,27.830731,-159.609817,-15.934333,50.140593,45.09418)
theta = c(0.001059,0.002666,-323.905445,-138.357572,-168.267418,-534.679427,-313.563635,-293.671233,-526.719469,-733.746404,-350.953381,-197.20832,-112.687272,52.274697,4.340953,-42.493677,15.503204,9.176326,56.519844)
theta = c(0.001059,0.002666,-322.883178,-135.127953,-164.805608,-533.032856,-312.309966,-292.231291,-525.867973,-733.127661,-349.918948,-196.09496,-110.219758,44.156148,-3.278183,-43.605837,12.390492,4.797252,52.448918)
theta = c(0.001059,0.002666,-311.484984,-92.311772,-114.807185,-515.000313,-294.971651,-273.848417,-516.09628,-725.893603,-337.345239,-183.195167,-74.27393,-3.974978,-12.758201,1.218279,-14.169452,-43.892657,5.550889)
theta = c(0.001059,0.002666,-286.851665,-68.560686,-90.224611,-484.881974,-267.694716,-247.442822,-501.97557,-716.719536,-317.473624,-152.392828,-19.274097,-4.290325,-15.891378,-7.718991,13.648211,-5.751732,10.232644)
theta = c(0.001059,0.002666,-251.750153,-59.927407,-83.459851,-442.36852,-228.401144,-210.374896,-483.153859,-705.172936,-287.893993,-103.596897,-1.743336,-1.557021,-12.827538,-13.954925,5.406765,9.189622,7.812447)
theta = c(0.001059,0.002666,-193.559365,-82.031813,-46.235102,-368.816387,-169.261625,-157.963209,-450.367193,-685.182512,-234.791096,-57.799936,-5.834821,1.954769,-3.890077,-15.307471,-14.836148,38.721619,4.613329)
theta = c(0.001059,0.002666,-29.479818,-25.544409,-33.229045,-52.105533,-89.984239,-130.651503,-135.041047,-110.833017,-53.938547,-37.734964,-7.589682,4.194742,-4.520605,-16.082078,14.399544,-0.950459,-70.797023)
theta = c(0.001059,0.002666,-42.759396,-39.471443,-66.519635,-76.081034,-119.21612,-154.701493,-172.268414,-189.288076,-156.979219,-34.745888,6.539258,35.383086,-9.958257,-78.889858,-1.330844,0.619098,1.688505)
theta = c(3.6e-05,-3e-06,-28.303136,-25.675351,-33.132617,-52.177853,-84.153861,-114.907509,-131.597203,-108.055736,-58.82482,-38.054529,-10.097285,5.304196,-3.056047,-13.941319,-1.623208,-4.933468,-2.373866)

est_params <- est_bi_differential_wls(empirical_values = EMPIRICAL_VALUES, location = locs_insample[1:1000, ], init_beta = theta[1],
                                      init_scale_horizontal = INIT_SCALE_HORIZONTAL, init_scale_vertical = INIT_SCALE_VERTICAL, init_scale_horizontal_fix = F, init_scale_vertical_fix = F,
                                      init_a1 = theta[1], init_b1 = theta[2], init_c1 = theta[2 + 1:nb1], init_a2 = 0, init_b2 = 0, init_c2 = INIT_C2, init_beta_fix = F, init_a1_fix = F, init_b1_fix = F, init_a2_fix = F, init_b2_fix = F,
                                      basis1 = basis1, basis2 = basis2, radius = earthRadiusKm, splines_degree = 2, knots1 = KNOTS1, knots2 = KNOTS2)


#best: 40209
