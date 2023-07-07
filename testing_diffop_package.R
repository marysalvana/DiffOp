library(DiffOp)

####### SIMULATION #######

library(dplyr)

x <- seq(0, 1, length.out = 10)
y <- seq(0, 1, length.out = 10)
loc2d <- expand.grid(x, y) %>% as.matrix()

depth <- seq(0, 1, length.out = 10)
loc3d <- cbind(rep(loc2d[, 1], each = length(depth)), rep(loc2d[, 2], each = length(depth)), depth)

earthRadiusKm = 6371

BETA = 0.5
SCALE_HORIZONTAL = 0.03
SCALE_VERTICAL = 0.3
A1 = A2 = 0.00001
B1 = B2 = 0.00001
C1 = sin((loc3d[, 3] + 0.1) * pi / 0.5)
C2 = cos((loc3d[, 3] + 0.1) * pi / 0.5)
D1 = D2 = 0

cov_mat <- cov_bi_differential(location = loc3d, beta = BETA, scale_horizontal = SCALE_HORIZONTAL, scale_vertical = SCALE_VERTICAL, a1 = A1, b1 = B1, c1 = C1, d1 = D1, a2 = A2, b2 = B2, c2 = C2, d2 = D2, radius = earthRadiusKm)

library(MASS)

set.seed(1234)
Z <- mvrnorm(1, mu = rep(0, ncol(cov_mat)), Sigma = cov_mat)

Z1 <- Z[1:nrow(loc3d)]
Z2 <- Z[nrow(loc3d) + 1:nrow(loc3d)]

subset <- which(loc3d[, 3] == 0)

library(fields)

par(mfrow = c(2, 1))
par(pty = 's')
par(mai = c(0.3, 0.1, 0.6, 0.5))
quilt.plot(loc3d[subset, 1], loc3d[subset, 2], Z1[subset], nx = length(x), ny = length(y), xlab = "Longitude", ylab = "", zlim = range(Z), xaxt = 'n')
mtext("Latitude", side = 2, line = 2)
mtext(expression(Z[1]), side = 2, line = 3, col = 4, cex = 1.5)
par(mai = c(0.8, 0.1, 0.1, 0.5))
quilt.plot(loc3d[subset, 1], loc3d[subset, 2], Z2[subset], nx = length(x), ny = length(y), xlab = "Longitude", ylab = "", zlim = range(Z))
mtext("Latitude", side = 2, line = 2)
mtext(expression(Z[2]), side = 2, line = 3, col = 4, cex = 1.5)


subset <- which(loc3d[, 2] == 0)

par(mfrow = c(2, 1))
par(pty = 's')
par(mai = c(0.3, 0.1, 0.6, 0.5))
quilt.plot(loc3d[subset, 1], loc3d[subset, 3], Z1[subset], nx = length(x), ny = length(depth), xlab = "Longitude", ylab = "", zlim = range(Z), ylim = c(1, 0), xaxt = 'n')
mtext("Depth", side = 2, line = 2)
mtext(expression(Z[1]), side = 2, line = 3, col = 4, cex = 1.5)
par(mai = c(0.8, 0.1, 0.1, 0.5))
quilt.plot(loc3d[subset, 1], loc3d[subset, 3], Z2[subset], nx = length(x), ny = length(depth), xlab = "Longitude", ylab = "", zlim = range(Z), ylim = c(1, 0))
mtext("Depth", side = 2, line = 2)
mtext(expression(Z[2]), side = 2, line = 3, col = 4, cex = 1.5)


subset <- which(loc3d[, 1] == 0)

par(mfrow = c(2, 1))
par(pty = 's')
par(mai = c(0.3, 0.1, 0.6, 0.5))
quilt.plot(loc3d[subset, 2], loc3d[subset, 3], Z1[subset], nx = length(y), ny = length(depth), xlab = "Latitude", ylab = "", zlim = range(Z), ylim = c(1, 0), xaxt = 'n')
mtext("Depth", side = 2, line = 2)
mtext(expression(Z[1]), side = 2, line = 3, col = 4, cex = 1.5)
par(mai = c(0.8, 0.1, 0.1, 0.5))
quilt.plot(loc3d[subset, 2], loc3d[subset, 3], Z2[subset], nx = length(y), ny = length(depth), xlab = "Latitude", ylab = "", zlim = range(Z), ylim = c(1, 0))
mtext("Depth", side = 2, line = 2)
mtext(expression(Z[2]), side = 2, line = 3, col = 4, cex = 1.5)

####### COMPUTING EMPIRICAL COVARIANCE #######

emp_cov <- compute_emp_cov(location = loc3d,
                           variable1_residuals = Z1, variable2_residuals = Z2,
                           bandwidth_horizontal = 0.1, bandwidth_vertical = 20,
                           radius = earthRadiusKm)

variance1 <- diag(cov_mat[1:nrow(loc3d), 1:nrow(loc3d)])
variance2 <- diag(cov_mat[nrow(loc3d) + 1:nrow(loc3d), nrow(loc3d) + 1:nrow(loc3d)])
covariance12 <- diag(cov_mat[1:nrow(loc3d), nrow(loc3d) + 1:nrow(loc3d)])
correlation12 <- covariance12 / sqrt(variance1 * variance2)

emp_variance1 <- diag(emp_cov[1:nrow(loc3d), 1:nrow(loc3d)])
emp_variance2 <- diag(emp_cov[nrow(loc3d) + 1:nrow(loc3d), nrow(loc3d) + 1:nrow(loc3d)])
emp_covariance12 <- diag(emp_cov[1:nrow(loc3d), nrow(loc3d) + 1:nrow(loc3d)])
emp_correlation12 <- emp_covariance12 / sqrt(emp_variance1 * emp_variance2)

par(mfrow = c(2, 3))

plot(variance1[1:10])
plot(variance2[1:10])
plot(correlation12[1:10])

plot(emp_variance1[1:10])
plot(emp_variance2[1:10])
plot(emp_correlation12[1:10])

plot(emp_cov[1, seq(1, 1000, by = 10)])


####### ESTIMATION #######

KNOTS1 <- c(0.1, 0.3, 0.5, 0.7, 0.9) #c(0.1, 0.5, 0.9)
KNOTS2 <- c(0.1, 0.3, 0.5, 0.7, 0.9) #c(0.1, 0.5, 0.9)

INIT_BETA = 1
INIT_SCALE_HORIZONTAL = log(0.02)
INIT_SCALE_VERTICAL = log(0.4)
INIT_A1 = INIT_A2 = 0.005
INIT_B1 = INIT_B2 = 0.005
INIT_D1 = INIT_D2 = 0

SPLINES_DEGREE = 2

set.seed(1235)
INIT_C1 <- runif(length(KNOTS1) + SPLINES_DEGREE + 1, -0.1, 0.1)

set.seed(1236)
INIT_C2 <- runif(length(KNOTS2) + SPLINES_DEGREE + 1, -0.1, 0.1)

est_params_wls <- est_bi_differential_wls(empirical_values = emp_cov, location = loc3d,
                                      init_beta = INIT_BETA,
                                      init_scale_horizontal = SCALE_HORIZONTAL,
                                      init_scale_vertical = SCALE_VERTICAL,
                                      init_a1 = A1, init_b1 = B1,
                                      init_c1_coef = INIT_C1, init_d1 = INIT_D1,
                                      init_a2 = A2, init_b2 = B2,
                                      init_c2_coef = INIT_C2, init_d2 = INIT_D2,
                                      a1_fix = TRUE, b1_fix = TRUE, a2_fix = TRUE, b2_fix = TRUE,
                                      beta_fix = TRUE, scale_horizontal_fix = TRUE, scale_vertical_fix = TRUE, d1_fix = TRUE, d2_fix = TRUE,
                                      radius = earthRadiusKm,
                                      splines_degree = SPLINES_DEGREE,
                                      knots1 = KNOTS1, knots2 = KNOTS2,
                                      w1 = 100, w2 = 100, w12 = 1000,
                                      iterlim = 1000, stepmax = 10, hessian = T)

#theta = c(1.09860651,8.85e-06,9.3e-06,0.55610915,1.06328949,0.71657971,-1.67041401,0.50020981,0.55009588,9.06e-06,9.1e-06,0.76672491,0.85972797,-1.51263424,0.19974999,1.29147789,0.76841933)

theta = c(-0.35699489,-0.306007,0.08675093,-0.03231822,0.01926797,0.00011174,-0.35695048,-0.30596901,0.08674003,-0.03231428,0.01926572,0.0003118)

est_params_wls <- est_bi_differential_wls(empirical_values = emp_cov, location = loc3d,
                                          init_beta = INIT_BETA,
                                          init_scale_horizontal = SCALE_HORIZONTAL,
                                          init_scale_vertical = SCALE_VERTICAL,
                                          init_a1 = A1, init_b1 = B1,
                                          init_c1_coef = theta[1:length(INIT_C1)], init_d1 = INIT_D1,
                                          init_a2 = A2, init_b2 = B2,
                                          init_c2_coef = theta[length(INIT_C1) + 1:length(INIT_C2)], init_d2 = INIT_D2,
                                          a1_fix = TRUE, b1_fix = TRUE, a2_fix = TRUE, b2_fix = TRUE,
                                          beta_fix = TRUE, scale_horizontal_fix = TRUE, scale_vertical_fix = TRUE, d1_fix = TRUE, d2_fix = TRUE,
                                          radius = earthRadiusKm,
                                          splines_degree = SPLINES_DEGREE,
                                          knots1 = KNOTS1, knots2 = KNOTS2,
                                          w1 = 100, w2 = 100, w12 = 0,
                                          iterlim = 1000, stepmax = 10, hessian = T)

wls_est_cov_mat <- cov_bi_differential(location = loc3d, beta = est_params_wls$est_beta, scale_horizontal = est_params_wls$est_scale_horizontal, scale_vertical = est_params_wls$est_scale_vertical, a1 = est_params_wls$est_a1, b1 = est_params_wls$est_b1, c1_coef = est_params_wls$est_c1_coef, d1 = est_params_wls$est_d1, a2 = est_params_wls$est_a2, b2 = est_params_wls$est_b2, c2_coef = est_params_wls$est_c2_coef, d2 = est_params_wls$est_d2, radius = earthRadiusKm, splines_degree = est_params_wls$splines_degree, knots1 = est_params_wls$knots1, knots2 = est_params_wls$knots2)

wls_est_variance1 <- diag(wls_est_cov_mat[1:nrow(loc3d), 1:nrow(loc3d)])
wls_est_variance2 <- diag(wls_est_cov_mat[nrow(loc3d) + 1:nrow(loc3d), nrow(loc3d) + 1:nrow(loc3d)])
wls_est_covariance12 <- diag(wls_est_cov_mat[1:nrow(loc3d), nrow(loc3d) + 1:nrow(loc3d)])
wls_est_correlation12 <- wls_est_covariance12 / sqrt(wls_est_variance1 * wls_est_variance2)

par(mfrow = c(2, 3))

plot(variance1[1:10])
plot(variance2[1:10])
plot(correlation12[1:10])

plot(wls_est_variance1[1:10])
plot(wls_est_variance2[1:10])
plot(wls_est_correlation12[1:10])

plot(emp_cov[1, 1000 + seq(1, 1000, by = 10)])
plot(wls_est_cov_mat[1, 1000 + seq(1, 1000, by = 10)])

plot(emp_cov[1, seq(1, 1000, by = 10)])
plot(wls_est_cov_mat[1, seq(1, 1000, by = 10)])

plot(emp_cov[1000 + 1, 1000 + seq(1, 1000, by = 10)])
plot(wls_est_cov_mat[1000 + 1, 1000 + seq(1, 1000, by = 10)])


est_params_mle <- est_bi_differential_mle(residuals = Z, location = loc3d,
                                      init_beta = theta[1],
                                      init_scale_horizontal = theta[2],
                                      init_scale_vertical = theta[3],
                                      init_a1 = theta[4], init_b1 = theta[5],
                                      init_c1_coef = theta[5 + 1:length(INIT_C1)], init_d1 = INIT_D1,
                                      init_a2 = theta[5 + length(INIT_C1) + 1], init_b2 = theta[5 + length(INIT_C1) + 2],
                                      init_c2_coef = theta[5 + length(INIT_C1) + 2 + 1:length(INIT_C2)], init_d2 = INIT_D2,
                                      d1_fix = TRUE, d2_fix = TRUE,
                                      radius = earthRadiusKm,
                                      splines_degree = SPLINES_DEGREE,
                                      knots1 = KNOTS1, knots2 = KNOTS2,
                                      iterlim = 1, hessian = T)

theta = c(-0.12107742,-4.89338811,-1.1452082,-17.4438933,-0.22523577,-0.89459747,-1.36523283,-2.91636576,-3.28533167,-2.07808692,-2.59283651,0.37166639,0.56941461,2.15740875,1.80154959,-0.11354508,0.67696366,1.84610893,1.93622192)

mle_est_cov_mat <- cov_bi_differential(location = loc3d, beta = est_params_mle$est_beta, scale_horizontal = est_params_mle$est_scale_horizontal, scale_vertical = est_params_mle$est_scale_vertical, a1 = est_params_mle$est_a1, b1 = est_params_mle$est_b1, c1_coef = est_params_mle$est_c1_coef, d1 = est_params_mle$est_d1, a2 = est_params_mle$est_a2, b2 = est_params_mle$est_b2, c2_coef = est_params_mle$est_c2_coef, d2 = est_params_mle$est_d2, radius = earthRadiusKm, splines_degree = est_params_mle$splines_degree, knots1 = est_params_mle$knots1, knots2 = est_params_mle$knots2)

mle_est_variance1 <- diag(mle_est_cov_mat[1:nrow(loc3d), 1:nrow(loc3d)])
mle_est_variance2 <- diag(mle_est_cov_mat[nrow(loc3d) + 1:nrow(loc3d), nrow(loc3d) + 1:nrow(loc3d)])
mle_est_covariance12 <- diag(mle_est_cov_mat[1:nrow(loc3d), nrow(loc3d) + 1:nrow(loc3d)])
mle_est_correlation12 <- mle_est_covariance12 / sqrt(mle_est_variance1 * mle_est_variance2)

par(mfrow = c(2, 3))

plot(variance1[1:10])
plot(variance2[1:10])
plot(correlation12[1:10])

plot(mle_est_variance1[1:10])
plot(mle_est_variance2[1:10])
plot(mle_est_correlation12[1:10])

plot(emp_cov[1, seq(1, 100, by = 10)])
plot(mle_est_cov_mat[1, seq(1, 100, by = 10)])

#return the theta vector also
print("DONE . . . ")

