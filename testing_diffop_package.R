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

KNOTS1 <- c(0.1, 0.5, 0.9)
KNOTS2 <- c(0.1, 0.5, 0.9)

INIT_BETA = 0.6
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
                                      init_scale_horizontal = INIT_SCALE_HORIZONTAL,
                                      init_scale_vertical = INIT_SCALE_VERTICAL,
                                      init_a1 = INIT_A1, init_b1 = INIT_B1,
                                      init_c1_coef = INIT_C1, init_d1 = INIT_D1,
                                      init_a2 = INIT_A2, init_b2 = INIT_B2,
                                      init_c2_coef = INIT_C2, init_d2 = INIT_D2,
                                      d1_fix = TRUE, d2_fix = TRUE,
                                      radius = earthRadiusKm,
                                      splines_degree = SPLINES_DEGREE,
                                      knots1 = KNOTS1, knots2 = KNOTS2,
                                      w1 = 100, w2 = 50000, w12 = 1000,
                                      iterlim = 1000, hessian = T)

theta = c(1.58e-06,-4.7813677,1.04131888,-3.1830562,0.04081191,0.00159106,-0.00666072,-0.12630797,0.1927597,-0.13017259,-0.06530689,-0.01515386,-0.0457482,0.01752583,0.02381038,-0.04111447,0.00381765,0.03179052,0.01964975)

est_params_wls <- est_bi_differential_wls(empirical_values = emp_cov, location = loc3d,
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
                                          w1 = 100, w2 = 50000, w12 = 1000,
                                          iterlim = 5, hessian = T)

theta = c(-0.05355685,-4.10759032,1.66002015,-5.67731641,-0.11125447,-0.02032645,-0.06730041,0.0420658,-0.07155357,0.03149789,0.0051676,-0.03456796,-0.00634583,0.00859718,0.01393778,-0.02207468,0.00269277,0.01529831,0.01034171)

est_cov_mat <- cov_bi_differential(location = loc3d, beta = est_params_wls$est_beta, scale_horizontal = est_params_wls$est_scale_horizontal, scale_vertical = est_params_wls$est_scale_vertical, a1 = est_params_wls$est_a1, b1 = est_params_wls$est_b1, c1_coef = est_params_wls$est_c1_coef, d1 = est_params_wls$est_d1, a2 = est_params_wls$est_a2, b2 = est_params_wls$est_b2, c2_coef = est_params_wls$est_c2_coef, d2 = est_params_wls$est_d2, radius = earthRadiusKm, splines_degree = est_params_wls$splines_degree, knots1 = est_params_wls$knots1, knots2 = est_params_wls$knots2)

est_variance1 <- diag(est_cov_mat[1:nrow(loc3d), 1:nrow(loc3d)])
est_variance2 <- diag(est_cov_mat[nrow(loc3d) + 1:nrow(loc3d), nrow(loc3d) + 1:nrow(loc3d)])
est_covariance12 <- diag(est_cov_mat[1:nrow(loc3d), nrow(loc3d) + 1:nrow(loc3d)])
est_correlation12 <- est_covariance12 / sqrt(est_variance1 * est_variance2)

par(mfrow = c(2, 3))

plot(variance1[1:10])
plot(variance2[1:10])
plot(correlation12[1:10])

plot(est_variance1[1:10])
plot(est_variance2[1:10])
plot(est_correlation12[1:10])

plot(emp_cov[1, seq(1, 100, by = 10)])
plot(est_cov_mat[1, seq(1, 100, by = 10)])


est_params <- est_bi_differential_mle(residuals = Z, location = loc3d, init_beta = est_params_wls$est_beta,
                                  init_scale_horizontal = est_params_wls$est_scale_horizontal,
                                  init_scale_vertical = est_params_wls$est_scale_vertical,
                                  init_a1 = est_params_wls$est_a1, init_b1 = est_params_wls$est_b1,
                                  init_c1_coef = est_params_wls$est_c1_coef, init_d1 = est_params_wls$est_d1,
                                  init_a2 = est_params_wls$est_a2, init_b2 = est_params_wls$est_b2,
                                  init_c2_coef = est_params_wls$est_c2_coef, init_d2 = est_params_wls$est_d2,
                                  d1_fix = TRUE, d2_fix = TRUE,
                                  radius = earthRadiusKm,
                                  splines_degree = est_params_wls$splines_degree,
                                  knots1 = est_params_wls$knots, knots2 = est_params_wls$knots2,
                                  iterlim = 1000, hessian = T)


print("DONE . . . ")

