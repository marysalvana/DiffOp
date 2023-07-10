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

INIT_BETA = 0
INIT_SCALE_HORIZONTAL = log(0.02)
INIT_SCALE_VERTICAL = log(0.2)
INIT_A1 = INIT_A2 = 0
INIT_B1 = INIT_B2 = 0
INIT_D1 = INIT_D2 = 0
INIT_C1_COEF = INIT_C2_COEF = 1

est_params_mle_simple <- est_bi_differential_mle(residuals = Z, location = loc3d,
                                          init_beta = INIT_BETA,
                                          init_scale_horizontal = INIT_SCALE_HORIZONTAL,
                                          init_scale_vertical = INIT_SCALE_VERTICAL,
                                          init_a1 = INIT_A1, init_b1 = INIT_B1,
                                          init_c1_coef = INIT_C1_COEF, init_d1 = INIT_D1,
                                          init_a2 = INIT_A2, init_b2 = INIT_B2,
                                          init_c2_coef = INIT_C2_COEF, init_d2 = INIT_D2,
                                          d1_fix = TRUE, d2_fix = TRUE,
                                          radius = earthRadiusKm,
                                          splines_degree = 0,
                                          iterlim = 1000, stepmax = 1, hessian = T)

#theta = c(-0.07594104,-5.09094054,-1.34307397,1.65761039,4.10280701,4.2595846,-4.60992802,-2.15566869,3.3232618)

#mle_est_cov_mat <- cov_bi_differential(location = loc3d, beta = -0.03795228, scale_horizontal = exp(-5.09094054), scale_vertical = exp(-1.34307397), a1 = 1.65761039 * 1e-3, b1 = 4.10280701 * 1e-3, c1_coef = 4.2595846, d1 = 0,
#                                       a2 = -4.60992802 * 1e-3, b2 = -2.15566869 * 1e-3, c2_coef = 3.3232618, d2 = 0, radius = earthRadiusKm, splines_degree = 0)

KNOTS1 <- c(0.1, 0.5, 0.9)
KNOTS2 <- c(0.1, 0.5, 0.9)

SPLINES_DEGREE = 2

set.seed(1235)
INIT_C1_COEF <- runif(length(KNOTS1) + SPLINES_DEGREE + 1, -0.1, 0.1)

set.seed(1236)
INIT_C2_COEF <- runif(length(KNOTS2) + SPLINES_DEGREE + 1, -0.1, 0.1)

#theta = c(-1.48154972,-3.03783575,-1.50604237,2.89812541,5.33625024,2.84475697,0.31430931,0.35702328,-1.32826959,-0.50023631,0.08196886,0.10911162)
#INIT_C1_COEF <- theta[1:length(INIT_C1_COEF)]
#INIT_C2_COEF <- theta[length(INIT_C1_COEF) + 1:length(INIT_C2_COEF)]

theta = c(0.08732791,2.36984792,5.77820957,-1.19797418,-3.21585704,-1.09787706,2.78784418,5.00726896,1.92599343,-5.44884928,-2.88333783,0.90521174,0.34422214,-1.32281772,-0.70832396,0.07114146,0.22650501)
theta = c(-0.62804552,2.96923077,6.55905817,-2.96441859,-3.47577949,-1.10837926,4.75604482,4.6003741,3.91656508,-6.64612809,-3.81412996,2.75397086,2.07701242,-4.23989661,-4.25764403,-4.40627351,-3.97551753)

est_params_mle <- est_bi_differential_mle(residuals = Z, location = loc3d,
                                          init_beta = theta[1],
                                          init_scale_horizontal = exp(-5.09094054),
                                          init_scale_vertical = exp(-1.34307397),
                                          init_a1 = theta[2], init_b1 = theta[3],
                                          init_c1_coef = theta[3 + 1:length(INIT_C1_COEF)], init_d1 = 0,
                                          init_a2 = theta[3 + length(INIT_C1_COEF) + 1], init_b2 = theta[3 + length(INIT_C1_COEF) + 2],
                                          init_c2_coef = theta[3 + length(INIT_C1_COEF) + 2 + 1:length(INIT_C2_COEF)], init_d2 = 0,
                                          scale_horizontal_fix = TRUE, scale_vertical_fix = TRUE,
                                          d1_fix = TRUE, d2_fix = TRUE,
                                          radius = earthRadiusKm,
                                          splines_degree = SPLINES_DEGREE,
                                          knots1 = KNOTS1, knots2 = KNOTS2,
                                          iterlim = 1000, stepmax = 1, hessian = T)

#theta = c(-4.65560477,-10.57604634,13.70004287,25.78831257,3.49253491,4.66280469,3.15236873,-0.63364572,-4.00411605,-5.60298597,-4.6558888,-5.21335133)
#theta = c(-4.05506591,-3.1916669,-1.59170792,5.46825639,4.9853047,3.94290889,4.39038617,2.5052013,-5.33414962,-5.07497427,-4.92331565,-4.68664234)
#negloglik:  -3669.4403

mle_est_cov_mat <- cov_bi_differential(location = loc3d, beta = -0.03795228, scale_horizontal = exp(-5.09094054), scale_vertical = exp(-1.34307397), a1 = 1.65761039 * 1e-3, b1 = 4.10280701 * 1e-3, c1_coef = theta[1:length(INIT_C1_COEF)], d1 = 0,
                                       a2 = -4.60992802 * 1e-3, b2 = -2.15566869 * 1e-3, c2_coef = theta[length(INIT_C1_COEF) + 1:length(INIT_C2_COEF)], d2 = 0, radius = earthRadiusKm, splines_degree = SPLINES_DEGREE, knots1 = KNOTS1, knots2 = KNOTS2)



est_params_mle <- est_bi_differential_mle(residuals = Z, location = loc3d,
                                          init_beta = theta[1],
                                          init_scale_horizontal = exp(-5.09094054),
                                          init_scale_vertical = exp(-1.34307397),
                                          init_a1 = theta[2], init_b1 = theta[3],
                                          init_c1_coef = theta[3 + 1:length(INIT_C1_COEF)], init_d1 = 0,
                                          init_a2 = theta[3 + length(INIT_C1_COEF) + 1], init_b2 = theta[3 + length(INIT_C1_COEF) + 2],
                                          init_c2_coef = theta[3 + length(INIT_C1_COEF) + 2 + 1:length(INIT_C2_COEF)], init_d2 = 0,
                                          c1_fix = TRUE, c2_fix = TRUE,
                                          d1_fix = TRUE, d2_fix = TRUE,
                                          radius = earthRadiusKm,
                                          splines_degree = SPLINES_DEGREE,
                                          knots1 = KNOTS1, knots2 = KNOTS2,
                                          iterlim = 1000, stepmax = 1, hessian = T)



est_params_wls <- est_bi_differential_wls(empirical_values = emp_cov, location = loc3d,
                                      init_beta = BETA,
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
                                      w1 = 1, w2 = 1, w12 = 1,
                                      iterlim = 1000, stepmax = 10, hessian = T)
#w12 = 0: minimum = 3.762514, with cov_mat
theta = c(-0.58824841,-1.12191566,-0.76319908,1.77320153,-0.53293114,-0.58714179,0.80815322,0.92454372,-1.61688095,0.22649478,1.35237963,0.80979308)

#w12 = 0: minimum = 120.9826, with emp_cov
theta = c(0.61850536,1.06475403,0.49700978,1.41275038,0.19323393,0.31166667,-0.13662041,-0.01967288,-0.51526506,-0.15759683,-0.24550915,-0.18268438)
#w12 = 0.1: minimum = 121.4423, with emp_cov
theta = c(0.61797399,1.06402943,0.49691149,1.41179833,0.19326209,0.31153443,-0.11824864,-0.01172988,-0.48771374,-0.14545836,-0.23689136,-0.17442548)
#w12 = 1: minimum = 123.5051, with emp_cov
theta = c(0.61910027,1.06540235,0.49718786,1.41199994,0.19333108,0.31161038,0.06939223,0.05087207,-0.06649795,0.01196319,-0.06364634,-0.03215431)

wls_est_cov_mat <- cov_bi_differential(location = loc3d, beta = est_params_wls$est_beta, scale_horizontal = est_params_wls$est_scale_horizontal, scale_vertical = est_params_wls$est_scale_vertical, a1 = est_params_wls$est_a1, b1 = est_params_wls$est_b1, c1_coef = est_params_wls$est_c1_coef, d1 = est_params_wls$est_d1, a2 = est_params_wls$est_a2, b2 = est_params_wls$est_b2, c2_coef = est_params_wls$est_c2_coef, d2 = est_params_wls$est_d2, radius = earthRadiusKm, splines_degree = est_params_wls$splines_degree, knots1 = est_params_wls$knots1, knots2 = est_params_wls$knots2)

wls_est_variance1 <- diag(wls_est_cov_mat[1:nrow(loc3d), 1:nrow(loc3d)])
wls_est_variance2 <- diag(wls_est_cov_mat[nrow(loc3d) + 1:nrow(loc3d), nrow(loc3d) + 1:nrow(loc3d)])
wls_est_covariance12 <- diag(wls_est_cov_mat[1:nrow(loc3d), nrow(loc3d) + 1:nrow(loc3d)])
wls_est_correlation12 <- wls_est_covariance12 / sqrt(wls_est_variance1 * wls_est_variance2)

par(mfrow = c(2, 3))

plot(emp_variance1[1:10])
plot(emp_variance2[1:10])
plot(emp_correlation12[1:10])

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
                                      init_beta = BETA,
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
                                      iterlim = 1000, stepmax = 10, hessian = T)

theta = c(1.02072349,1.70698695,4.4786531,7.22314096,-1.40800868,13.91150388,10.40307889,-5.24999602,-3.55195695,-2.24217283,-1.63818859,-1.11613754)

est_params_mle <- est_bi_differential_mle(residuals = Z, location = loc3d,
                                          init_beta = 1.1,
                                          init_scale_horizontal = log(SCALE_HORIZONTAL),
                                          init_scale_vertical = log(SCALE_VERTICAL),
                                          init_a1 = A1, init_b1 = B1,
                                          init_c1_coef = theta[1:length(INIT_C1)], init_d1 = INIT_D1,
                                          init_a2 = A2, init_b2 = B2,
                                          init_c2_coef = theta[length(INIT_C1) + 1:length(INIT_C2)], init_d2 = INIT_D2,
                                          d1_fix = TRUE, d2_fix = TRUE,
                                          radius = earthRadiusKm,
                                          splines_degree = SPLINES_DEGREE,
                                          knots1 = KNOTS1, knots2 = KNOTS2,
                                          iterlim = 1000, stepmax = 1, hessian = T)

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

plot(emp_cov[2, seq(1, 100, by = 10)])
plot(mle_est_cov_mat[5, seq(1, 100, by = 10)])

#return the theta vector also
print("DONE . . . ")

