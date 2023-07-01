NEGLOGLIK1 <- function(theta, empirical_values, basis1, basis2, radius){

  print(paste("theta = c(", paste(round(theta, 8), collapse=","), ")", sep = ''))

  nb1 <- ncol(basis1)
  nb2 <- ncol(basis2)

  BETA <- 0 #start with 0 to get the marginals correctly first, then slowly keep increasing to 1. For B4, go up to the minimum rho that will minimize the WLS

  SCALE_HORIZONTAL <- 0.024995 #use the I1 estimate from overleaf
  SCALE_VERTICAL <- 0.005251 # start with 0.5, then keep decreasing to 0.024813

  A1 <- 0.01 * (1 / (1 + exp(--0.071505)))
  B1 <- 0.02 * (1 / (1 + exp(-0.60286))) - 0.02/2
  C1_coef <- theta#[2 + 1:nb1]

  A2 <- 0
  B2 <- 0
  C2_coef <- rep(0.1, nb1)

  D1 <- 0
  D2 <- 0

  if(theta[1] < -10){
    return(Inf)
  }

  C1_param <- basis1 %*% matrix(C1_coef, ncol = 1)
  C2_param <- basis2 %*% matrix(C2_coef, ncol = 1)

  cov_mat <- cov_bi_differential(location = location, beta = BETA,
                                 scale_horizontal = SCALE_HORIZONTAL, scale_vertical = SCALE_VERTICAL,
                                 a1 = A1, b1 = B1, c1 = C1_param, d1 = D1, a2 = A2, b2 = B2, c2 = C2_param, d2 = D2,
                                 radius = radius)
  #variance1 <- diag(cov_mat[1:nrow(location), 1:nrow(location)])
  #variance2 <- diag(cov_mat[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)])

  #emp_variance1 <- diag(empirical_values[1:nrow(location), 1:nrow(location)])
  #emp_variance2 <- diag(empirical_values[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)])

  covariance1 <- cov_mat[1:nrow(location), 1:nrow(location)]
  covariance2 <- cov_mat[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)]
  covariance12 <- cov_mat[1:nrow(location), nrow(location) + 1:nrow(location)]
  #correlation12 <- covariance12 / outer(sqrt(variance1), sqrt(variance2), '*')
  correlation12 <- diag(covariance12) / sqrt(diag(covariance1) * diag(covariance2))

  emp_covariance1 <- empirical_values[1:nrow(location), 1:nrow(location)]
  emp_covariance2 <- empirical_values[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)]
  emp_covariance12 <- empirical_values[1:nrow(location), nrow(location) + 1:nrow(location)]
  #emp_correlation12 <- emp_covariance12 / outer(sqrt(emp_variance1), sqrt(emp_variance2), '*')
  emp_correlation12 <- diag(emp_covariance12) / sqrt(diag(emp_covariance1) * diag(emp_covariance2))

  plot_bi_differential(location = location[1:50, ], est_beta = BETA,
                       est_scale_horizontal = SCALE_HORIZONTAL, est_scale_vertical = SCALE_VERTICAL,
                       est_a1 = A1, est_b1 = B1, est_c1 = C1_coef, est_d1 = D1, est_a2 = A2, est_b2 = B2, est_c2 = C2_coef, est_d2 = D2,
                       basis1 = basis1[1:50, ], basis2 = basis2[1:50, ], radius = radius)

  out <- sum((diag(emp_covariance1) - diag(covariance1))^2 * 100 + (diag(emp_covariance2) - diag(covariance2))^2 * 50000 + (emp_correlation12 - correlation12)^2 * 100)
  #out <- sum((diag(emp_covariance1) - diag(covariance1))^2 * 10 + (diag(emp_covariance2) - diag(covariance2))^2 * 1000 + (diag(emp_correlation12) - diag(correlation12))^2 * 0)

  return(out)
  #}
}

NEGLOGLIK2 <- function(theta, empirical_values, basis1, basis2, radius){

  print(paste("theta = c(", paste(round(theta, 8), collapse=","), ")", sep = ''))

  nb1 <- ncol(basis1)
  nb2 <- ncol(basis2)

  BETA <- 0 #start with 0 to get the marginals correctly first, then slowly keep increasing to 1. For B4, go up to the minimum rho that will minimize the WLS

  SCALE_HORIZONTAL <- 0.024995 #use the I1 estimate from overleaf
  SCALE_VERTICAL <- 0.005251 # start with 0.5, then keep decreasing to 0.024813

  A1 <- 0.01 * (1 / (1 + exp(--0.071505)))
  B1 <- 0.02 * (1 / (1 + exp(-0.60286))) - 0.02/2
  C1_coef <- c(-9.19221095,22.91040418,-34.30315502,138.3237635,110.85729951,125.42288729,-29.69320553,15.78731606,-4.54615698)

  A2 <- 0.02 * (1 / (1 + exp(-0.139186))) - 0.02/2
  B2 <- 0.02 * (1 / (1 + exp(-0.000412))) - 0.02/2
  C2_coef <- theta#[2 + 1:nb2]

  D1 <- 0
  D2 <- 0

  if(theta[1] < -10){
    return(Inf)
  }

  C1_param <- basis1 %*% matrix(C1_coef, ncol = 1)
  C2_param <- basis2 %*% matrix(C2_coef, ncol = 1)

  cov_mat <- cov_bi_differential(location = location, beta = BETA,
                                 scale_horizontal = SCALE_HORIZONTAL, scale_vertical = SCALE_VERTICAL,
                                 a1 = A1, b1 = B1, c1 = C1_param, d1 = D1, a2 = A2, b2 = B2, c2 = C2_param, d2 = D2,
                                 radius = radius)

  covariance1 <- cov_mat[1:nrow(location), 1:nrow(location)]
  covariance2 <- cov_mat[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)]
  covariance12 <- cov_mat[1:nrow(location), nrow(location) + 1:nrow(location)]
  #correlation12 <- covariance12 / outer(sqrt(variance1), sqrt(variance2), '*')
  correlation12 <- diag(covariance12) / sqrt(diag(covariance1) * diag(covariance2))

  emp_covariance1 <- empirical_values[1:nrow(location), 1:nrow(location)]
  emp_covariance2 <- empirical_values[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)]
  emp_covariance12 <- empirical_values[1:nrow(location), nrow(location) + 1:nrow(location)]
  #emp_correlation12 <- emp_covariance12 / outer(sqrt(emp_variance1), sqrt(emp_variance2), '*')
  emp_correlation12 <- diag(emp_covariance12) / sqrt(diag(emp_covariance1) * diag(emp_covariance2))

  plot_bi_differential(location = location[1:50, ], est_beta = BETA,
                       est_scale_horizontal = SCALE_HORIZONTAL, est_scale_vertical = SCALE_VERTICAL,
                       est_a1 = A1, est_b1 = B1, est_c1 = C1_coef, est_d1 = D1, est_a2 = A2, est_b2 = B2, est_c2 = C2_coef, est_d2 = D2,
                       basis1 = basis1[1:50, ], basis2 = basis2[1:50, ], radius = radius)

  out <- sum((diag(emp_covariance1) - diag(covariance1))^2 * 100 + (diag(emp_covariance2) - diag(covariance2))^2 * 50000 + (emp_correlation12 - correlation12)^2 * 100)

  return(out)
}

NEGLOGLIK3 <- function(theta, empirical_values, basis1, basis2, radius){

  print(paste("theta = c(", paste(round(theta, 8), collapse=","), ")", sep = ''))

  nb1 <- ncol(basis1)
  nb2 <- ncol(basis2)

  BETA <- 1 #start with 0 to get the marginals correctly first, then slowly keep increasing to 1. For B4, go up to the minimum rho that will minimize the WLS

  if(REF_LOC_IND == 2){
    #SCALE_HORIZONTAL <- 0.251897#use the I1 estimate from overleaf
    #SCALE_VERTICAL <- 0.017972
  }else if(REF_LOC_IND == 3){
    #SCALE_HORIZONTAL <- 0.220331
    #SCALE_VERTICAL <- 0.009407
  }else if(REF_LOC_IND == 5){
    #SCALE_HORIZONTAL <- 0.048388
    #SCALE_VERTICAL <- 0.014254
  }else if(REF_LOC_IND == 6){
    #SCALE_HORIZONTAL <- 0.114639
    #SCALE_VERTICAL <- 0.008547
  }

  SCALE_HORIZONTAL <- 0.181846
  SCALE_VERTICAL <- 0.016428 # start with 0.5, then keep decreasing to the true SCALE_VERTICAL value

  A1 <- 0.01 * (1 / (1 + exp(-theta[1])))
  B1 <- 0.02 * (1 / (1 + exp(-theta[2]))) - 0.02/2
  #C1_coef <- c(38.849119,29.097003,63.582509,181.502465,12.4488,16.071857,0.615395,-0.440928,-2.497189)
  #C1_coef <- c(8.258865,-7.10906,51.38904,152.453364,-34.652917,6.645691,-1.692241,1.427228,-0.577232)
  #C1_coef <- c(5.214247,-9.29537,52.994906,152.864557,-35.414152,6.746395,-1.511015,1.176445,-0.289005)

  A2 <- -3.6e-05#0.02 * (1 / (1 + exp(-theta[3]))) - 0.02/2
  B2 <- 2.7e-05#0.02 * (1 / (1 + exp(-theta[4]))) - 0.02/2
  #C2_coef <- c(21.653443,14.253736,24.879136,11.220403,-1.483879,1.256915,-0.594413,0.763171,-0.508831)
  #C2_coef <- c(14.289976,12.129137,12.005572,10.930104,-2.963638,0.943934,-0.277868,0.247353,-0.094527)
  #C2_coef <- c(14.97932,13.474747,14.219633,12.913269,-2.851019,0.479333,-0.003195,0.25214,-0.081682)

  D1 <- 0
  D2 <- 0

  if(abs(theta[1]) > 10){
    #return(Inf)
  }

  C1_param <- basis1 %*% matrix(C1_coef, ncol = 1)
  C2_param <- basis2 %*% matrix(C2_coef, ncol = 1)

  cov_mat <- cov_bi_differential(location = location, beta = BETA,
                                 scale_horizontal = SCALE_HORIZONTAL, scale_vertical = SCALE_VERTICAL,
                                 a1 = A1, b1 = B1, c1 = C1_param, d1 = D1, a2 = A2, b2 = B2, c2 = C2_param, d2 = D2,
                                 radius = radius)
  #variance1 <- diag(cov_mat[1:nrow(location), 1:nrow(location)])
  #variance2 <- diag(cov_mat[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)])

  #emp_variance1 <- diag(empirical_values[1:nrow(location), 1:nrow(location)])
  #emp_variance2 <- diag(empirical_values[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)])

  covariance1 <- cov_mat[1:nrow(location), 1:nrow(location)]
  covariance2 <- cov_mat[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)]
  covariance12 <- cov_mat[1:nrow(location), nrow(location) + 1:nrow(location)]
  #correlation12 <- covariance12 / outer(sqrt(variance1), sqrt(variance2), '*')
  correlation12 <- diag(covariance12) / sqrt(diag(covariance1) * diag(covariance2))

  emp_covariance1 <- empirical_values[1:nrow(location), 1:nrow(location)]
  emp_covariance2 <- empirical_values[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)]
  emp_covariance12 <- empirical_values[1:nrow(location), nrow(location) + 1:nrow(location)]
  #emp_correlation12 <- emp_covariance12 / outer(sqrt(emp_variance1), sqrt(emp_variance2), '*')
  emp_correlation12 <- diag(emp_covariance12) / sqrt(diag(emp_covariance1) * diag(emp_covariance2))

  plot_bi_differential(location = location[1:50, ], est_beta = BETA,
                       est_scale_horizontal = SCALE_HORIZONTAL, est_scale_vertical = SCALE_VERTICAL,
                       est_a1 = A1, est_b1 = B1, est_c1 = C1_coef, est_d1 = D1, est_a2 = A2, est_b2 = B2, est_c2 = C2_coef, est_d2 = D2,
                       basis1 = basis1[1:50, ], basis2 = basis2[1:50, ], radius = radius)

  out <- sum((diag(emp_covariance1) - diag(covariance1))^2 * 100 + (diag(emp_covariance2) - diag(covariance2))^2 * 50000 + (emp_correlation12 - correlation12)^2 * 100)

  return(out)
  #}
}

NEGLOGLIK4 <- function(theta, emp_covariance1, emp_covariance2, emp_covariance, basis1, basis2, radius){

  print(paste("theta = c(", paste(round(theta, 8), collapse=","), ")", sep = ''))

  nb1 <- ncol(basis1)
  nb2 <- ncol(basis2)

  BETA <- 1 #start with 0 to get the marginals correctly first, then slowly keep increasing to 1. For B4, go up to the minimum rho that will minimize the WLS

  if(REF_LOC_IND == 1){
    #SCALE_HORIZONTAL <- 0.051631#use the I1 estimate from overleaf
    #SCALE_VERTICAL <- 0.007955
  }else if(REF_LOC_IND == 2){
    #SCALE_HORIZONTAL <- 0.181705#use the I1 estimate from overleaf
    #SCALE_VERTICAL <- 0.016264
  }else if(REF_LOC_IND == 3){
    #SCALE_HORIZONTAL <- 0.220331
    #SCALE_VERTICAL <- 0.009407
  }else if(REF_LOC_IND == 4){
    #SCALE_HORIZONTAL <- 0.02506
    #SCALE_VERTICAL <- 0.005254
  }else if(REF_LOC_IND == 5){
    #SCALE_HORIZONTAL <- 0.021257
    #SCALE_VERTICAL <- 0.047564
  }else if(REF_LOC_IND == 6){
    #SCALE_HORIZONTAL <- 0.076977
    #SCALE_VERTICAL <- 0.013413
  }

  SCALE_HORIZONTAL <- 0.231705#0.278977#0.07744#0.5 * (1 / (1 + exp(-theta[1]))) #0.11#0.251874
  SCALE_VERTICAL <- 0.018841# start with 0.5, then keep decreasing to the true SCALE_VERTICAL value

  A1 <- 0.01 * (1 / (1 + exp(-theta[1])))
  B1 <- 0.02 * (1 / (1 + exp(-theta[2]))) - 0.02/2
  C1_coef <- theta[2 + 1:nb1]

  A2 <- 0.02 * (1 / (1 + exp(-theta[2 + nb1 + 1]))) - 0.02/2
  B2 <- 0.02 * (1 / (1 + exp(-theta[2 + nb1 + 2]))) - 0.02/2
  C2_coef <- theta[2 + nb1 + 2 + 1:nb2]

  D1 <- 0
  D2 <- 0

  if(abs(theta[1]) > 5 | abs(theta[2]) > 5){
    #return(Inf)
  }

  C1_param <- basis1 %*% matrix(C1_coef, ncol = 1)
  C2_param <- basis2 %*% matrix(C2_coef, ncol = 1)

  cov_mat <- cov_bi_differential(location = location, beta = BETA,
                                 scale_horizontal = SCALE_HORIZONTAL, scale_vertical = SCALE_VERTICAL,
                                 a1 = A1, b1 = B1, c1 = C1_param, d1 = D1, a2 = A2, b2 = B2, c2 = C2_param, d2 = D2,
                                 radius = radius)
  #variance1 <- diag(cov_mat[1:nrow(location), 1:nrow(location)])
  #variance2 <- diag(cov_mat[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)])

  #emp_variance1 <- diag(empirical_values[1:nrow(location), 1:nrow(location)])
  #emp_variance2 <- diag(empirical_values[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)])

  covariance1 <- cov_mat[1:nrow(location), 1:nrow(location)]
  covariance2 <- cov_mat[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)]
  covariance12 <- cov_mat[1:nrow(location), nrow(location) + 1:nrow(location)]
  #correlation12 <- covariance12 / outer(sqrt(variance1), sqrt(variance2), '*')
  correlation12 <- diag(covariance12) / sqrt(diag(covariance1) * diag(covariance2))

  #emp_correlation12 <- emp_covariance12 / outer(sqrt(emp_variance1), sqrt(emp_variance2), '*')
  emp_correlation12 <- diag(emp_covariance) / sqrt(diag(emp_covariance1) * diag(emp_covariance2))


  plot_bi_differential(location = location[1:50, ], est_beta = BETA,
                       est_scale_horizontal = SCALE_HORIZONTAL, est_scale_vertical = SCALE_VERTICAL,
                       est_a1 = A1, est_b1 = B1, est_c1 = C1_coef, est_d1 = D1, est_a2 = A2, est_b2 = B2, est_c2 = C2_coef, est_d2 = D2,
                       basis1 = basis1[1:50, ], basis2 = basis2[1:50, ], radius = radius, emp_covariance1 = emp_covariance1[1:50, 1:50], emp_covariance2 = emp_covariance2[1:50, 1:50], emp_correlation12 = emp_correlation12[1:50], splines_degree = 2)

  #out <- sum((diag(emp_covariance1) - diag(covariance1))^2 * 100 + (diag(emp_covariance2) - diag(covariance2))^2 * 50000 + (emp_correlation12 - correlation12)^2 * 100)
  out <- sum((diag(emp_covariance1) - diag(covariance1))^2 * 1000 + (diag(emp_covariance2) - diag(covariance2))^2 * 500000 + (emp_correlation12 - correlation12)^2 * 100)
  #out <- sum((diag(emp_covariance1) - diag(covariance1))^2 * 1000)
  #out <- sum((diag(emp_covariance2) - diag(covariance2))^2 * 1000)
  #out <- sum((emp_correlation12 - correlation12)^2 * 100)
  #out <- sum((diag(emp_covariance1_sub) - diag(covariance1))^2 * 500000000 + (diag(emp_covariance2_sub) - diag(covariance2))^2 * 100000000000 + (emp_correlation12 - correlation12)^2 * 1000000000)

  #plot(diag(emp_covariance1)[1:50], location[1:50, 3], type = 'l', xlab = 'Variance', ylab = 'Depth', ylim = c(max(location[, 3]), 0))

  #for(ll in 2:10){
    #lines(diag(emp_covariance1)[(ll - 1) * 50 + 1:50], location[(ll - 1) * 50 + 1:50, 3], lwd = 2, col = 2)
  #}

  #plot(diag(covariance1)[1:50], location[1:50, 3], type = 'l', xlab = 'Variance', ylab = 'Depth', ylim = c(max(location[, 3]), 0))

  #for(ll in 1:10){
    #lines(diag(covariance1)[(ll - 1) * 50 + 1:50], location[(ll - 1) * 50 + 1:50, 3], lwd = 2, col = 4)
  #}

  return(out)
  #}
}


################### SYNTHETIC MARGINALS FROM I3 ESTIMATES   ###################

theta0 = c(-19.798219,-4.798902,-62.215354,-122.38621,28.253029,1.006228,-0.235633,-0.267975,-0.097607,14.206361,11.169431,11.410739,12.127465,-2.225656,0.644952,-0.164696,-1e-06,-0.109094)

BETA = 0

SCALE_HORIZONTAL <- 0.251897
SCALE_VERTICAL <- 0.017519

A1 <- 0.000108
B1 <- 6.3e-05
C1_coef <- theta0[1:nb1]

A2 <- -2e-06
B2 <- 2e-06
C2_coef <- theta0[nb1 + 1:nb2]

C1_param <- basis1 %*% matrix(C1_coef, ncol = 1)
C2_param <- basis2 %*% matrix(C2_coef, ncol = 1)

D1 = D2 = 0

cov_mat_marginals <- cov_bi_differential(location = location, beta = BETA,
                                         scale_horizontal = SCALE_HORIZONTAL, scale_vertical = SCALE_VERTICAL,
                                         a1 = A1, b1 = B1, c1 = C1_param, d1 = D1, a2 = A2, b2 = B2, c2 = C2_param, d2 = D2,
                                         radius = radius)

###################   ###################

NEGLOGLIK5 <- function(theta, empirical_values, basis1, basis2, radius){

  print(paste("theta = c(", paste(round(theta, 8), collapse=","), ")", sep = ''))

  nb1 <- ncol(basis1)
  nb2 <- ncol(basis2)

  BETA <- 0.7 #start with 0 to get the marginals correctly first, then slowly keep increasing to 1. For B4, go up to the minimum rho that will minimize the WLS

  if(REF_LOC_IND == 1){
    #SCALE_HORIZONTAL <- 0.051631#use the I1 estimate from overleaf
    #SCALE_VERTICAL <- 0.007955
  }else if(REF_LOC_IND == 2){
    #SCALE_HORIZONTAL <- 0.181705#use the I1 estimate from overleaf
    #SCALE_VERTICAL <- 0.016264
  }else if(REF_LOC_IND == 3){
    #SCALE_HORIZONTAL <- 0.220331
    #SCALE_VERTICAL <- 0.009407
  }else if(REF_LOC_IND == 4){
    #SCALE_HORIZONTAL <- 0.02506
    #SCALE_VERTICAL <- 0.005254
  }else if(REF_LOC_IND == 5){
    #SCALE_HORIZONTAL <- 0.021257
    #SCALE_VERTICAL <- 0.047564
  }else if(REF_LOC_IND == 6){
    #SCALE_HORIZONTAL <- 0.076977
    #SCALE_VERTICAL <- 0.013413
  }

  SCALE_HORIZONTAL <- 0.111561#0.181846
  SCALE_VERTICAL <- 0.017345#0.014252 # start with 0.5, then keep decreasing to the true SCALE_VERTICAL value

  #theta0 = c(7.96862592,1.36275107,47.54151561,129.13662873,14.44771646,38.0088246,1.99247066,9.34189849,-29.04960955,12.07961682,14.29899901,8.12003918,15.27748862,0.63911336,0.89385433,-7.39657169,6.87994068,5.62571323)

  A1 <- 0.00012#0.01 * (1 / (1 + exp(-theta[1])))
  B1 <- 0.000343#0.02 * (1 / (1 + exp(-theta[2]))) - 0.02/2
  C1_coef <- theta[1:nb1]

  A2 <- -5.6e-05#0.02 * (1 / (1 + exp(-theta[3]))) - 0.02/2
  B2 <- 1.8e-05#0.02 * (1 / (1 + exp(-theta[4]))) - 0.02/2
  C2_coef <- theta[nb1 + 1:nb2]

  D1 <- 0
  D2 <- 0

  if(abs(theta[1]) > 10){
    #return(Inf)
  }

  C1_param <- basis1 %*% matrix(C1_coef, ncol = 1)
  C2_param <- basis2 %*% matrix(C2_coef, ncol = 1)

  cov_mat <- cov_bi_differential(location = location, beta = BETA,
                                 scale_horizontal = SCALE_HORIZONTAL, scale_vertical = SCALE_VERTICAL,
                                 a1 = A1, b1 = B1, c1 = C1_param, d1 = D1, a2 = A2, b2 = B2, c2 = C2_param, d2 = D2,
                                 radius = radius)

  covariance1_marginals <- cov_mat_marginals[1:nrow(location), 1:nrow(location)]
  covariance2_marginals <- cov_mat_marginals[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)]

  covariance1 <- cov_mat[1:nrow(location), 1:nrow(location)]
  covariance2 <- cov_mat[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)]
  covariance12 <- cov_mat[1:nrow(location), nrow(location) + 1:nrow(location)]
  correlation12 <- diag(covariance12) / sqrt(diag(covariance1) * diag(covariance2))

  emp_covariance1 <- empirical_values[1:nrow(location), 1:nrow(location)]
  emp_covariance2 <- empirical_values[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)]
  emp_covariance12 <- empirical_values[1:nrow(location), nrow(location) + 1:nrow(location)]
  emp_correlation12 <- diag(emp_covariance12) / sqrt(diag(emp_covariance1) * diag(emp_covariance2))

  plot_bi_differential(location = location[1:50, ], est_beta = BETA,
                       est_scale_horizontal = SCALE_HORIZONTAL, est_scale_vertical = SCALE_VERTICAL,
                       est_a1 = A1, est_b1 = B1, est_c1 = C1_coef, est_d1 = D1, est_a2 = A2, est_b2 = B2, est_c2 = C2_coef, est_d2 = D2,
                       basis1 = basis1[1:50, ], basis2 = basis2[1:50, ], radius = radius)

  out <- sum((diag(covariance1_marginals) - diag(covariance1))^2 * 100 + (diag(covariance2_marginals) - diag(covariance2))^2 * 50000 + (emp_correlation12 - correlation12)^2 * 100)
  #out <- sum((diag(covariance1_marginals) - diag(covariance1))^2 * 0 + (diag(covariance2_marginals) - diag(covariance2))^2 * 0 + (emp_correlation12 - correlation12)^2 * 100)
  #out <- sum((diag(covariance1_marginals) - diag(covariance1))^2 * 100 + (diag(covariance2_marginals) - diag(covariance2))^2 * 50000 + (emp_correlation12 - correlation12)^2 * 10000)
  #out <- sum((diag(covariance1_marginals) - diag(covariance1))^2 * 0 + (diag(covariance2_marginals) - diag(covariance2))^2 * 500 + (emp_correlation12 - correlation12)^2 * 10)

  return(out)
  #}
}

NEGLOGLIK6 <- function(theta, empirical_values, basis1, basis2, radius){

  print(paste("theta = c(", paste(round(theta, 8), collapse=","), ")", sep = ''))

  nb1 <- ncol(basis1)
  nb2 <- ncol(basis2)

  BETA <- 1 #start with 0 to get the marginals correctly first, then slowly keep increasing to 1. For B4, go up to the minimum rho that will minimize the WLS

  if(REF_LOC_IND == 1){
    #SCALE_HORIZONTAL <- 0.051631#use the I1 estimate from overleaf
    #SCALE_VERTICAL <- 0.007955
  }else if(REF_LOC_IND == 2){
    #SCALE_HORIZONTAL <- 0.181705#use the I1 estimate from overleaf
    #SCALE_VERTICAL <- 0.016264
  }else if(REF_LOC_IND == 3){
    #SCALE_HORIZONTAL <- 0.220331
    #SCALE_VERTICAL <- 0.009407
  }else if(REF_LOC_IND == 4){
    #SCALE_HORIZONTAL <- 0.02506
    #SCALE_VERTICAL <- 0.005254
  }else if(REF_LOC_IND == 5){
    #SCALE_HORIZONTAL <- 0.021257
    #SCALE_VERTICAL <- 0.047564
  }else if(REF_LOC_IND == 6){
    #SCALE_HORIZONTAL <- 0.076977
    #SCALE_VERTICAL <- 0.013413
  }

  SCALE_HORIZONTAL <- 0.251897#0.181846
  SCALE_VERTICAL <- 0.017519#0.014252 # start with 0.5, then keep decreasing to the true SCALE_VERTICAL value

  A1 <- 0.000108
  B1 <-  6.3e-05
  C1_coef <- theta[1:nb1]

  A2 <- -2e-06
  B2 <- 2e-06
  C2_coef <- theta[nb1 + 1:nb2]

  D1 <- 0
  D2 <- 0

  if(abs(theta[1]) > 10){
    #return(Inf)
  }

  C1_param <- basis1 %*% matrix(C1_coef, ncol = 1)
  C2_param <- basis2 %*% matrix(C2_coef, ncol = 1)

  cov_mat <- cov_bi_differential(location = location, beta = BETA,
                                 scale_horizontal = SCALE_HORIZONTAL, scale_vertical = SCALE_VERTICAL,
                                 a1 = A1, b1 = B1, c1 = C1_param, d1 = D1, a2 = A2, b2 = B2, c2 = C2_param, d2 = D2,
                                 radius = radius)


  covariance1_marginals <- cov_mat_marginals[1:nrow(location), 1:nrow(location)]
  covariance2_marginals <- cov_mat_marginals[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)]

  covariance1 <- cov_mat[1:nrow(location), 1:nrow(location)]
  covariance2 <- cov_mat[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)]
  covariance12 <- cov_mat[1:nrow(location), nrow(location) + 1:nrow(location)]
  correlation12 <- diag(covariance12) / sqrt(diag(covariance1) * diag(covariance2))

  emp_covariance1 <- empirical_values[1:nrow(location), 1:nrow(location)]
  emp_covariance2 <- empirical_values[nrow(location) + 1:nrow(location), nrow(location) + 1:nrow(location)]
  emp_covariance12 <- empirical_values[1:nrow(location), nrow(location) + 1:nrow(location)]
  emp_correlation12 <- diag(emp_covariance12) / sqrt(diag(emp_covariance1) * diag(emp_covariance2))

  plot_bi_differential(location = location[1:50, ], est_beta = BETA,
                       est_scale_horizontal = SCALE_HORIZONTAL, est_scale_vertical = SCALE_VERTICAL,
                       est_a1 = A1, est_b1 = B1, est_c1 = C1_coef, est_d1 = D1, est_a2 = A2, est_b2 = B2, est_c2 = C2_coef, est_d2 = D2,
                       basis1 = basis1[1:50, ], basis2 = basis2[1:50, ], radius = radius)

  out <- sum((diag(covariance1_marginals) - diag(covariance1))^2 * 100 + (diag(covariance2_marginals) - diag(covariance2))^2 * 50000 + (emp_correlation12 - correlation12)^2 * 10000)
  #out <- sum((diag(covariance1_marginals) - diag(covariance1))^2 * 100 + (diag(covariance2_marginals) - diag(covariance2))^2 * 50000 + (emp_correlation12 - correlation12)^2 * 10)
  #out <- sum((diag(covariance1_marginals) - diag(covariance1))^2 * 100 + (diag(covariance2_marginals) - diag(covariance2))^2 * 5000000 + (emp_correlation12 - correlation12)^2 * 10)

  return(out)
  #}
}


plot_bi_differential <- function(location, est_beta, est_scale_horizontal, est_scale_vertical, est_a1, est_b1, est_c1, est_d1 = NULL, est_a2, est_b2, est_c2, est_d2 = NULL, radius, basis1, nb1 = ncol(basis1), basis2, nb2 = ncol(basis2), splines_degree = 4, emp_covariance1, emp_covariance2, emp_correlation12){

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
  lines(diag(emp_covariance1), location[, 3], lwd = 2, col = 2)

  plot(variance2, location[, 3], type = 'l', xlab = 'Variance', ylab = 'Depth', ylim = c(max(location[, 3]), 0))
  lines(diag(emp_covariance2), location[, 3], lwd = 2, col = 2)

  plot(covariance12 / sqrt(variance1 * variance2), location[, 3], type = 'l', xlab = 'Colocated Correlation', ylab = 'Depth', ylim = c(max(location[, 3]), 0))
  lines(emp_correlation12, location[, 3], lwd = 2, col = 2)

}


