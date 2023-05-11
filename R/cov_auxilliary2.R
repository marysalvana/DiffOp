NEGLOGLIK <- function(theta, empirical_values, basis1, basis2, radius){

  print(paste("theta = c(", paste(round(theta, 8), collapse=","), ")", sep = ''))

  nb1 <- ncol(basis1)
  nb2 <- ncol(basis2)

  BETA <- 0 #start with 0 to get the marginals correctly first, then slowly keep increasing to 1. For B4, go up to the minimum rho that will minimize the WLS

  SCALE_HORIZONTAL <- 0.021257#0.018255#0.052659 #use the I1 estimate from overleaf
  SCALE_VERTICAL <- 0.047564#0.026259 # start with 0.5, then keep decreasing to 0.024813

  A1 <- 0.01 * (1 / (1 + exp(-theta[1])))
  B1 <- 0.02 * (1 / (1 + exp(-theta[2]))) - 0.02/2
  C1_coef <- theta[2 + 1:nb1]

  A2 <- 0.02 * (1 / (1 + exp(-theta[2 + nb1 + 1]))) - 0.02/2
  B2 <-0.02 * (1 / (1 + exp(-theta[2 + nb1 + 2]))) - 0.02/2
  C2_coef <- theta[2 + nb1 + 2 + 1:nb2]

  #A2 <- 0#2 * (1 / (1 + exp(-theta[2 + nb1 + 1]))) - 1
  #B2 <- 0#2 * (1 / (1 + exp(-theta[2 + nb1 + 2]))) - 1
  #C2_coef <- rep(0.1, nb1)#theta[2 + nb1 + 2 + 1:nb2]

  #A1 <- 0#(1 / (1 + exp(-theta[1])))
  #B1 <- 0#2 * (1 / (1 + exp(-theta[2]))) - 1
  #C1_coef <- rep(0.1, nb2)#theta[2 + 1:nb1]

  #A2 <- 2 * (1 / (1 + exp(-theta[1]))) - 1
  #B2 <- 2 * (1 / (1 + exp(-theta[2]))) - 1
  #C2_coef <- theta[2 + 1:nb2]

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






