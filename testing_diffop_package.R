library(DiffOp, lib.loc = "/project/jun/msalvana/R/x86_64-pc-linux-gnu-library/4.2/")
#library(DiffOp)

SYNTHETIC_DATA = F
PLOTTING = F
PREDICTION = F
STEP = 1

MODEL = 'B4'

if(SYNTHETIC_DATA){

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

  if(PLOTTING){
    subset <- which(loc3d[, 3] == 0)

    library(fields)

    pdf(file = paste('/Users/laisalvana/Library/CloudStorage/GoogleDrive-yourlainess@gmail.com/My Drive/Work/Research/UH/bivariate_argo/figures/ex1j.pdf', sep = ''), width = 6.5, height = 10)

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

    dev.off()

    subset <- which(loc3d[, 2] == 0)

    pdf(file = paste('/Users/laisalvana/Library/CloudStorage/GoogleDrive-yourlainess@gmail.com/My Drive/Work/Research/UH/bivariate_argo/figures/ex1k.pdf', sep = ''), width = 6.5, height = 10)

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

    dev.off()

    subset <- which(loc3d[, 1] == 0)

    pdf(file = paste('/Users/laisalvana/Library/CloudStorage/GoogleDrive-yourlainess@gmail.com/My Drive/Work/Research/UH/bivariate_argo/figures/ex1l.pdf', sep = ''), width = 6.5, height = 10)

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

    dev.off()
  }

  ####### ESTIMATION #######

  INIT_BETA = 0
  INIT_SCALE_HORIZONTAL = log(0.02)
  INIT_SCALE_VERTICAL = log(0.2)
  INIT_A1 = INIT_A2 = 0
  INIT_B1 = INIT_B2 = 0
  INIT_D1 = INIT_D2 = 0

  INNER_KNOTS1 <- c(0.1, 0.5, 0.9)
  INNER_KNOTS2 <- c(0.1, 0.5, 0.9)

  SPLINES_DEGREE = 2

  set.seed(1235)
  INIT_C1_COEF <- runif(length(INNER_KNOTS1) + SPLINES_DEGREE + 1, -0.1, 0.1)

  set.seed(1236)
  INIT_C2_COEF <- runif(length(INNER_KNOTS2) + SPLINES_DEGREE + 1, -0.1, 0.1)

  if(STEP == 1){

    est_params_mle <- est_bi_differential_mle(residuals = Z, location = loc3d,
                                              init_beta = INIT_BETA,
                                              init_scale_horizontal = INIT_SCALE_HORIZONTAL,
                                              init_scale_vertical = INIT_SCALE_VERTICAL,
                                              init_a1 = INIT_A1, init_b1 = INIT_B1,
                                              init_c1_coef = INIT_C1_COEF, init_d1 = INIT_D1,
                                              init_a2 = INIT_A2, init_b2 = INIT_B2,
                                              init_c2_coef = INIT_C2_COEF, init_d2 = INIT_D2,
                                              a1_scaling = 1e-3, b1_scaling = 1e-3,
                                              a2_scaling = 1e-3, b2_scaling = 1e-3,
                                              d1_fix = TRUE, d2_fix = TRUE,
                                              radius = earthRadiusKm,
                                              splines_degree = SPLINES_DEGREE,
                                              inner_knots1 = INNER_KNOTS1,
                                              inner_knots2 = INNER_KNOTS2,
                                              iterlim = 5, stepmax = 1, hessian = T)

  }else if(STEP == 2){

    theta = c(-0.61533429,-4.15832011,-1.38957984,-3.80559365,0.01252671,-6.40473156,-5.95366813,-0.73421591,2.82499177,-0.72270707,-1.26147224,-0.04835948,-0.02282755,-2.13890809,-1.63745121,2.73516846,-0.26072885,-3.13790731,-2.53402997)

    start_time = Sys.time()

    est_params_mle2 <- est_bi_differential_mle(residuals = Z, location = loc3d,
                                               init_beta = theta[1],
                                               init_scale_horizontal = theta[2],
                                               init_scale_vertical = theta[3],
                                               init_a1 = theta[4], init_b1 = theta[5],
                                               init_c1_coef = theta[5 + 1:length(INIT_C1_COEF)], init_d1 = 0,
                                               init_a2 = theta[5 + length(INIT_C1_COEF) + 1], init_b2 = theta[5 + length(INIT_C1_COEF) + 2],
                                               init_c2_coef = theta[5 + length(INIT_C1_COEF) + 2 + 1:length(INIT_C2_COEF)], init_d2 = 0,
                                               a1_scaling = 1e-3, b1_scaling = 1e-3,
                                               a2_scaling = 1e-3, b2_scaling = 1e-3,
                                               c1_fix = TRUE, c2_fix = TRUE, d1_fix = TRUE, d2_fix = TRUE,
                                               radius = earthRadiusKm,
                                               splines_degree = SPLINES_DEGREE,
                                               inner_knots1 = INNER_KNOTS1,
                                               inner_knots2 = INNER_KNOTS2,
                                               iterlim = 1000, stepmax = 1, hessian = T)

    print(est_params_mle2)

    end_time = Sys.time()

    TOTAL_TIME <- as.numeric(end_time - start_time, units = "secs")

    print(TOTAL_TIME)

  }


  ####### PREDICTION #######

  if(PREDICTION){
    x_new <- x[-length(x)] + (x[-1] - x[-length(x)]) / 2
    loc2d_new1 <- expand.grid(x_new, y) %>% as.matrix()

    y_new <- y[-length(y)] + (y[-1] - y[-length(y)]) / 2
    loc2d_new2 <- expand.grid(c(x, x_new), y_new) %>% as.matrix()

    loc2d_new <- rbind(loc2d_new1, loc2d_new2)

    loc3d_new <- cbind(rep(loc2d_new[, 1], each = length(depth)), rep(loc2d_new[, 2], each = length(depth)), depth)

    Z_pred <- predict_bi_differential(residuals = Z, location = loc3d,
                                      location_new = loc3d_new, est_beta = 0.345513,
                                      est_scale_horizontal = 0.01541487,
                                      est_scale_vertical = 0.2473264,
                                      est_a1 = 2.534356e-05, est_b1 = -4.751549e-06,
                                      est_c1_coef = c(-6.0423734, -5.8585368, -0.7910191, 2.9407869, -0.7534488, -1.3051477), est_d1 = 0,
                                      est_a2 = -4.914751e-05, est_b2 = -2.395199e-05,
                                      est_c2_coef = c(-2.1237812, -1.7196022, 2.9021517, -0.3172772, -3.0745846, -2.5363932), est_d2 = 0,
                                      radius = earthRadiusKm, splines_degree = SPLINES_DEGREE,
                                      inner_knots1 = INNER_KNOTS1, inner_knots2 = INNER_KNOTS2)

    Z1_pred <- Z_pred[1:nrow(loc3d_new)]
    Z2_pred <- Z_pred[nrow(loc3d_new) + 1:nrow(loc3d_new)]

    loc3d_full <- rbind(loc3d, loc3d_new)
    Z1_full <- c(Z1, Z1_pred)
    Z2_full <- c(Z2, Z2_pred)

    Z_full <- c(Z, Z_pred)

    if(PLOTTING){

      library(fields)

      subset <- which(loc3d_full[, 3] == 0)

      pdf(file = paste('/Users/laisalvana/Library/CloudStorage/GoogleDrive-yourlainess@gmail.com/My Drive/Work/Research/UH/bivariate_argo/figures/ex1j_new.pdf', sep = ''), width = 6.5, height = 10)

      par(mfrow = c(2, 1))
      par(pty = 's')
      par(mai = c(0.3, 0.1, 0.6, 0.5))
      quilt.plot(loc3d_full[subset, 1], loc3d_full[subset, 2], Z1_full[subset], nx = length(x) + length(x_new), ny = length(y) + length(y_new), xlab = "Longitude", ylab = "", zlim = range(Z_full), xaxt = 'n')
      mtext("Latitude", side = 2, line = 2)
      mtext(expression(Z[1]), side = 2, line = 3, col = 4, cex = 1.5)
      par(mai = c(0.8, 0.1, 0.1, 0.5))
      quilt.plot(loc3d_full[subset, 1], loc3d_full[subset, 2], Z2_full[subset], nx = length(x) + length(x_new), ny = length(y) + length(y_new), xlab = "Longitude", ylab = "", zlim = range(Z_full))
      mtext("Latitude", side = 2, line = 2)
      mtext(expression(Z[2]), side = 2, line = 3, col = 4, cex = 1.5)

      dev.off()

      subset <- which(loc3d_full[, 2] == 0)

      pdf(file = paste('/Users/laisalvana/Library/CloudStorage/GoogleDrive-yourlainess@gmail.com/My Drive/Work/Research/UH/bivariate_argo/figures/ex1k_new.pdf', sep = ''), width = 6.5, height = 10)

      par(mfrow = c(2, 1))
      par(pty = 's')
      par(mai = c(0.3, 0.1, 0.6, 0.5))
      quilt.plot(loc3d_full[subset, 1], loc3d_full[subset, 3], Z1_full[subset], nx = length(x) + length(x_new), ny = length(depth), xlab = "Longitude", ylab = "", zlim = range(Z_full), ylim = c(1, 0), xaxt = 'n')
      mtext("Depth", side = 2, line = 2)
      mtext(expression(Z[1]), side = 2, line = 3, col = 4, cex = 1.5)
      par(mai = c(0.8, 0.1, 0.1, 0.5))
      quilt.plot(loc3d_full[subset, 1], loc3d_full[subset, 3], Z2_full[subset], nx = length(x) + length(x_new), ny = length(depth), xlab = "Longitude", ylab = "", zlim = range(Z_full), ylim = c(1, 0))
      mtext("Depth", side = 2, line = 2)
      mtext(expression(Z[2]), side = 2, line = 3, col = 4, cex = 1.5)

      dev.off()

      subset <- which(loc3d_full[, 1] == 0)

      pdf(file = paste('/Users/laisalvana/Library/CloudStorage/GoogleDrive-yourlainess@gmail.com/My Drive/Work/Research/UH/bivariate_argo/figures/ex1l_new.pdf', sep = ''), width = 6.5, height = 10)

      par(mfrow = c(2, 1))
      par(pty = 's')
      par(mai = c(0.3, 0.1, 0.6, 0.5))
      quilt.plot(loc3d_full[subset, 2], loc3d_full[subset, 3], Z1_full[subset], nx = length(y) + length(y_new), ny = length(depth), xlab = "Latitude", ylab = "", zlim = range(Z_full), ylim = c(1, 0), xaxt = 'n')
      mtext("Depth", side = 2, line = 2)
      mtext(expression(Z[1]), side = 2, line = 3, col = 4, cex = 1.5)
      par(mai = c(0.8, 0.1, 0.1, 0.5))
      quilt.plot(loc3d_full[subset, 2], loc3d_full[subset, 3], Z2_full[subset], nx = length(y) + length(y_new), ny = length(depth), xlab = "Latitude", ylab = "", zlim = range(Z_full), ylim = c(1, 0))
      mtext("Depth", side = 2, line = 2)
      mtext(expression(Z[2]), side = 2, line = 3, col = 4, cex = 1.5)

      dev.off()
    }
  }

}else{

  data("argo_ref_loc1")

  profile_no_unique <- unique(argo_ref_loc1$ProfileNumber)

  profile_no_closest <- which(argo_ref_loc1$DistanceKmFromRefLoc == min(argo_ref_loc1$DistanceKmFromRefLoc))

  PLOTTING = F

  if(PLOTTING){

    plot(0, 0, xlim = range(argo_ref_loc1$TemperatureResiduals), ylim = c(2000, 0), xlab = 'Temperature Residuals', ylab = 'Depth (meters)', type = 'n')

    for(prof in profile_no_unique){
      subset <- which(argo_ref_loc1$ProfileNumber == prof)
      lines(argo_ref_loc1$TemperatureResiduals[subset], argo_ref_loc1$Pressure[subset], lwd = 0.5, col = 'gray')
    }

    lines(argo_ref_loc1$TemperatureResiduals[profile_no_closest], argo_ref_loc1$Pressure[profile_no_closest], lwd = 1.5, col = 'black')

    plot(0, 0, xlim = range(argo_ref_loc1$SalinityResiduals), ylim = c(2000, 0), xlab = 'Salinity Residuals', ylab = 'Depth (meters)', type = 'n')

    for(prof in profile_no_unique){
      subset <- which(argo_ref_loc1$ProfileNumber == prof)
      lines(argo_ref_loc1$SalinityResiduals[subset], argo_ref_loc1$Pressure[subset], lwd = 0.5, col = 'gray')
    }

    lines(argo_ref_loc1$SalinityResiduals[profile_no_closest], argo_ref_loc1$Pressure[profile_no_closest], lwd = 1.5, col = 'black')

  }

  #ind_pred <- profile_no_closest
  ind_pred <- c(1:50, seq(1051, 2500, by = 1))

  loc3d <- cbind(argo_ref_loc1$Longitude, argo_ref_loc1$Latitude, argo_ref_loc1$Pressure)
  locs_insample <- loc3d[-ind_pred, ]
  locs_outsample <- loc3d[ind_pred, ]

  Z_insample <- c(argo_ref_loc1$TemperatureResiduals[-ind_pred], argo_ref_loc1$SalinityResiduals[-ind_pred])
  Z_outsample <- c(argo_ref_loc1$TemperatureResiduals[ind_pred], argo_ref_loc1$SalinityResiduals[ind_pred])

  earthRadiusKm = 6371

  emp_cov <- compute_emp_cov(location = locs_insample,
                             variable1_residuals = Z_insample[1:nrow(locs_insample)],
                             variable2_residuals = Z_insample[nrow(locs_insample) + 1:nrow(locs_insample)],
                             bandwidth_horizontal = 0.009, bandwidth_vertical = 0.03,
                             radius = earthRadiusKm)

  #emp_variance1 <- diag(emp_cov[1:nrow(locs_insample), 1:nrow(locs_insample)])
  #emp_variance2 <- diag(emp_cov[nrow(locs_insample) + 1:nrow(locs_insample), nrow(locs_insample) + 1:nrow(locs_insample)])
  #emp_covariance12 <- diag(emp_cov[1:nrow(locs_insample), nrow(locs_insample) + 1:nrow(locs_insample)])
  #emp_correlation12 <- emp_covariance12 / sqrt(emp_variance1 * emp_variance2)

  INNER_KNOTS1 <- c(50, 100, 300, 500, 700, 1000)
  INNER_KNOTS2 <- c(50, 100, 300, 500, 700, 1000)

  SPLINES_DEGREE <- 2

  INIT_BETA = 0
  INIT_SCALE_HORIZONTAL = -5.7331239134
  INIT_SCALE_VERTICAL = 2.564949 #0.26357412
  INIT_A1 = 4.07291428
  INIT_B1 = 3.59165191
  INIT_C1_COEF = c(-0.87653831,0.50746214,0.57719275,0.61521185,0.37327857,0.14017077,-0.02700237,0.0004893,0.071279) * 1e-2
  INIT_A2 = 2.55032082
  INIT_B2 = -2.54202768
  INIT_C2_COEF = c(0.04952738,0.22467233,-0.16542861,-0.04447287,-0.0159132,0.00096767,-0.00318231,0.00404039,-0.00412418) * 1e-1
  INIT_D1 = INIT_D2 = 0

  #est_params_wls <- est_bi_differential_wls(empirical_values = emp_cov, location = locs_insample,
  #                                          init_beta = INIT_BETA,
  #                                          init_scale_horizontal = INIT_SCALE_HORIZONTAL,
  #                                          init_scale_vertical = INIT_SCALE_VERTICAL,
  #                                          init_a1 = INIT_A1, init_b1 = INIT_B1,
  #                                          init_c1_coef = INIT_C1_COEF, init_d1 = INIT_D1,
  #                                          init_a2 = INIT_A2, init_b2 = INIT_B2,
  #                                          init_c2_coef = INIT_C2_COEF, init_d2 = INIT_D2,
  #                                          a1_scaling = 1e-3, b1_scaling = 1e-3,
  #                                          a2_scaling = 1e-3, b2_scaling = 1e-3,
  #                                          c1_coef_scaling = 100, c2_coef_scaling = 10,
  #                                          d1_fix = TRUE, d2_fix = TRUE,
  #                                          radius = earthRadiusKm,
  #                                          splines_degree = SPLINES_DEGREE,
  #                                          inner_knots1 = INNER_KNOTS1, inner_knots2 = INNER_KNOTS2,
  #                                          w1 = 100, w2 = 50000, w12 = 10,
  #                                          iterlim = 1000, stepmax = 1, hessian = T)


  INNER_KNOTS1 <- c(100, 500, 1000)
  INNER_KNOTS2 <- c(100, 500, 1000)

  SPLINES_DEGREE = 2

  set.seed(1235)
  INIT_C1_COEF <- runif(length(INNER_KNOTS1) + SPLINES_DEGREE + 1, -0.1, 0.1)

  set.seed(1236)
  INIT_C2_COEF <- runif(length(INNER_KNOTS2) + SPLINES_DEGREE + 1, -0.1, 0.1)

  if(MODEL == 'B4'){

    if(STEP == 1){
      est_params_mle <- est_bi_differential_mle(residuals = Z_insample, location = locs_insample,
                                                init_beta = INIT_BETA,
                                                init_scale_horizontal = INIT_SCALE_HORIZONTAL,
                                                init_scale_vertical = INIT_SCALE_VERTICAL,
                                                init_a1 = INIT_A1, init_b1 = INIT_B1,
                                                init_c1_coef = INIT_C1_COEF, init_d1 = INIT_D1,
                                                init_a2 = INIT_A2, init_b2 = INIT_B2,
                                                init_c2_coef = INIT_C2_COEF, init_d2 = INIT_D2,
                                                a1_scaling = 1e-3, b1_scaling = 1e-3,
                                                a2_scaling = 1e-3, b2_scaling = 1e-3,
                                                c1_coef_scaling = 1, c2_coef_scaling = 1,
                                                d1_fix = TRUE, d2_fix = TRUE,
                                                radius = earthRadiusKm,
                                                splines_degree = SPLINES_DEGREE,
                                                inner_knots1 = INNER_KNOTS1, inner_knots2 = INNER_KNOTS2,
                                                iterlim = 1000, stepmax = 1, hessian = T)

      print(est_params_mle)

    }else if(STEP == 2){

      est_params_mle2 <- est_bi_differential_mle(residuals = Z_insample, location = locs_insample,
                                                 init_beta = theta[1],
                                                 init_scale_horizontal = theta[2],
                                                 init_scale_vertical = theta[3],
                                                 init_a1 = theta[4], init_b1 = theta[5],
                                                 init_c1_coef = theta[5 + 1:length(INIT_C1_COEF)], init_d1 = 0,
                                                 init_a2 = theta[5 + length(INIT_C1_COEF) + 1], init_b2 = theta[5 + length(INIT_C1_COEF) + 2],
                                                 init_c2_coef = theta[5 + length(INIT_C1_COEF) + 2 + 1:length(INIT_C2_COEF)], init_d2 = 0,
                                                 a1_scaling = 1e-3, b1_scaling = 1e-3,
                                                 a2_scaling = 1e-3, b2_scaling = 1e-3,
                                                 c1_coef_scaling = 100, c2_coef_scaling = 10,
                                                 d1_fix = TRUE, d2_fix = TRUE,
                                                 radius = earthRadiusKm,
                                                 splines_degree = SPLINES_DEGREE,
                                                 inner_knots1 = INNER_KNOTS1,
                                                 inner_knots2 = INNER_KNOTS2,
                                                 iterlim = 1000, stepmax = 1, hessian = T)

      print(est_params_mle2)

    }

  }else if(MODEL == 'B3'){

    if(STEP == 1){

      start_time = Sys.time()

      est_params_mle <- est_bi_differential_mle(residuals = Z_insample, location = locs_insample,
                                                init_beta = 1,
                                                init_scale_horizontal = exp(-4.5),
                                                init_scale_vertical = exp(-0.12490163),
                                                init_a1 = INIT_A1, init_b1 = INIT_B1,
                                                init_c1_coef = INIT_C1_COEF, init_d1 = INIT_D1,
                                                init_a2 = INIT_A2, init_b2 = INIT_B2,
                                                init_c2_coef = INIT_C2_COEF, init_d2 = INIT_D2,
                                                vertical_scale_scaling = 1e-1,
                                                a1_scaling = 1e-3, b1_scaling = 1e-3,
                                                a2_scaling = 1e-3, b2_scaling = 1e-3,
                                                c1_coef_scaling = 1, c2_coef_scaling = 1,
                                                beta_fix = TRUE, scale_horizontal_fix = TRUE,
                                                scale_vertical_fix = TRUE,
                                                d1_fix = TRUE, d2_fix = TRUE,
                                                radius = earthRadiusKm,
                                                splines_degree = SPLINES_DEGREE,
                                                inner_knots1 = INNER_KNOTS1, inner_knots2 = INNER_KNOTS2,
                                                iterlim = 1000, stepmax = 1, hessian = T)
      print(est_params_mle)

      end_time = Sys.time()

      TOTAL_TIME <- as.numeric(end_time - start_time, units = "secs")

      print(TOTAL_TIME)

    }else if(STEP == 2){

      theta = c(-4.5, -0.12490163, 0.87215125,-1.23722843,0.64253503,0.91872335,0.77176668,-0.01379812,0.08096983,-0.01675213,-0.52312691,0.37771119,0.13158042,0.19086266,0.02772741,0.01351954,-0.00075524,0.00723768)

      est_params_mle2 <- est_bi_differential_mle(residuals = Z_insample, location = locs_insample,
                                                 init_beta = 0,
                                                 init_scale_horizontal = theta[1],
                                                 init_scale_vertical = theta[2],
                                                 init_a1 = theta[3], init_b1 = theta[4],
                                                 init_c1_coef = theta[4 + 1:length(INIT_C1_COEF)], init_d1 = 0,
                                                 init_a2 = theta[4 + length(INIT_C1_COEF) + 1], init_b2 = theta[4 + length(INIT_C1_COEF) + 2],
                                                 init_c2_coef = theta[4 + length(INIT_C1_COEF) + 2 + 1:length(INIT_C2_COEF)], init_d2 = 0,
                                                 a1_scaling = 1e-3, b1_scaling = 1e-3,
                                                 a2_scaling = 1e-3, b2_scaling = 1e-3,
                                                 beta_fix = TRUE, d1_fix = TRUE, d2_fix = TRUE,
                                                 radius = earthRadiusKm,
                                                 splines_degree = SPLINES_DEGREE,
                                                 inner_knots1 = INNER_KNOTS1,
                                                 inner_knots2 = INNER_KNOTS2,
                                                 iterlim = 1000, stepmax = 1, hessian = T)
      print(est_params_mle2)
    }
  }

  if(PLOTTING){

    mle_est_cov_mat <- cov_bi_differential(location = locs_insample,
                                           beta = est_params_mle2$est_beta,
                                           scale_horizontal = est_params_mle2$est_scale_horizontal,
                                           scale_vertical = est_params_mle2$est_scale_vertical,
                                           a1 = est_params_mle2$est_a1,
                                           b1 = est_params_mle2$est_b1,
                                           c1_coef = est_params_mle2$est_c1_coef,
                                           d1 = est_params_mle2$est_d1,
                                           a2 = est_params_mle2$est_a2,
                                           b2 = est_params_mle2$est_b2,
                                           c2_coef = est_params_mle2$est_c2_coef,
                                           d2 = est_params_mle2$est_d2,
                                           radius = earthRadiusKm,
                                           splines_degree = est_params_mle2$splines_degree,
                                           inner_knots1 = est_params_mle2$inner_knots1,
                                           inner_knots2 = est_params_mle2$inner_knots2)

    mle_est_variance1 <- diag(mle_est_cov_mat[1:nrow(locs_insample), 1:nrow(locs_insample)])
    mle_est_variance2 <- diag(mle_est_cov_mat[nrow(locs_insample) + 1:nrow(locs_insample), nrow(locs_insample) + 1:nrow(locs_insample)])
    mle_est_covariance12 <- diag(mle_est_cov_mat[1:nrow(locs_insample), nrow(locs_insample) + 1:nrow(locs_insample)])
    mle_est_correlation12 <- mle_est_covariance12 / sqrt(mle_est_variance1 * mle_est_variance2)

    par(mfrow = c(2, 3))

    plot(mle_est_variance1[1:50])
    plot(mle_est_variance2[1:50])
    plot(mle_est_correlation12[1:50])

  }
}

####### SCRATCH #######

SCRATCH = F

if(SCRATCH){

  par(mfrow = c(2, 3))

  plot(variance1[1:10])
  plot(variance2[1:10])
  plot(correlation12[1:10])

  plot(emp_variance1[1:10])
  plot(emp_variance2[1:10])
  plot(emp_correlation12[1:10])

  plot(emp_cov[1, seq(1, 1000, by = 10)])

  INIT_BETA = 0
  INIT_SCALE_HORIZONTAL = log(0.02)
  INIT_SCALE_VERTICAL = log(0.2)
  INIT_A1 = INIT_A2 = 0
  INIT_B1 = INIT_B2 = 0
  INIT_D1 = INIT_D2 = 0
  INIT_C1_COEF = INIT_C2_COEF = 1

  est_params_mle <- est_bi_differential_mle(residuals = Z, location = loc3d,
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

  #theta = c(-4.65560477,-10.57604634,13.70004287,25.78831257,3.49253491,4.66280469,3.15236873,-0.63364572,-4.00411605,-5.60298597,-4.6558888,-5.21335133)
  #theta = c(-4.05506591,-3.1916669,-1.59170792,5.46825639,4.9853047,3.94290889,4.39038617,2.5052013,-5.33414962,-5.07497427,-4.92331565,-4.68664234)
  #negloglik:  -3669.4403

  mle_est_cov_mat <- cov_bi_differential(location = loc3d, beta = est_params_mle$est_beta, scale_horizontal = est_params_mle$est_scale_horizontal, scale_vertical = est_params_mle$est_scale_vertical, a1 = est_params_mle$est_a1, b1 = est_params_mle$est_b1, c1_coef = est_params_mle$est_c1_coef, d1 = 0,
                                         a2 = est_params_mle$est_a2, b2 = est_params_mle$est_b2, c2_coef = est_params_mle$est_c2_coef, d2 = 0, radius = earthRadiusKm, splines_degree = SPLINES_DEGREE, knots1 = KNOTS1, knots2 = KNOTS2)

  theta = c(-0.536965721, 8.515860565, 10.760186134, -3.307735726, -3.799200253, -1.717745441,
            4.614814480, 3.955609777, 4.351554653, -0.004724896, -0.016816880, 4.565644420,
            3.274769413, -5.380340068, -2.490770406, -4.280324693, -3.357557703)

  est_params_mle2 <- est_bi_differential_mle(residuals = Z, location = loc3d,
                                             init_beta = theta[1],
                                             init_scale_horizontal = -5.09094054,
                                             init_scale_vertical = -1.34307397,
                                             init_a1 = theta[2], init_b1 = theta[3],
                                             init_c1_coef = theta[3 + 1:length(INIT_C1_COEF)], init_d1 = 0,
                                             init_a2 = theta[3 + length(INIT_C1_COEF) + 1], init_b2 = theta[3 + length(INIT_C1_COEF) + 2],
                                             init_c2_coef = theta[3 + length(INIT_C1_COEF) + 2 + 1:length(INIT_C2_COEF)], init_d2 = 0,
                                             d1_fix = TRUE, d2_fix = TRUE,
                                             radius = earthRadiusKm,
                                             splines_degree = SPLINES_DEGREE,
                                             knots1 = KNOTS1, knots2 = KNOTS2,
                                             iterlim = 1000, stepmax = 1, hessian = T)

  theta = c(-0.76204644,-4.55804423,-1.94605509,0.12655982,0.84346228,-1.77218938,-2.58028244,-4.76804631,11.30294595,19.61923519,13.42143752,0.01686747,-0.0091102,4.45767872,3.74930212,-5.95174338,-13.13890413,-6.97093002,-5.27942716)
  theta = c(-0.6632838,-4.68657405,-1.96919145,-0.0181081,-0.00779444,-2.29079773,-3.37330351,-6.36289091,15.65553785,17.68589565,15.42971718,0.01913544,-0.02844181,4.50154414,3.6582661,-5.87692042,-10.13408209,-6.36045902,-4.62118685)
  theta = c(-0.66590267,-4.79841894,-2.08258705,-0.02122895,-0.00879531,-2.86779729,-4.22278888,-7.97421044,19.62227991,22.18458008,19.31297487,0.02191214,-0.03393269,5.63829784,4.58674391,-7.37019452,-12.71733809,-7.97630343,-5.79275871)

  est_params_mle2 <- est_bi_differential_mle(residuals = Z, location = loc3d,
                                             init_beta = theta[1],
                                             init_scale_horizontal = theta[2],
                                             init_scale_vertical = theta[3],
                                             init_a1 = theta[4], init_b1 = theta[5],
                                             init_c1_coef = theta[5 + 1:length(INIT_C1_COEF)], init_d1 = 0,
                                             init_a2 = theta[5 + length(INIT_C1_COEF) + 1], init_b2 = theta[5 + length(INIT_C1_COEF) + 2],
                                             init_c2_coef = theta[5 + length(INIT_C1_COEF) + 2 + 1:length(INIT_C2_COEF)], init_d2 = 0,
                                             d1_fix = TRUE, d2_fix = TRUE,
                                             radius = earthRadiusKm,
                                             splines_degree = SPLINES_DEGREE,
                                             knots1 = KNOTS1, knots2 = KNOTS2,
                                             iterlim = 1000, stepmax = 1, hessian = T)

  mle_est_cov_mat <- cov_bi_differential(location = loc3d, beta = 0.1329132, scale_horizontal = 0.01393392, scale_vertical = 0.2279596, a1 = exp(-3.35678979) * 1e-3, b1 = 0.01580268 * 1e-3, c1_coef = c(-6.13353992,-6.45820406,-1.0103012,3.56694107,-0.91411356,-1.59157996), d1 = 0,
                                         a2 = -0.05821896 * 1e-3, b2 = -0.02910226 * 1e-3, c2_coef = c(-2.58166409,-2.03010153,3.41535242,-0.38895313,-3.52456627,-2.8773148), d2 = 0, radius = earthRadiusKm, splines_degree = SPLINES_DEGREE, knots1 = KNOTS1, knots2 = KNOTS2)


  theta2 = c(-0.54476181,-4.84498315,-1.43957751,5.17121931,6.46133916,0.03203775,0.21407682)
  theta2 = c(0.02620172,-4.74685329,-1.48004991,1.00673748,1.88466654,0.0007775,0.06572587)
  theta2 = c(-0.30003348,-4.64710744,-1.46041448,-0.12023626,0.07423777,-0.01005004,0.00749788)

  est_params_mle_fin <- est_bi_differential_mle(residuals = Z, location = loc3d,
                                                init_beta = theta2[1],
                                                init_scale_horizontal = theta2[2],
                                                init_scale_vertical = theta2[3],
                                                init_a1 = theta2[4], init_b1 = theta2[5],
                                                init_c1_coef = theta[3 + 1:length(INIT_C1_COEF)], init_d1 = 0,
                                                init_a2 = theta2[6], init_b2 = theta2[7],
                                                init_c2_coef = theta[3 + length(INIT_C1_COEF) + 2 + 1:length(INIT_C2_COEF)], init_d2 = 0,
                                                c1_fix = TRUE, c2_fix = TRUE,
                                                d1_fix = TRUE, d2_fix = TRUE,
                                                radius = earthRadiusKm,
                                                splines_degree = SPLINES_DEGREE,
                                                knots1 = KNOTS1, knots2 = KNOTS2,
                                                iterlim = 1000, stepmax = 1, hessian = T)

  ##############################################################

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


}
