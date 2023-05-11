library(devtools)
library(usethis)

usethis::use_r("DAT")

devtools::document()
devtools::check()
devtools::install()

#To generate PDF version of manual
#On the terminal: R CMD check DiffOp

library(DiffOp)
library(dplyr)

x <- seq(0, 1, length.out = 10)
y <- seq(0, 1, length.out = 10)

loc2d <- expand.grid(x, y) %>% as.matrix()
depth <- seq(0, 1, length.out = 10)
loc3d <- cbind(rep(loc2d[, 1], each = length(depth)), rep(loc2d[, 2], each = length(depth)), depth)

earthRadiusKm = 6371

BETA = 0.5
SCALE_HORIZONTAL = 0.03 #0.002998395 #0.01
SCALE_VERTICAL = 0.3 #0.01644395 #0.5
A1 = A2 = 0.00001
B1 = B2 = 0.00001
C1 = sin((loc3d[1:10, 3] + 0.1) * pi / 0.5)
C2 = cos((loc3d[1:10, 3] + 0.1) * pi / 0.5)
D1 = D2 = 0

cov_mat <- cov_bi_differential(location = loc3d, beta = BETA,
                               scale_horizontal = SCALE_HORIZONTAL, scale_vertical = SCALE_VERTICAL,
                               a1 = A1, b1 = B1, c1 = C1, d1 = D1, a2 = A2, b2 = B2, c2 = C2, d2 = D2,
                               radius = earthRadiusKm)

## SIMULATION

library(MASS)

set.seed(1235)
Z <- mvrnorm(1, mu = rep(0, ncol(cov_mat)), Sigma = cov_mat)
Z1 <- Z[1:nrow(loc3d)]
Z2 <- Z[nrow(loc3d) + 1:nrow(loc3d)]

## PLOTTING SIMULATION

library(fields)

subset <- which(loc3d[, 3] == 0)

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

## ESTIMATION

KNOTS1 <- seq(0, max(loc3d[, 3]), length.out = 10)
KNOTS2 <- seq(0, max(loc3d[, 3]), length.out = 10)

basis1 <- bsplineBasis(loc3d[, 3], 2, KNOTS1)
nb1 <- ncol(basis1)
basis2 <- bsplineBasis(loc3d[, 3], 2, KNOTS2)
nb2 <- ncol(basis2)

INIT_BETA = 0.7
INIT_SCALE_HORIZONTAL = 0.002998395 #0.01
INIT_SCALE_VERTICAL = 0.01644395 #0.5
INIT_A1 = INIT_A2 = 0
INIT_B1 = INIT_B2 = 0
INIT_D1 = INIT_D2 = 0

set.seed(1235)
INIT_C1 <- runif(nb1, -5, 5)

set.seed(1236)
INIT_C2 <- runif(nb2, -5, 5)

#Estimate d
est_params <- est_bi_differential(residuals = Z, location = loc3d, init_beta = INIT_BETA,
                                  init_scale_horizontal = INIT_SCALE_HORIZONTAL, init_scale_vertical = INIT_SCALE_VERTICAL,
                                  init_a1 = INIT_A1, init_b1 = INIT_B1, init_c1 = INIT_C1, init_d1 = INIT_D1, init_a2 = INIT_A2, init_b2 = INIT_B2, init_c2 = INIT_C2, init_d2 = INIT_D2,
                                  basis1 = basis1, basis2 = basis2, radius = earthRadiusKm)

#No d
est_params <- est_bi_differential(residuals = Z, location = loc3d, init_beta = INIT_BETA,
                                  init_scale_horizontal = INIT_SCALE_HORIZONTAL, init_scale_vertical = INIT_SCALE_VERTICAL,
                                  init_a1 = INIT_A1, init_b1 = INIT_B1, init_c1 = INIT_C1, init_a2 = INIT_A2, init_b2 = INIT_B2, init_c2 = INIT_C2,
                                  basis1 = basis1, basis2 = basis2, radius = earthRadiusKm)

## REAL DATA

#load("/Users/laisalvana/Library/CloudStorage/GoogleDrive-yourlainess@gmail.com/My Drive/Work/Research/UH/bivariate_argo/data/jan_march_residuals.RData")
#plotSurfaceMap(variable = 'temp')

data("argo_ref_loc1")
earthRadiusKm = 6371

REF_LOC_IND = 1

if(REF_LOC_IND == 1){
  SCALE_HORIZONTAL = 0.002998395
  SCALE_VERTICAL = 0.01644395
}else if(REF_LOC_IND == 2){
  SCALE_HORIZONTAL = 0.4847018
  SCALE_VERTICAL = 0.1075883
}else if(REF_LOC_IND == 3){
  SCALE_HORIZONTAL = 0.07790723
  SCALE_VERTICAL = 0.2272741
}else if(REF_LOC_IND == 4){
  SCALE_HORIZONTAL = 0.003586815
  SCALE_VERTICAL = 0.009939062
}else if(REF_LOC_IND == 5){
  SCALE_HORIZONTAL = 0.2794646
  SCALE_VERTICAL = 0.1496154
}else if(REF_LOC_IND == 6){
  SCALE_HORIZONTAL = 0.3654966
  SCALE_VERTICAL = 0.1015624
}

ind <- 51:1000
ind_pred <- 1:50

locs_insample <- DAT@locs[ind, 1:3]
locs_outsample <- DAT@locs[ind_pred, 1:3]
loc3d <- rbind(locs_insample, locs_outsample)

Z_insample <- cbind(DAT@measurements[ind, ], DAT@measurements[2000 + ind, ])
Z_outsample <- cbind(DAT@measurements[ind_pred, ], DAT@measurements[2000 + ind_pred, ])
Z <- c(rbind(Z_insample, Z_outsample))

KNOTS1 <- seq(0, max(loc3d[, 3]), length.out = 3)
KNOTS2 <- seq(0, max(loc3d[, 3]), length.out = 3)

basis1 <- bsplineBasis(loc3d[, 3], 2, KNOTS1)
nb1 <- ncol(basis1)
basis2 <- bsplineBasis(loc3d[, 3], 2, KNOTS2)
nb2 <- ncol(basis2)

INIT_BETA = 0.7 #try 0, 0.7
INIT_SCALE_HORIZONTAL = log(SCALE_HORIZONTAL)
INIT_SCALE_VERTICAL = log(SCALE_VERTICAL)
INIT_A1 = INIT_A2 = 0 #0.001
INIT_B1 = INIT_B2 = 0 #0.001
INIT_D1 = INIT_D2 = 0

set.seed(1235)
INIT_C1 <- runif(nb1, -5, 5)

set.seed(1236)
INIT_C2 <- runif(nb2, -5, 5)

est_params <- est_bi_differential(residuals = Z, location = loc3d, init_beta = INIT_BETA,
                                  init_scale_horizontal = INIT_SCALE_HORIZONTAL, init_scale_vertical = INIT_SCALE_VERTICAL, init_scale_horizontal_fix = F, init_scale_vertical_fix = F,
                                  init_a1 = INIT_A1, init_b1 = INIT_B1, init_c1 = INIT_C1, init_a2 = INIT_A2, init_b2 = INIT_B2, init_c2 = INIT_C2,
                                  basis1 = basis1, basis2 = basis2, radius = earthRadiusKm)

################################################################


ref_lat <- c(40, 0, -40, 40, 0, -40)
ref_long <- c(-175, -175, -175, -30, -30, -30)

loc3d_eval <- cbind(ref_long[1], ref_lat[1], seq(0, 2000, length.out = 100))

new_basis1 <- bsplineBasis(loc3d_eval[, 3], 2, KNOTS1)
new_basis2 <- bsplineBasis(loc3d_eval[, 3], 2, KNOTS2)

plot_bi_differential(location = loc3d_eval, est_beta = theta[1],
                    est_scale_horizontal = theta[2], est_scale_vertical = theta[3],
                    est_a1 = theta[4], est_b1 = theta[5], est_c1 = theta[5 + 1:nb1], est_a2 = theta[nb1 + 6], est_b2 = theta[nb1 + 7], est_c2 = theta[nb1 + 7 + 1:nb2],
                    basis1 = basis1, basis2 = basis2, radius = earthRadiusKm)
abline(h = KNOTS1, col = 'brown')

