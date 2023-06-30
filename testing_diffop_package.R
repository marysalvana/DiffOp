
library(DiffOp)

library(dplyr)


x <- seq(0, 1, length.out = 20)
y <- seq(0, 1, length.out = 20)
loc2d <- expand.grid(x, y) %>% as.matrix()

depth <- seq(0, 1, length.out = 10)
loc3d <- cbind(rep(loc2d[, 1], each = length(depth)), rep(loc2d[, 2], each = length(depth)), depth)

earthRadiusKm = 6371

BETA = 0.5
SCALE_HORIZONTAL = 0.03
SCALE_VERTICAL = 0.3
A1 = A2 = 0.00001
B1 = B2 = 0.00001
C1 = sin((loc3d[1:10, 3] + 0.1) * pi / 0.5)
C2 = cos((loc3d[1:10, 3] + 0.1) * pi / 0.5)
D1 = D2 = 0

#cov_mat <- cov_bi_differential(location = loc3d, beta = BETA, scale_horizontal = SCALE_HORIZONTAL, scale_vertical = SCALE_VERTICAL, a1 = A1, b1 = B1, c1 = C1, d1 = D1, a2 = A2, b2 = B2, c2 = C2, d2 = D2, radius = earthRadiusKm)

#library(MASS)

start_time = Sys.time()

set.seed(1234)
#Z <- mvrnorm(1, mu = rep(0, ncol(cov_mat)), Sigma = cov_mat)

#save(Z, file = 'Z.RData')

end_time = Sys.time()

TOTAL_TIME <- as.numeric(end_time - start_time, units = "secs")

print(TOTAL_TIME)

args <- commandArgs(trailingOnly = TRUE)
NUM_PROCESSORS = as.numeric(args[1])

cov_mat <- cov_bi_differential_parallel(location = loc3d, beta = BETA, scale_horizontal = SCALE_HORIZONTAL, scale_vertical = SCALE_VERTICAL, a1 = A1, b1 = B1, c1 = C1, d1 = D1, a2 = A2, b2 = B2, c2 = C2, d2 = D2, radius = earthRadiusKm, num_processors = NUM_PROCESSORS)

print("DONE . . . ")

