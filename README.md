## DiffOp

An R package that implements a flexible multivariate
  nonstationary cross-covariance function model based on the differential
  operators approach defined in 3-dimensional (3D) space. A defining
  feature of the model is that it can accommodate nonstationarity in the
  variances, colocated correlations, and other spatial features of a
  multivariate Gaussian random field, and is particularly flexible along
  the vertical dimension. With the implemented cross-covariance function
  model, users can simulate synthetic multivariate Gaussian random fields.
  Inference for model parameters is done via maximum likelihood estimation.
  Predictions at locations with no observations can also be performed.

---

### License

This package is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License, version 3, as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the GNU
General Public License for more details.

A copy of the GNU General Public License, version 3, is available at
<https://www.r-project.org/Licenses/GPL-3>

---

## Installation

To install the package:

* Clone the DiffOp package from [here](https://github.com/marysalvana/DiffOp) by running the following commands on the terminal: `git clone https://github.com/marysalvana/DiffOp.git`.
* Build the package by running the following commands on the terminal: `R CMD build DiffOp`. A `.tar.gz` file will appear with filename `DiffOp_1.0.0.tar.gz`.
* Install the package by running the following commands on the terminal: `R CMD INSTALL DiffOp_1.0.0.tar.gz` or `R CMD INSTALL DiffOp_1.0.0.tar.gz --library=/home/salvanmo/R/x86_64-pc-linux-gnu-library/3.6/` when need to install in specific directory. The package is ready to be used in `R` and will be loaded when you run the code `library(DiffOp)` in `R`. If the installation directory is different from the default, change `.libPaths()` in `R` by running the command: `.libPaths("/home/salvanmo/R/x86_64-pc-linux-gnu-library/3.6/")`.

---

## Sample R codes using the DiffOp package
### Synthetic data generation and model fitting

```
library(DiffOp)

library(dplyr)

x <- seq(0, 1, length.out = 10)
y <- seq(0, 1, length.out = 10)
loc2d <- expand.grid(x, y) %>% as.matrix()

depth <- seq(0, 1, length.out = 10)
loc3d <- cbind(rep(loc2d[, 1], each = length(depth)),
               rep(loc2d[, 2], each = length(depth)), depth)

earthRadiusKm = 6371
BETA = 0.5
SCALE_HORIZONTAL = 0.03
SCALE_VERTICAL = 0.3
A1 = A2 = 0.00001
B1 = B2 = 0.00001
C1 <- sin((loc3d[, 3] + 0.1) * pi / 0.5)
C2 <- cos((loc3d[, 3] + 0.1) * pi / 0.5)
D1 = D2 = 0

cov_mat <- cov_bi_differential(location = loc3d, beta = BETA,
                               scale_horizontal = SCALE_HORIZONTAL, scale_vertical = SCALE_VERTICAL,
                               a1 = A1, b1 = B1, c1 = C1, d1 = D1, a2 = A2, b2 = B2, c2 = C2, d2 = D2,
                               radius = earthRadiusKm)

library(MASS)

set.seed(1234)
Z <- mvrnorm(1, mu = rep(0, ncol(cov_mat)), Sigma = cov_mat)
Z1 <- Z[1:nrow(loc3d)]
Z2 <- Z[nrow(loc3d) + 1:nrow(loc3d)]

INIT_BETA = 0
INIT_SCALE_HORIZONTAL = log(0.02)
INIT_SCALE_VERTICAL = log(0.2)
INIT_A1 = INIT_A2 = 0
INIT_B1 = INIT_B2 = 0
INIT_D1 = INIT_D2 = 0

INNER_KNOTS1 <- c(0.1, 0.5, 0.9)
INNER_KNOTS2 <- c(0.1, 0.5, 0.9)

SPLINES_DEGREE = 2

no_of_c1_coef = length(INNER_KNOTS1) + SPLINES_DEGREE + 1
no_of_c2_coef = length(INNER_KNOTS2) + SPLINES_DEGREE + 1

set.seed(1235)
INIT_C1_COEF <- runif(no_of_c1_coef, -0.1, 0.1)

set.seed(1236)
INIT_C2_COEF <- runif(no_of_c2_coef, -0.1, 0.1)

#Fitting a stationary model

est_params_mle <- est_bi_differential_mle(residuals = Z,
                                          location = loc3d, init_beta = 0,
                                          init_scale_horizontal = log(0.1),
                                          init_scale_vertical = log(0.1),
                                          init_a1 = INIT_A1, init_b1 = INIT_B1,
                                          init_c1_coef = rep(0, 6), init_d1 = 1,
                                          init_a2 = INIT_A2, init_b2 = INIT_B2,
                                          init_c2_coef = rep(0, 6), init_d2 = 1,
                                          a1_scaling = 1e-3, b1_scaling = 1e-3,
                                          a2_scaling = 1e-3, b2_scaling = 1e-3,
                                          beta_fix = T, #scale_horizontal_fix = T, scale_vertical_fix = T,
                                          a1_fix = T, b1_fix = T, a2_fix = T, b2_fix = T,
                                          c1_fix = T, c2_fix = T,
                                          radius = earthRadiusKm,
                                          splines_degree = SPLINES_DEGREE,
                                          inner_knots1 = INNER_KNOTS1,
                                          inner_knots2 = INNER_KNOTS2,
                                          iterlim = 1000, stepmax = 1, hessian = F)

#Fitting a nonstationary model, such that we fix the scale_horizontal and scale_vertical
#  parameters to the estimates above. 

est_params_mle <- est_bi_differential_mle(residuals = Z,
                                          location = loc3d, init_beta = 0,
                                          init_scale_horizontal = log(0.004154042),
                                          init_scale_vertical = log(0.3120303),
                                          init_a1 = INIT_A1, init_b1 = INIT_B1,
                                          init_c1_coef = rep(0, 6), init_d1 = INIT_D1,
                                          init_a2 = INIT_A2, init_b2 = INIT_B2,
                                          init_c2_coef = rep(0, 6), init_d2 = INIT_D2,
                                          a1_scaling = 1e-3, b1_scaling = 1e-3,
                                          a2_scaling = 1e-3, b2_scaling = 1e-3,
                                          scale_horizontal_fix = T, scale_vertical_fix = T,
                                          a1_fix = T, b1_fix = T, a2_fix = T, b2_fix = T,
                                          d1_fix = T, d2_fix = T, radius = earthRadiusKm,
                                          splines_degree = SPLINES_DEGREE,
                                          inner_knots1 = INNER_KNOTS1,
                                          inner_knots2 = INNER_KNOTS2,
                                          iterlim = 1000, stepmax = 1, hessian = F)
```

### Argo data analysis and model fitting
```
data("argo_ref_loc1")

ind_pred <- 1:50

loc3d <- cbind(argo_ref_loc1$Longitude, argo_ref_loc1$Latitude,
               argo_ref_loc1$Pressure)
locs_insample <- loc3d[-ind_pred, ]
locs_outsample <- loc3d[ind_pred, ]

Z_insample <- c(argo_ref_loc1$TemperatureResiduals[-ind_pred],
                argo_ref_loc1$SalinityResiduals[-ind_pred])
Z_outsample <- c(argo_ref_loc1$TemperatureResiduals[ind_pred],
                 argo_ref_loc1$SalinityResiduals[ind_pred])

earthRadiusKm = 6371

INIT_BETA = 0.9
INIT_SCALE_HORIZONTAL = log(0.1)
INIT_SCALE_VERTICAL = log(0.1)
INIT_A1 = INIT_A2 = 0
INIT_B1 = INIT_B2 = 0
INIT_D1 = INIT_D2 = 0

est_params_mle_step1 <- est_bi_differential_mle(residuals = Z_insample,
                                          location = locs_insample, init_beta = 1,
                                          init_scale_horizontal = INIT_SCALE_HORIZONTAL,
                                          init_scale_vertical = INIT_SCALE_VERTICAL,
                                          init_a1 = INIT_A1, init_b1 = INIT_B1,
                                          init_c1_coef = 1, init_d1 = 0,
                                          init_a2 = INIT_A2, init_b2 = INIT_B2,
                                          init_c2_coef = 1, init_d2 = 0,
                                          a1_scaling = 1e-3, b1_scaling = 1e-3,
                                          a2_scaling = 1e-3, b2_scaling = 1e-3,
                                          beta_fix = T, scale_horizontal_fix = F, scale_vertical_fix = F,
                                          a1_fix = F, b1_fix = F, a2_fix = F, b2_fix = F,
                                          c1_fix = F, c2_fix = F,
                                          d1_fix = T, d2_fix = T,
                                          radius = earthRadiusKm,
                                          splines_degree = 0,
                                          iterlim = 1000, stepmax = 1, hessian = F)

INNER_KNOTS1 <- seq(100, 1000, by = 100)
INNER_KNOTS2 <- seq(100, 1000, by = 100)

SPLINES_DEGREE = 1

no_of_c1_coef = length(INNER_KNOTS1) + SPLINES_DEGREE + 1
no_of_c2_coef = length(INNER_KNOTS2) + SPLINES_DEGREE + 1

set.seed(1234)
INIT_C1_COEF <- runif(no_of_c1_coef, -4, 4)
set.seed(1235)
INIT_C2_COEF <- runif(no_of_c2_coef, -4, 4)

est_params_mle_step2 <- est_bi_differential_mle(residuals = Z_insample,
                                                location = locs_insample, init_beta = 0,
                                                init_scale_horizontal = exp(-3.11236305),
                                                init_scale_vertical = -3.75414498,
                                                init_a1 = 0, init_b1 = 0,
                                                init_c1_coef = INIT_C1_COEF, init_d1 = 0,
                                                init_a2 = 0, init_b2 = 0,
                                                init_c2_coef = INIT_C2_COEF, init_d2 = 0,
                                                beta_fix = T,
                                                scale_horizontal_fix = T, scale_vertical_fix = F,
                                                a1_fix = T, b1_fix = T, a2_fix = T, b2_fix = T,
                                                d1_fix = T, d2_fix = T, radius = earthRadiusKm,
                                                splines_degree = SPLINES_DEGREE,
                                                inner_knots1 = INNER_KNOTS1,
                                                inner_knots2 = INNER_KNOTS2,
                                                iterlim = 1000, stepmax = 1, hessian = F)

#When optimization stops because of the prompt "Maximum step size exceeded 5 consecutive times.", just run the MLE again from the last parameter values as starting point. 

theta = c(0.23007337,-6.32609215,2.26555344,0.53293427,0.65316327,0.47183072,0.74479939,-1.58372243,-0.88918217,0.09310385,0.1114813,0.16532703,-0.07239835,-0.11218514,-0.14882934,-0.13306325,0.128212,0.04194192,-0.0648904,-0.02101291,-0.01829493,0.07032543,-0.03063398,-0.01355103,0.00168753)

for(ll in 1:100){
  theta <- est_params_mle_step2$theta
  est_params_mle_step2 <- est_bi_differential_mle(residuals = Z_insample,
                                                  location = locs_insample, init_beta = 0,
                                                  init_scale_horizontal = exp(-3.11236305),
                                                  init_scale_vertical = theta[1],
                                                  init_a1 = 0, init_b1 = 0,
                                                  init_c1_coef = theta[1 + 1:no_of_c1_coef], init_d1 = 0,
                                                  init_a2 = 0, init_b2 = 0,
                                                  init_c2_coef = theta[1 + no_of_c1_coef + 1:no_of_c2_coef], init_d2 = 0,
                                                  beta_fix = T,
                                                  scale_horizontal_fix = T, scale_vertical_fix = F,
                                                  a1_fix = T, b1_fix = T, a2_fix = T, b2_fix = T,
                                                  d1_fix = T, d2_fix = T, radius = earthRadiusKm,
                                                  splines_degree = SPLINES_DEGREE,
                                                  inner_knots1 = INNER_KNOTS1,
                                                  inner_knots2 = INNER_KNOTS2,
                                                  iterlim = 1000, stepmax = 1, hessian = F)
}

#Plotting the c1 and c2 functions using MLE values when BETA is not estimated

basis1 <- bsplineBasis(locs_insample[, 3], SPLINES_DEGREE, INNER_KNOTS1)
nb1 <- ncol(basis1)
basis2 <- bsplineBasis(locs_insample[, 3], SPLINES_DEGREE, INNER_KNOTS1)
nb2 <- ncol(basis2)

c1_coef = theta[1 + 1:no_of_c1_coef]
c2_coef = theta[1 + no_of_c1_coef + 1:no_of_c2_coef]

c1 <- basis1 %*% matrix(c1_coef, ncol = 1)
c2 <- basis2 %*% matrix(c2_coef, ncol = 1)
plot(c1[1:50])
plot(c2[1:50])

#Plotting the marginal variances and colocated correlations using MLE values when BETA is not estimated
cov_mat <- cov_bi_differential(location = locs_insample, beta = 1,
                               scale_horizontal = exp(-3.11236305), 
                               scale_vertical = exp(theta[1]),
                               a1 = 0, b1 = 0,
                               c1_coef = c1_coef, d1 = 0,
                               a2 = 0, b2 = 0,
                               c2_coef = c2_coef, d2 = 0,
                               radius = earthRadiusKm, splines_degree = SPLINES_DEGREE,
                               inner_knots1 = INNER_KNOTS1, inner_knots2 = INNER_KNOTS2)

plot(diag(cov_mat[1:50, 1:50]))
plot(diag(cov_mat[2450 + 1:50, 2450 + 1:50]))
plot(diag(cov_mat[1:50, 2450 + 1:50]) / sqrt(diag(cov_mat[1:50, 1:50]) * diag(cov_mat[2450 + 1:50, 2450 + 1:50])))

est_params_mle_step3 <- est_bi_differential_mle(residuals = Z_insample,
                                                location = locs_insample, init_beta = 1,
                                                init_scale_horizontal = exp(-3.11236305),
                                                init_scale_vertical = theta[1],
                                                init_a1 = 0, init_b1 = 0,
                                                init_c1_coef = theta[1 + 1:no_of_c1_coef], init_d1 = 0,
                                                init_a2 = 0, init_b2 = 0,
                                                init_c2_coef = theta[1 + no_of_c1_coef + 1:no_of_c2_coef], init_d2 = 0,
                                                beta_fix = T,
                                                scale_horizontal_fix = T, scale_vertical_fix = F,
                                                a1_fix = T, b1_fix = T, a2_fix = T, b2_fix = T,
                                                d1_fix = T, d2_fix = T, radius = earthRadiusKm,
                                                splines_degree = SPLINES_DEGREE,
                                                inner_knots1 = INNER_KNOTS1,
                                                inner_knots2 = INNER_KNOTS2,
                                                iterlim = 1000, stepmax = 1, hessian = F)

theta = c(0.21772112,-6.31371709,2.26079917,0.65465982,0.69462169,0.48813674,0.74496263,-1.58737371,-0.90051984,0.09055384,0.13868718,0.1675015,-0.07370263,-0.80400007,0.28395557,0.03499888,0.05125029,0.03182939,0.02715794,-0.04775821,0.12488736,-0.013939,-0.02279462,-0.02816342,0.01209895)

for(ll in 1:100){
  theta <- est_params_mle_step3$theta
  est_params_mle_step3 <- est_bi_differential_mle(residuals = Z_insample,
                                                  location = locs_insample, init_beta = 1,
                                                  init_scale_horizontal = exp(-3.11236305),
                                                  init_scale_vertical = theta[1],
                                                  init_a1 = 0, init_b1 = 0,
                                                  init_c1_coef = theta[1 + 1:no_of_c1_coef], init_d1 = 0,
                                                  init_a2 = 0, init_b2 = 0,
                                                  init_c2_coef = theta[1 + no_of_c1_coef + 1:no_of_c2_coef], init_d2 = 0,
                                                  beta_fix = T,
                                                  scale_horizontal_fix = T, scale_vertical_fix = F,
                                                  a1_fix = T, b1_fix = T, a2_fix = T, b2_fix = T,
                                                  d1_fix = T, d2_fix = T, radius = earthRadiusKm,
                                                  splines_degree = SPLINES_DEGREE,
                                                  inner_knots1 = INNER_KNOTS1,
                                                  inner_knots2 = INNER_KNOTS2,
                                                  iterlim = 1000, stepmax = 1, hessian = F)
  
}

c1_coef = theta[1 + 1:no_of_c1_coef]
c2_coef = theta[1 + no_of_c1_coef + 1:no_of_c2_coef]

c1 <- basis1 %*% matrix(c1_coef, ncol = 1)
c2 <- basis2 %*% matrix(c2_coef, ncol = 1)
plot(c1[1:50])
plot(c2[1:50])

#Plotting the marginal variances and colocated correlations using MLE values when BETA is not estimated
cov_mat <- cov_bi_differential(location = locs_insample, beta = 1,
                               scale_horizontal = exp(-3.11236305), 
                               scale_vertical = exp(theta[1]),
                               a1 = 0, b1 = 0,
                               c1_coef = c1_coef, d1 = 0,
                               a2 = 0, b2 = 0,
                               c2_coef = c2_coef, d2 = 0,
                               radius = earthRadiusKm, splines_degree = SPLINES_DEGREE,
                               inner_knots1 = INNER_KNOTS1, inner_knots2 = INNER_KNOTS2)

plot(diag(cov_mat[1:50, 1:50]))
plot(diag(cov_mat[2450 + 1:50, 2450 + 1:50]))
plot(diag(cov_mat[1:50, 2450 + 1:50]) / sqrt(diag(cov_mat[1:50, 1:50]) * diag(cov_mat[2450 + 1:50, 2450 + 1:50])))

```

---

## How to install the DiffOp package
1. Install the necessary software:
    - `R`:
        * For Apple silicon (M1/M2) Macs, make sure you installed the `R` binary for Apple silicon arm64. You can download [R-4.3.1-arm64.pkg](https://cran.r-project.org/bin/macosx/) from CRAN.
        * In University of Houston (UH) Carya cluster, simply load the `R` module by running the following command on the terminal:
        ```
        module load R/4.2.0-foss-2021b
        ```
    - `OpenMPI`:
        * For Apple silicon (M1/M2) Macs, `OpenMPI` can be installed using `brew` by running the following command on the terminal: 
        ```
        brew install open-mpi
        ```
        * In UH Carya cluster, loading the `R` module above automatically loads `OpenMPI`. No need to run any additional commands.
    - `gfortran`:
        * For Apple silicon (M1/M2) Macs, `gfortran` can be installed using `brew` by running the following command on the terminal: 
        ```
        brew install gcc
        ```
        * In UH Carya cluster, loading the `R` module above automatically loads `gcc`. No need to run any additional commands.
    - `GSL`:
        * For Apple silicon (M1/M2) Macs, `GSL` can be installed using `brew` by running the following command on the terminal: 
        ```
        brew install gsl
        ```
        * In UH Carya cluster, loading the `R` module above automatically loads `GSL`. No need to run any additional commands.

2. Create a Makevars file and define the paths to the source files of the software in Step #1:
     ```
     mkdir ~/.R
     vi ~/.R/Makevars
     ```
     Inside `~/.R/Makevars`, add the following:
     ```
     export CPATH=/opt/homebrew/Cellar/gsl/2.7.1/include/
     export LIBRARY_PATH=/opt/homebrew/Cellar/gsl/2.7.1/lib/
     export LD_LIBRARY_PATH=/opt/homebrew/Cellar/gsl/2.7.1/lib/:$LD_LIBRARY_PATH

     FC = /usr/local/gfortran/bin/gfortran
     F77 = /usr/local/gfortran/bin/gfortran
     FLIBS = -L/usr/local/gfortran/lib
     ```

3. Download the necessary `R` packages from Github:
     ```
     git clone https://github.com/RBigData/pbdMPI.git
     git clone https://github.com/RBigData/pbdSLAP.git
     git clone https://github.com/marysalvana/pbdBASE.git
     git clone https://github.com/RBigData/pbdDMAT.git
     ```

4. Install the `R` packages in Step #2 in the following order:
     ```
     R CMD build pbdMPI	#This will create the `pbdMPI_0.4-7.tar.gz` file which you will INSTALL next.
     R CMD INSTALL pbdMPI_0.4-7.tar.gz

     R CMD build pbdSLAP	#This will create the `pbdSLAP_0.3-2.tar.gz` file which you will INSTALL next.
     R CMD INSTALL pbdSLAP_0.3-2.tar.gz

     R CMD build pbdBASE	#This will create the `pbdBASE_0.5-3.tar.gz` file which you will INSTALL next.
     R CMD INSTALL pbdBASE_0.5-3.tar.gz

     R CMD build pbdDMAT	#This will create the `pbdDMAT_0.5-2.tar.gz` file which you will INSTALL next.
     R CMD INSTALL pbdDMAT_0.5-2.tar.gz
     ```

5. Download the github repository:
     ```
     git clone https://github.com/marysalvana/DiffOp.git
     ```

6. Build the R package:
     ```
     R CMD build DiffOp
     ```
   This will create the `DiffOp_1.0.0.tar.gz` file.

7. Install the R package:
     * For general laptops, run the command:
     ```
     R CMD INSTALL DiffOp_1.0.0.tar.gz
     ```
     * In UH Carya cluster, run the command:
     ```
     R CMD INSTALL DiffOp_1.0.0.tar.gz --library=/project/jun/msalvana/R/x86_64-pc-linux-gnu-library/4.2/
     ```
---

## Sample job script to run in UH Carya the R codes using the DiffOp package

    ```
    #!/bin/bash
    #SBATCH -p batch
    #SBATCH -J DiffOp

    #SBATCH --time=10:00:00
    #SBATCH -o %j.out
    #SBATCH -e %j.err
    #SBATCH --mem=100GB

    #SBATCH --cpus-per-task=20
    #SBATCH -n 1

    echo "Starting at `date`"
    echo "Running on hosts: $SLURM_NODELIST"
    echo "Running on $SLURM_NNODES nodes."
    echo "Running $SLURM_NTASKS tasks."
    echo "Current working directory is `pwd`"

    module load R/4.2.0-foss-2021b

    # set path variable
    path=`pwd -P`

    echo "FITTING MODEL ON REAL DATASET"

    mpiexec -np $SLURM_NTASKS Rscript ./testing_diffop_package.R 
    ```

Save the job script in a file named `job.sub`. 


## How to run the R code

1. Submit the job script:
     ```
     sbatch job.sub
     ```
This command will generate a job ID number that you can use to monitor the job.

2. View the R codes output:
     ```
     tail -f job_id_number.out
     ```
   or 
     ```
     vi job_id_number.out
     ```

3. To check error logs from R code runs:
     ```
     vi job_id_number.err
     ```

## Authors

DiffOp is authored and maintained by:
* [Mary Lai Salvana](https://marylaisalvana.com)
* [Mikyoung Jun](https://sites.google.com/view/mikyoung-jun/home?authuser=0)
