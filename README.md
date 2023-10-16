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

---

## Sample R codes using the DiffOp package

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

set.seed(1235)
INIT_C1_COEF <- runif(length(INNER_KNOTS1) + SPLINES_DEGREE + 1, -0.1, 0.1)

set.seed(1236)
INIT_C2_COEF <- runif(length(INNER_KNOTS2) + SPLINES_DEGREE + 1, -0.1, 0.1)

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
                                          d1_fix = TRUE, d2_fix = TRUE, radius = earthRadiusKm,
                                          splines_degree = SPLINES_DEGREE,
                                          inner_knots1 = INNER_KNOTS1,
                                          inner_knots2 = INNER_KNOTS2,
                                          iterlim = 1000, stepmax = 1, hessian = F)
```

Save the R codes in a file named `testing_diffop_package.R`.

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
