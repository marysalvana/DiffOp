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

1. Download the github repository:
     ```
     git clone https://github.com/marysalvana/DiffOp.git
     ```

2. Load the necessary modules:
   + In University of Houston Carya:
     ```
     module load R/4.2.0-foss-2021b
     ```

3. Build the R package:
     ```
     R CMD build DiffOp
     ```
This will create the `DiffOp_1.0.0.tar.gz` file.

4. Install the R package:
     ```
     R CMD INSTALL DiffOp_1.0.0.tar.gz
     ```
