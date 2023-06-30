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
   + In University of Houston (UH) Carya cluster:
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
    #SBATCH -n 4

    echo "Starting at `date`"
    echo "Running on hosts: $SLURM_NODELIST"
    echo "Running on $SLURM_NNODES nodes."
    echo "Running $SLURM_NTASKS tasks."
    echo "Current working directory is `pwd`"

    module load R/4.2.0-foss-2021b

    # set path variable
    path=`pwd -P`

    echo "FITTING MODEL ON REAL DATASET"

    mpiexec -np $SLURM_NTASKS Rscript ./testing_diffop_package.R  $SLURM_NTASKS
    ```

Save the job script in a file named `job.sub`. 

---

## Sample R codes using the DiffOp package



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

## Authors

DiffOp is authored and maintained by:
* [Mary Lai Salvana](https://marylaisalvana.com)
* [Mikyoung Jun](https://sites.google.com/view/mikyoung-jun/home?authuser=0)
