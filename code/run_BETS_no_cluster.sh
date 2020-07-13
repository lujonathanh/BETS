#!/bin/bash
set -e

### A. Package folder, prep the scripts that will be run
cd code/
source ./package_params_cpipeline.sh
./package_for_cluster_cpipeline.sh
cd $FOLDER
./prep_jobs_rand_cv.sh
./prep_jobs_bootstrap.sh

### B. Set hyperparameters
# Set the list of scripts to run from.
export scriptlist=cv_parallel_script_list.txt
#1. If on your own computer, do
./run_all_parallel_no_cluster.sh
#   If submitting jobs to cluster, do `./run_all_parallel_wait.sh`
# Wait for the jobs to complete.
# Set the hyperparameter for the fit.
./set_hyper.sh

### C. Fit the model on the original data.
source ./package_params_cpipeline.sh
export scriptlist=fit_parallel_script_list.txt
#1. If on your own computer, do
./run_all_parallel_no_cluster.sh
#   If submitting jobs to cluster, do `./run_all_parallel_wait.sh`
# Wait for the jobs to complete.
./finish-effect.sh

### D. Perform stability selection (from bootstrap samples).
source ./package_params_cpipeline.sh
export scriptlist=bootstrap_parallel_script_list.txt
#1. If on your own computer, do
./run_all_parallel_no_cluster.sh
#   If submitting jobs to cluster, do `./run_all_parallel_wait.sh`
# Wait for the jobs to complete.
export scriptlist=finish-effect-bootstrap_parallel_script_list.txt
#1. If on your own computer, do
./run_all_parallel_no_cluster.sh
#   If submitting jobs to cluster, do `./run_all_parallel_wait.sh`
# Wait for the jobs to complete.
# Combine the bootstrap elastic net fits.
./get_result_bootstrap_lite.sh
# Combine the significant networks for each bootstrap sample.
./get_result_bootstrap-fdr-0.05-effect_lite.sh

### E. Format the output.
# Put all the timing results together now that it's done.
./summarize_time.sh
# Organize the results.
./downstream_prep.sh
# All the results are now under `run_l-fdr`