# BETS Tutorial

In overview, BETS:

* takes as input (or multiple replicates) gene matrix and a randomized version of the gene matrix, where each gene's expression has been permuted over time.
* creates a "run" folder where the input and code are copied.
* runs the scripts in parallel
* outputs a network with FDR thresholding performed based off of bootstrap frequencies.


# Preparing for a Run

1. Set the parameters at `code/package_params_cpipeline.sh`
1. [Optional] If you want to run on your own computing cluster, modify `code/run_all_parallel_wait.sh`
1. Package for the cluster
  * `cd code/`
  * `source ./package_params_cpipeline.sh`
  * `./package_for_cluster_cpipeline.sh`
  * `cd $FOLDER`
  * `./prep_jobs_rand_cv.sh`
  * `./prep_jobs_bootstrap.sh`

# Running BETS (in detail)

## 1. Set hyperparameters
1. Set the list of scripts to run from. `export scriptlist=cv_parallel_script_list.txt`
1. If on your own computer, do `./run_all_parallel_no_cluster.sh`, otherwise do `./run_all_parallel_wait.sh`
1. Wait for the jobs to complete.
1. Set the hyperparameter for the fit. `./set_hyper.sh`

## 2. Fit the model on the original data.
1. `source ./package_params_cpipeline.sh`
1. `export scriptlist=fit_parallel_script_list.txt`
1. If on your own computer, do `./run_all_parallel_no_cluster.sh`, otherwise do `./run_all_parallel_wait.sh`
1. Wait for the jobs to complete.
1. `./finish-effect.sh`

## 3. Perform bootstrap stability selection.
1. `source ./package_params_cpipeline.sh`
1. `export scriptlist=bootstrap_parallel_script_list.txt`
1. If on your own computer, do `./run_all_parallel_no_cluster.sh`, otherwise do `./run_all_parallel_wait.sh`
1. Wait for the jobs to complete.
1. `export scriptlist=finish-effect-bootstrap-all.sh`
1. If on your own computer, do `./run_all_parallel_no_cluster.sh`, otherwise do `./run_all_parallel_wait.sh`
1. Wait for the jobs to complete.
1. `./get_result_bootstrap_lite.sh`
1. `./get_result_bootstrap-fdr-0.05-effect_lite.sh`

## 4. Format the output.
1. Put all the timing results together now that it's done. '`./summarize_time.sh`
1. Write all results to `run_l-fdr`. `./downstream_prep.sh`
