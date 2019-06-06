# BETS Tutorial

In overview, BETS:

* takes as input (or multiple replicates of) gene matrix and a randomized version of the gene matrix, where each gene's expression has been permuted over time.
* creates a "run" folder where the input and code are copied.
* runs the scripts in parallel
* outputs a network with FDR thresholding performed based off of bootstrap frequencies.

# Input

BETS takes in a tab-delimited gene by timepoint file. It should have the first column be the genenames, and the columns after be the ordered timepoints (left is earlier, right is later). See `data/DREAM/insilico/0mean/insilico_size100_1_0mean_TS-rep-1.txt` as an example.

`gene   time1   time2 time3 ...
geneA   0.1 0.1 -0.2    ...
geneB   0.3 -0.2    -0.1    ...
...`

If there are multiple replicates, you instead give a 


# Running BETS (in detail)

## 1. Preparing for a Run

1. Get to the BETS directory from the command line.
1. Set the parameters at `code/package_params_cpipeline.sh`
1. (OPTIONAL) If you want to run on your own computing cluster, modify `code/run_all_parallel_wait.sh`
1. Package for the cluster
  * `cd code/`
  * `source ./package_params_cpipeline.sh`
  * `./package_for_cluster_cpipeline.sh`
  * `cd $FOLDER`
  * `./prep_jobs_rand_cv.sh`
  * `./prep_jobs_bootstrap.sh`

## 2. Set hyperparameters
1. Set the list of scripts to run from. `export scriptlist=cv_parallel_script_list.txt`
1. If on your own computer, do `./run_all_parallel_no_cluster.sh`  
   Otherwise do `./run_all_parallel_wait.sh`
1. Wait for the jobs to complete.
1. Set the hyperparameter for the fit. `./set_hyper.sh`

## 3. Fit the model on the original data.
1. `source ./package_params_cpipeline.sh`
1. `export scriptlist=fit_parallel_script_list.txt`
1. If on your own computer, do `./run_all_parallel_no_cluster.sh`, otherwise do `./run_all_parallel_wait.sh`
1. Wait for the jobs to complete.
1. `./finish-effect.sh`

## 4. Perform stability selection (from bootstrap samples).
1. `source ./package_params_cpipeline.sh`
1. `export scriptlist=bootstrap_parallel_script_list.txt`
1. If on your own computer, do `./run_all_parallel_no_cluster.sh`, otherwise do `./run_all_parallel_wait.sh`
1. Wait for the jobs to complete.
1. `export scriptlist=finish-effect-bootstrap-all.sh`
1. If on your own computer, do `./run_all_parallel_no_cluster.sh`, otherwise do `./run_all_parallel_wait.sh`
1. Wait for the jobs to complete.
1. Combine the bootstrap elastic net fits. `./get_result_bootstrap_lite.sh`
1. Combine the significant networks for each bootstrap sample. `./get_result_bootstrap-fdr-0.05-effect_lite.sh`

## 5. Format the output.
1. Put all the timing results together now that it's done. '`./summarize_time.sh`
1. Organize the results. `./downstream_prep.sh`
1. All the results are now under `run_l-fdr`

# Questions?

Reach out at the [Google Group](https://groups.google.com/forum/#!forum/bets-support)!
