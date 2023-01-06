# BETS Tutorial

In overview, BETS:

* takes as input (or multiple replicates of) gene matrix and a randomized version of the gene matrix, where each gene's expression has been permuted over time.
* creates a "run" folder where the input and code are copied.
* runs the scripts in parallel
* outputs a network with FDR thresholding performed based off of bootstrap frequencies.

# Input

BETS takes in a tab-delimited gene by timepoint file. It should have the first column be the genenames, and the columns after be the ordered timepoints (left is earlier, right is later). See `data/DREAM/insilico/0mean/insilico_size100_1_0mean_TS-rep-1.txt` as an example.

`gene   time1   time2 time3 ...  `

`geneA   0.1 0.1 -0.2    ...  `

`geneB   0.3 -0.2    -0.1    ...  `

`...`

If there are multiple replicates, you instead provide a text file with list of each of the individual replicate files. See `data/DREAM/insilico_size100_1/0mean/reps.txt` as an example: it lists individual replicate files like `data/DREAM/insilico/0mean/insilico_size100_1_0mean_TS-rep-1.txt` and `data/DREAM/insilico/0mean/insilico_size100_1_0mean_TS-rep-2.txt`. Each individual replicate files have the same format as above.

BETS treats replicates as independent samples, so please make sure *your replicates are measured at the same timepoints and list genes in the same order*.

# Running BETS (in detail)

## 0. Make sure you have the libraries installed in Python 3

`python3 -m pip install -r requirements.txt`

## 1. Modify key run parameters

1. Get to the BETS directory from the command line.
1. Set the parameters at `code/package_params_cpipeline.sh`
1. (OPTIONAL) If you want to run on a computing cluster, modify `code/run_all_parallel_wait.sh` so that it submits jobs appropriately.

## 2. Package the run folder, cd into it

  * `cd code/`
  * `source ./package_params_cpipeline.sh`
  * `./package_for_cluster_cpipeline.sh`
  * `cd $FOLDER`

## 3. Perform the Run.

Assuming you are in the run folder, a shortcut to run all of the following parts of step 3A-E automatically on your computer without having to intervene step-by-step is:

`./run_BETS_no_cluster.sh`


### 3A. Prep the scripts that will be run
  * `source ./package_params_cpipeline.sh`
  * `./prep_jobs_rand_cv.sh`
  * `./prep_jobs_bootstrap.sh`

### 3B. Set hyperparameters
1. Set the list of scripts to run from. `export scriptlist=cv_parallel_script_list.txt`
1. If on your own computer, do `./run_all_parallel_no_cluster.sh`  
   If submitting jobs to cluster, do `./run_all_parallel_wait.sh`
1. Wait for the jobs to complete. 
   (you can ignore the `line 3: module: command not found` errors)
1. Set the hyperparameter for the fit. `./set_hyper.sh`

### 3C. Fit the model on the original data.
1. `source ./package_params_cpipeline.sh`
1. `export scriptlist=fit_parallel_script_list.txt`
1. If on your own computer, do `./run_all_parallel_no_cluster.sh`. If submitting jobs to cluster, do`./run_all_parallel_wait.sh`
1. Wait for the jobs to complete.
   (you can ignore the `line 3: module: command not found` errors)
1. `./finish-effect.sh`

### 3D. Perform stability selection (from bootstrap samples).
1. `source ./package_params_cpipeline.sh`
1. `export scriptlist=bootstrap_parallel_script_list.txt`
1. If on your own computer, do `./run_all_parallel_no_cluster.sh`. If submitting jobs to cluster, do `./run_all_parallel_wait.sh`
1. Wait for the jobs to complete. 
   (you can ignore the `line 3: module: command not found` errors)
1. `export scriptlist=finish-effect-bootstrap_parallel_script_list.txt`
1. If on your own computer, do `./run_all_parallel_no_cluster.sh`. If submitting jobs to cluster, do `./run_all_parallel_wait.sh`
1. Wait for the jobs to complete.
   (you can ignore the `line 3: module: command not found` errors)
1. Combine the bootstrap elastic net fits. `./get_result_bootstrap_lite.sh`
1. Combine the significant networks for each bootstrap sample. `./get_result_bootstrap-fdr-0.05-effect_lite.sh`

### 3E. Format the output.
1. Put all the timing results together now that it's done. `./summarize_time.sh`
1. Organize the results. `./downstream_prep.sh`
1. All the results are now under `run_l-fdr`

## 4. Run with permuted data.
1. Edit `package_params_cpipeline.sh`:

  A. replace your `DATAFILE` with a version of `DATAFILE` where every gene's temporal profile has been independently shuffled across time, separately for distinct replicates. In this example, change

`export DATAFILE=../data/DREAM/insilico_size100_1/0mean/reps.txt`

to

`export DATAFILE=../data/DREAM/insilico_size100_1/0mean/reps-urand.txt`

(note this permuted data set is distinct from `RANDDATAFILE=../data/DREAM/insilico_size100_1/0mean/reps-rand.txt`.  Both are generated in the same way, but with different random seeds.)

  B. add as a suffix of `_urand` to GENES. In this example, change
  
`export GENES=insilico_size100_1`

to 

`export GENES=insilico_size100_1_urand`

2. Run sections 1-3, above, as before.

## 5. Perform false discovery thresholding.
1. Change into the run directory with the original data as in Step 3.

In this example, change into `BETS/runs/insilico_size100_1-0mean-reps-enet-2-g`

2. Modify `get_umbrella_null_results.py`

`urand_run_folder = "../URANDFOLDERNAME/run_l-fdr"`

In this example, change to `urand_run_folder = "../insilico_size100_1_urand-0mean-reps-enet-2-g/run_l-fdr"`

3. Do `python3 get_umbrella_null_results.py`

4. The FDR-thresholded network is now available at `analysis/bootstrap-fdr-0.2-network.txt` and the bootstrap frequency threshold is available at `analysis/bootstrap-fdr-0.2-network-umbrella-results.csv`

In this example, you should get a network with 3364 edges.

# Questions?

Reach out at the [Google Group](https://groups.google.com/forum/#!forum/bets-support)!
