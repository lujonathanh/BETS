In our manuscript, we analyzed possible associations in a separate over-expression data set. We used a VAR model that regressed each effect gene’s expression level on its previous expression level and the causal gene’s
previous expression level, assuming normal noise $\epsilon_t \sim \mathcal{N} (0, 1)$:

$$X^g_t = c^g X^{g}_{t-1} + d^{g', g}X^{g'}_{t-1} + \epsilon_t$$

No regularization was included, and ordinary least squares was used to fit the equation. Note that the expression $X^{g'}_
of a causal gene $g'$ is fit as a single predictor without the other expression. Lag $1$ is used due to the larger time gaps in the over-expression.

Here we descripe how to fit this pairwise model. We assume you've completed step 0 in `BETS_tutorial.md`

# 1. Preparing for a Run

1. Get to the BETS directory from the command line.
1. Set the parameters at `code/package_params_cpipeline.sh`
1. (OPTIONAL) If you want to run on your own computing cluster, modify `code/run_all_parallel_wait.sh`
1. Package for the cluster
  * `cd code/`
  * `source ./package_params_cpipeline.sh`
  * `./package_for_cluster_cpipeline.sh`
  * `cd $FOLDER`
  * `./prep_jobs_rand_cv.sh`
  * `./prep_jobs_pairwise.sh`
  
# 2. Fit the pairwise jobs.
1. `source ./package_params_cpipeline.sh`
1. `export scriptlist=pairwise_parallel_script_list.txt`
1. If on your own computer, do `./run_all_parallel_no_cluster.sh`, otherwise do `./run_all_parallel_wait.sh`
1. Wait for the jobs to complete.
1. Finish with `./pairwise-finish-scripts/finish.sh`

## 3. Format the output.
1. Put all the timing results together now that it's done. `./summarize_time.sh`
1. Organize the results. `./downstream_prep.sh`
1. All the results are now under `run_pairwise`
