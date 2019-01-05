############
FOLDERS
############

plots/:
-- betas/: Causal coefficients from original & randomized data
---- beta_nonzero_lag-k.png: overlaid betas from normal and randomized, for lag k
---- beta_abs_lag-k.png: overlaid absolute-betas from normal and randomized, for lag k

-- fdr-*/: expression time series for significant gene pairs
---- ____-ENSG_lag-k-preds-l.png: plot of the lth set of significant predictors for gene ENSG. Note they may be lag-1 or lag-2 causal.

-- hyper/: hyperparameters VS held-out r^2, mse, d.o.f. of coefficients

fdr-*/:
___acoefs-lag-k-fdr-*.txt: cause-by-effect matrix of significant coefs for the kth lag. Note if total lag was L, 1 <= k <= L.
__acoefs-lag-k-fdr-*-network.txt: Above in network format, each row is one edge
__threshes.p: pickle file with the threshold used to achieve FDR
__acoefs-lag-k-fdr-*-README.txt: Get README of the network statistics.


hyper/:
results from hyperparameter fits


############
FILES
############

fit_all_summary_normal.txt: summarize normal fit network

fit_all_summary_random.txt: summarize random fit network

fit_all_summary_fdr.txt: summarize post-FDR network

__acoefs-lag-k.txt: cause-by-effect matrix of coefs from original data for the kth lag. Note if total lag was L, 1 <= k <= L.
__acoefsr-lag-k.txt: above, but for randomized data

__fit_result_df.txt: fit results for original data for the kth lag. Note if total lag was L, 1 <= k <= L.
__fit_result_dfr.txt: above, but for randomized data

__coefs.p: Pickle file of the fitted coefficients. NOTE: NOT ALIGNED,  NOT FOR DOWNSTREAM USE.
__coefsr.p: above for randomized data. NOTE: NOT ALIGNED,  NOT FOR DOWNSTREAM USE.

__intercepts.p: Pickle file of the fitted


####################
SUBMISSION PROCEDURE
####################

Procedure
* ./package_for_cluster_cpipeline.sh
* ./prep_jobs_rand_cv.sh
* source ./package_params_cpipeline.sh
* export scriptlist=cv_parallel_script_list.txt
* ./run_all_parallel_wait.sh
* (wait)
* ./set_hyper.sh
* source ./package_params_cpipeline.sh
* export scriptlist=fit_parallel_script_list.txt
* ./run_all_parallel_wait.sh
* (wait)
* ./finish.sh