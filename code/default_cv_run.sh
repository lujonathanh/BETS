source ./package_params_cpipeline.sh
./prep_jobs_rand_cv.sh
export scriptlist=cv_parallel_script_list.txt
./run_all_parallel_wait.sh