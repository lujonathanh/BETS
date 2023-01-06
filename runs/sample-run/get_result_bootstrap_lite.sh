START=$(date)
time python3 get_result_bootstrap.py -df reps.txt -lr 1 -osf bootstrap -rsf bootstrap/bootstrap-results -o insilico_size100_1-0mean-reps-enet-2-g_lite -l 2 -tn e -b all_bootstrap_coefs.txt -da 1 -dl 1 -uabrd 0
time python3 get_intercept_bootstrap.py -b all_bootstrap_intercepts.txt -rsf bootstrap/bootstrap-results -o insilico_size100_1-0mean-reps-enet-2-g
END=$(date)
echo get_result_bootstrap_lite.sh,$START,$END,$SECONDS >> timing/bootstrap_result_time.csv
