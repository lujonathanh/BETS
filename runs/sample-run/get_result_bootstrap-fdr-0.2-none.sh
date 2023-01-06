START=$(date)
set -e
time python3 get_result_bootstrap.py -df reps.txt -lr 1 -osf bootstrap -rsf bootstrap/bootstrap-results-fdr-0.2-none -o insilico_size100_1-0mean-reps-enet-2-g-fdr-0.2-none -l 2 -tn e -b all_bootstrap_coefs_fdr-0.2-none.txt -da 0 -tbf bootstrap-transpose-fdr-0.2-none -uabrd 1
END=$(date)
echo get_result_bootstrap-fdr-0.2-none.sh,$START,$END,$SECONDS >> timing/bootstrap_result_time.csv
