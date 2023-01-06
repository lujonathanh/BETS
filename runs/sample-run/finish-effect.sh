#!/bin/bash
START=$(date)
set -e
time python3 integrate_outputs_rand_row.py -i output_matr_list.txt -o int_matr_list.txt  -t a  && time python3 integrate_outputs_rand_row.py -i output_df_list.txt -o int_df_list.txt -t d 
time python3 get_result_coef.py -df reps.txt -rdf reps-rand.txt -lr 1 -bh hyper/best_hyper.p -o insilico_size100_1-0mean-reps-enet-2-g -cf fit/insilico_size100_1-0mean-reps-enet-2-g_coefs.p -if fit/insilico_size100_1-0mean-reps-enet-2-g_intercepts.p -cfr fit/insilico_size100_1-0mean-reps-enet-2-g_coefsr.p -fr fit/insilico_size100_1-0mean-reps-enet-2-g_fit_result_df.txt -frr fit/insilico_size100_1-0mean-reps-enet-2-g_fit_result_dfr.txt -l 2 -sb e -tn enet
END=$(date)
echo finish-effect.sh,$START,$END,$SECONDS >> timing/result_time.csv
