#!/bin/bash
START=$(date)
time python3 get_result_coef.py -df reps.txt -rdf reps-rand.txt -lr 1 -bh hyper/best_hyper.p -o insilico_size100_1-0mean-reps-enet-2-g -cf fit/insilico_size100_1-0mean-reps-enet-2-g_coefs.p -if fit/insilico_size100_1-0mean-reps-enet-2-g_intercepts.p -cfr fit/insilico_size100_1-0mean-reps-enet-2-g_coefsr.p -fr fit/insilico_size100_1-0mean-reps-enet-2-g_fit_result_df.txt -frr fit/insilico_size100_1-0mean-reps-enet-2-g_fit_result_dfr.txt -l 2 -sb n -tn enet -pcf 1 
END=$(date)
echo plot_coef.sh,$START,$END,$SECONDS >> timing/result_time.csv
