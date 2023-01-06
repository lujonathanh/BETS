set -e
START=$(date)
time python3 integrate_outputs_rand_row.py -i bootstrap/40/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-40_output_matr_list.txt -o bootstrap/40/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-40_int_matr_list.txt -t m -a 1  && time python3 integrate_outputs_rand_row.py -i bootstrap/40/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-40_output_df_list.txt -t d -o bootstrap/40/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-40_int_df_list.txt
time python3 get_result_coef.py -df reps.txt -rdf reps-rand.txt -lr 1 -bh hyper/best_hyper.p -o insilico_size100_1-0mean-reps-enet-2-g-bootstrap-40 -cf bootstrap/40/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-40_coefs.p -if bootstrap/40/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-40_intercepts.p -cfr bootstrap/40/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-40_coefsr.p -fr bootstrap/40/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-40_fit_result_df.txt -frr bootstrap/40/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-40_fit_result_dfr.txt -l 2 -sb e -tn enet -of bootstrap/40
END=$(date)
echo bootstrap-finish-scripts/effect/finish-effect-bootstrap-40.sh,$START,$END,$SECONDS >> timing/bootstrap_finish_time.csv
