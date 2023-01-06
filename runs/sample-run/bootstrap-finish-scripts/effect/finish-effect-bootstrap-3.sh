set -e
START=$(date)
time python3 integrate_outputs_rand_row.py -i bootstrap/3/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-3_output_matr_list.txt -o bootstrap/3/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-3_int_matr_list.txt -t m -a 1  && time python3 integrate_outputs_rand_row.py -i bootstrap/3/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-3_output_df_list.txt -t d -o bootstrap/3/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-3_int_df_list.txt
time python3 get_result_coef.py -df reps.txt -rdf reps-rand.txt -lr 1 -bh hyper/best_hyper.p -o insilico_size100_1-0mean-reps-enet-2-g-bootstrap-3 -cf bootstrap/3/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-3_coefs.p -if bootstrap/3/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-3_intercepts.p -cfr bootstrap/3/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-3_coefsr.p -fr bootstrap/3/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-3_fit_result_df.txt -frr bootstrap/3/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-3_fit_result_dfr.txt -l 2 -sb e -tn enet -of bootstrap/3
END=$(date)
echo bootstrap-finish-scripts/effect/finish-effect-bootstrap-3.sh,$START,$END,$SECONDS >> timing/bootstrap_finish_time.csv
