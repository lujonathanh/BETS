set -e
START=$(date)
time python3 integrate_outputs_rand_row.py -i bootstrap/46/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-46_output_matr_list.txt -o bootstrap/46/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-46_int_matr_list.txt -t m -a 1  && time python3 integrate_outputs_rand_row.py -i bootstrap/46/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-46_output_df_list.txt -t d -o bootstrap/46/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-46_int_df_list.txt
time python3 get_result_coef.py -df reps.txt -rdf reps-rand.txt -lr 1 -bh hyper/best_hyper.p -o insilico_size100_1-0mean-reps-enet-2-g-bootstrap-46 -cf bootstrap/46/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-46_coefs.p -if bootstrap/46/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-46_intercepts.p -cfr bootstrap/46/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-46_coefsr.p -fr bootstrap/46/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-46_fit_result_df.txt -frr bootstrap/46/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-46_fit_result_dfr.txt -l 2 -sb e -tn enet -of bootstrap/46
END=$(date)
echo bootstrap-finish-scripts/effect/finish-effect-bootstrap-46.sh,$START,$END,$SECONDS >> timing/bootstrap_finish_time.csv
