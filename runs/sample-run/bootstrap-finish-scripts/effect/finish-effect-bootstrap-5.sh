set -e
START=$(date)
time python3 integrate_outputs_rand_row.py -i bootstrap/5/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-5_output_matr_list.txt -o bootstrap/5/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-5_int_matr_list.txt -t m -a 1  && time python3 integrate_outputs_rand_row.py -i bootstrap/5/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-5_output_df_list.txt -t d -o bootstrap/5/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-5_int_df_list.txt
time python3 get_result_coef.py -df reps.txt -rdf reps-rand.txt -lr 1 -bh hyper/best_hyper.p -o insilico_size100_1-0mean-reps-enet-2-g-bootstrap-5 -cf bootstrap/5/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-5_coefs.p -if bootstrap/5/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-5_intercepts.p -cfr bootstrap/5/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-5_coefsr.p -fr bootstrap/5/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-5_fit_result_df.txt -frr bootstrap/5/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-5_fit_result_dfr.txt -l 2 -sb e -tn enet -of bootstrap/5
END=$(date)
echo bootstrap-finish-scripts/effect/finish-effect-bootstrap-5.sh,$START,$END,$SECONDS >> timing/bootstrap_finish_time.csv
