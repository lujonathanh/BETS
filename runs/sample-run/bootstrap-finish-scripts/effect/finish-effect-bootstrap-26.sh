set -e
START=$(date)
time python3 integrate_outputs_rand_row.py -i bootstrap/26/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-26_output_matr_list.txt -o bootstrap/26/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-26_int_matr_list.txt -t m -a 1  && time python3 integrate_outputs_rand_row.py -i bootstrap/26/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-26_output_df_list.txt -t d -o bootstrap/26/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-26_int_df_list.txt
time python3 get_result_coef.py -df reps.txt -rdf reps-rand.txt -lr 1 -bh hyper/best_hyper.p -o insilico_size100_1-0mean-reps-enet-2-g-bootstrap-26 -cf bootstrap/26/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-26_coefs.p -if bootstrap/26/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-26_intercepts.p -cfr bootstrap/26/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-26_coefsr.p -fr bootstrap/26/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-26_fit_result_df.txt -frr bootstrap/26/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-26_fit_result_dfr.txt -l 2 -sb e -tn enet -of bootstrap/26
END=$(date)
echo bootstrap-finish-scripts/effect/finish-effect-bootstrap-26.sh,$START,$END,$SECONDS >> timing/bootstrap_finish_time.csv
