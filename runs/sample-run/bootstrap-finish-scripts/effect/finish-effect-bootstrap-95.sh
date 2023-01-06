set -e
START=$(date)
time python3 integrate_outputs_rand_row.py -i bootstrap/95/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-95_output_matr_list.txt -o bootstrap/95/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-95_int_matr_list.txt -t m -a 1  && time python3 integrate_outputs_rand_row.py -i bootstrap/95/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-95_output_df_list.txt -t d -o bootstrap/95/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-95_int_df_list.txt
time python3 get_result_coef.py -df reps.txt -rdf reps-rand.txt -lr 1 -bh hyper/best_hyper.p -o insilico_size100_1-0mean-reps-enet-2-g-bootstrap-95 -cf bootstrap/95/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-95_coefs.p -if bootstrap/95/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-95_intercepts.p -cfr bootstrap/95/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-95_coefsr.p -fr bootstrap/95/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-95_fit_result_df.txt -frr bootstrap/95/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-95_fit_result_dfr.txt -l 2 -sb e -tn enet -of bootstrap/95
END=$(date)
echo bootstrap-finish-scripts/effect/finish-effect-bootstrap-95.sh,$START,$END,$SECONDS >> timing/bootstrap_finish_time.csv
