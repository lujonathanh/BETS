set -e
START=$(date)
time python3 integrate_outputs_rand_row.py -i bootstrap/6/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-6_output_matr_list.txt -o bootstrap/6/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-6_int_matr_list.txt -t m -a 1  && time python3 integrate_outputs_rand_row.py -i bootstrap/6/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-6_output_df_list.txt -t d -o bootstrap/6/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-6_int_df_list.txt
time python3 get_result_coef.py -df reps.txt -rdf reps-rand.txt -lr 1 -bh hyper/best_hyper.p -o insilico_size100_1-0mean-reps-enet-2-g-bootstrap-6 -cf bootstrap/6/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-6_coefs.p -if bootstrap/6/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-6_intercepts.p -cfr bootstrap/6/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-6_coefsr.p -fr bootstrap/6/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-6_fit_result_df.txt -frr bootstrap/6/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-6_fit_result_dfr.txt -l 2 -sb e -tn enet -of bootstrap/6
END=$(date)
echo bootstrap-finish-scripts/effect/finish-effect-bootstrap-6.sh,$START,$END,$SECONDS >> timing/bootstrap_finish_time.csv
