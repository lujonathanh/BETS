set -e
START=$(date)
time python3 integrate_outputs_rand_row.py -i bootstrap/56/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-56_output_matr_list.txt -o bootstrap/56/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-56_int_matr_list.txt -t m -a 1  && time python3 integrate_outputs_rand_row.py -i bootstrap/56/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-56_output_df_list.txt -t d -o bootstrap/56/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-56_int_df_list.txt
time python3 get_result_coef.py -df reps.txt -rdf reps-rand.txt -lr 1 -bh hyper/best_hyper.p -o insilico_size100_1-0mean-reps-enet-2-g-bootstrap-56 -cf bootstrap/56/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-56_coefs.p -if bootstrap/56/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-56_intercepts.p -cfr bootstrap/56/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-56_coefsr.p -fr bootstrap/56/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-56_fit_result_df.txt -frr bootstrap/56/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-56_fit_result_dfr.txt -l 2 -sb e -tn enet -of bootstrap/56
END=$(date)
echo bootstrap-finish-scripts/effect/finish-effect-bootstrap-56.sh,$START,$END,$SECONDS >> timing/bootstrap_finish_time.csv
