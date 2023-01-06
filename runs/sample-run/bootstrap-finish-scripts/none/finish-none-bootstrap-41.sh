set -e
START=$(date)
time python3 integrate_outputs_rand_row.py -i bootstrap/41/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-41_output_matr_list.txt -o bootstrap/41/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-41_int_matr_list.txt -t m -a 1  && time python3 integrate_outputs_rand_row.py -i bootstrap/41/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-41_output_df_list.txt -t d -o bootstrap/41/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-41_int_df_list.txt
time python3 get_result_coef.py -df reps.txt -rdf reps-rand.txt -lr 1 -bh hyper/best_hyper.p -o insilico_size100_1-0mean-reps-enet-2-g-bootstrap-41 -cf bootstrap/41/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-41_coefs.p -if bootstrap/41/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-41_intercepts.p -cfr bootstrap/41/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-41_coefsr.p -fr bootstrap/41/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-41_fit_result_df.txt -frr bootstrap/41/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-41_fit_result_dfr.txt -l 2 -sb n -tn enet -of bootstrap/41
END=$(date)
echo bootstrap-finish-scripts/none/finish-none-bootstrap-41.sh,$START,$END,$SECONDS >> timing/bootstrap_finish_time.csv
