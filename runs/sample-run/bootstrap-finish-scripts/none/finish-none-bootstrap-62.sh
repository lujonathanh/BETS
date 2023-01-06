set -e
START=$(date)
time python3 integrate_outputs_rand_row.py -i bootstrap/62/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-62_output_matr_list.txt -o bootstrap/62/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-62_int_matr_list.txt -t m -a 1  && time python3 integrate_outputs_rand_row.py -i bootstrap/62/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-62_output_df_list.txt -t d -o bootstrap/62/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-62_int_df_list.txt
time python3 get_result_coef.py -df reps.txt -rdf reps-rand.txt -lr 1 -bh hyper/best_hyper.p -o insilico_size100_1-0mean-reps-enet-2-g-bootstrap-62 -cf bootstrap/62/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-62_coefs.p -if bootstrap/62/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-62_intercepts.p -cfr bootstrap/62/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-62_coefsr.p -fr bootstrap/62/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-62_fit_result_df.txt -frr bootstrap/62/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-62_fit_result_dfr.txt -l 2 -sb n -tn enet -of bootstrap/62
END=$(date)
echo bootstrap-finish-scripts/none/finish-none-bootstrap-62.sh,$START,$END,$SECONDS >> timing/bootstrap_finish_time.csv
