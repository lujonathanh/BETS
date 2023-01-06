set -e
START=$(date)
time python3 integrate_outputs_rand_row.py -i bootstrap/12/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-12_output_matr_list.txt -o bootstrap/12/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-12_int_matr_list.txt -t m -a 1  && time python3 integrate_outputs_rand_row.py -i bootstrap/12/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-12_output_df_list.txt -t d -o bootstrap/12/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-12_int_df_list.txt
time python3 get_result_coef.py -df reps.txt -rdf reps-rand.txt -lr 1 -bh hyper/best_hyper.p -o insilico_size100_1-0mean-reps-enet-2-g-bootstrap-12 -cf bootstrap/12/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-12_coefs.p -if bootstrap/12/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-12_intercepts.p -cfr bootstrap/12/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-12_coefsr.p -fr bootstrap/12/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-12_fit_result_df.txt -frr bootstrap/12/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-12_fit_result_dfr.txt -l 2 -sb n -tn enet -of bootstrap/12
END=$(date)
echo bootstrap-finish-scripts/none/finish-none-bootstrap-12.sh,$START,$END,$SECONDS >> timing/bootstrap_finish_time.csv