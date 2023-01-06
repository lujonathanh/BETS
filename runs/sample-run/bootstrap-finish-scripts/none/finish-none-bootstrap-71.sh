set -e
START=$(date)
time python3 integrate_outputs_rand_row.py -i bootstrap/71/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-71_output_matr_list.txt -o bootstrap/71/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-71_int_matr_list.txt -t m -a 1  && time python3 integrate_outputs_rand_row.py -i bootstrap/71/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-71_output_df_list.txt -t d -o bootstrap/71/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-71_int_df_list.txt
time python3 get_result_coef.py -df reps.txt -rdf reps-rand.txt -lr 1 -bh hyper/best_hyper.p -o insilico_size100_1-0mean-reps-enet-2-g-bootstrap-71 -cf bootstrap/71/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-71_coefs.p -if bootstrap/71/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-71_intercepts.p -cfr bootstrap/71/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-71_coefsr.p -fr bootstrap/71/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-71_fit_result_df.txt -frr bootstrap/71/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-71_fit_result_dfr.txt -l 2 -sb n -tn enet -of bootstrap/71
END=$(date)
echo bootstrap-finish-scripts/none/finish-none-bootstrap-71.sh,$START,$END,$SECONDS >> timing/bootstrap_finish_time.csv
