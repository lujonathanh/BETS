#!/bin/bash
START=$(date)
module load anaconda3
time python3 fit_bootstrap.py -d reps.txt -rd reps-rand.txt -lr 1 -o bootstrap/54/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-54-row-0 -bh hyper/best_hyper.p -t e -l 2 -rl rows/insilico_size100_1-0mean-reps-enet-2-g-row-0.p -n g -s 54 -oa 1
END=$(date)
echo bootstrap-fit-scripts/54/insilico_size100_1-0mean-reps-enet-2-g-bootstrap-54-row-0.sh,$START,$END,$SECONDS >> timing/bootstrap_fit_time.csv
