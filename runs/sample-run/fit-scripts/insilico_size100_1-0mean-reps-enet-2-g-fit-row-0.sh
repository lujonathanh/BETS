#!/bin/bash
START=$(date)
module load anaconda3
time python3 fit_all.py -d reps.txt -rd reps-rand.txt -lr 1 -o fit/insilico_size100_1-0mean-reps-enet-2-g-fit-row-0 -bh hyper/best_hyper.p -t e -l 2 -rl rows/insilico_size100_1-0mean-reps-enet-2-g-row-0.p -n g -oa 0
END=$(date)
echo fit-scripts/insilico_size100_1-0mean-reps-enet-2-g-fit-row-0.sh,$START,$END,$SECONDS >> timing/fit_time.csv
