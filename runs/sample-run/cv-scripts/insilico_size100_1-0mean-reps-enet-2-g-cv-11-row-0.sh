#!/bin/bash
START=$(date)
module load anaconda3
time python3 cross_validate.py -d reps.txt -lr 1 -o hyper/insilico_size100_1-0mean-reps-enet-2-g-cv-11-row-0-result.txt -hl hyper/insilico_size100_1-0mean-reps-enet-2-g-hyper-11.p -t e -l 2 -rl rows/insilico_size100_1-0mean-reps-enet-2-g-row-0.p
END=$(date)
echo cv-scripts/insilico_size100_1-0mean-reps-enet-2-g-cv-11-row-0.sh,$START,$END,$SECONDS >> timing/hyper_time.csv
