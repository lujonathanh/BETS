#!/bin/bash
module load python/2.7
module load python/2.7/scipy-mkl
module load python/2.7/numpy-mkl
module load anaconda
time python fit_all.py -d small_er-reg-cutoff-2-plus-GR-0mean1var-reps.txt -rd small_er-reg-cutoff-2-plus-GR-0mean1var-reps-rand.txt -lr 1 -o fit/test-confirm-old-result-er-norm-ridge-2-g-fit-row-0 -bh hyper/best_hyper.p -t r -l 2 -rl test-confirm-old-result-er-norm-ridge-2-g-row-0.p -n g
