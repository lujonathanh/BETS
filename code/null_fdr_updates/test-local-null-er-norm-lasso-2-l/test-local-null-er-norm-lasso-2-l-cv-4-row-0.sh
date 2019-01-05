#!/bin/bash
module load python/2.7
module load python/2.7/scipy-mkl
module load python/2.7/numpy-mkl
module load anaconda
time python cross_validate.py -d small_er-reg-cutoff-2-plus-GR-0mean1var-reps.txt -lr 1 -o hyper/test-local-null-er-norm-lasso-2-l-cv-4-row-0-result.txt -hl hyper/test-local-null-er-norm-lasso-2-l-hyper-4.p -t l -l 2 -rl test-local-null-er-norm-lasso-2-l-row-0.txt
