#!/bin/bash
module load python/2.7
module load python/2.7/scipy-mkl
module load python/2.7/numpy-mkl
module load anaconda
time python cross_validate.py -d small_er-reg-cutoff-2-plus-GR-0mean1var-reps.txt -lr 1 -o hyper/test-new-union-method-er-norm-lasso-2-cv-5-row-0-result.txt -hl hyper/test-new-union-method-er-norm-lasso-2-hyper-5.p -t l -l 2 -rl test-new-union-method-er-norm-lasso-2-row-0.txt
