#!/bin/bash
module load python/2.7
module load python/2.7/scipy-mkl
module load python/2.7/numpy-mkl
module load anaconda
time python cross_validate.py -d small_er-reg-cutoff-2-plus-GR-0mean1var-reps.txt -lr 1 -o hyper/test-confirm-old-result-parallelize-er-norm-ridge-2-g-cv-2-row-0-result.txt -hl hyper/test-confirm-old-result-parallelize-er-norm-ridge-2-g-hyper-2.p -t r -l 2 -rl test-confirm-old-result-parallelize-er-norm-ridge-2-g-row-0.p
