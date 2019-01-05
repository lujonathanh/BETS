__author__ = 'jlu96'

import os

input_file = "../data/GeneExpressionData/featurecounts.genes.TPM.selected_reps.ln.surrogate_variables_corrected.protein_coding.txt"
name = "Protein-Coding_TPM"
cause_type = "g"
out_file_prefix = name






# INDICES OF GENES TO CONSIDER
start_index = 0 # None
end_index = 100 # None


# GRANGER CAUSALITY PARAMS
p_threshold = 0.05
model_order_min = 2
model_order_max = 3 # max modle order, inclusive

#MULTIPROCESSING
use_processes = True

try:
    procnum = int(os.environ['NPROCS'])
except KeyError:
    procnum = 32