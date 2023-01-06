__author__ = 'jlu96'

import prep_jobs as pj
import sys
import os
import pandas as pd
import math
import collections
import geneTSmunging as gtm
import csv

def get_parser():
    # Parse arguments
    import argparse

    description = 'Prepare cluster jobs by partitioning tests by rows and hyper-parameters.'
    parser = argparse.ArgumentParser(description=description)


    parser.add_argument('-d', '--data_file', required=True)

    parser.add_argument('-lr', '--load_reps', required=True, type=int)

    parser.add_argument('-o', '--output_name', required=True)

    parser.add_argument('-sn', '--script_num', type=int, default=3,
                        help="Number of scripts, total, that need to be run")

    parser.add_argument('-p', '--parallel_num', type=int, default=0)

    parser.add_argument('-l', '--lag', type=int, required=True)

    return parser

def run(args):

    data_file = args.data_file

    if args.load_reps:
        genes, geneTS = gtm.load_basic_rep_file_list(data_file)
        #dfs, genes, geneTS, df, __, __ = gtm.load_rep_file_list(data_file)
    else:
        df = pd.read_csv(data_file, sep="\t")
        genes, geneTS = gtm.get_gene_TS(df)
    n = len(genes)




    # Make row files
    # Split up the rows according to number of input scripts
    partition_rows = pj.partition_inputs(list(range(n)), args.script_num)

    row_filenames = []


    print("*************")
    print("ROWS")
    print("*************")

    for partition_row, i in zip(partition_rows, list(range(len(partition_rows)))):

        row_filename = os.path.join("rows", args.output_name + "-row-" + str(i) + ".p")
        row_filenames.append(row_filename)

    print("Reading rows from format: ", row_filename)

    print("*************")
    print("PAIRWISE")
    print("*************")


    # Run the actual fit
    # Need an integration
    if not os.path.exists("pairwise"):
        os.makedirs("pairwise")

    # For the pairwise individual fit scripts
    if not os.path.exists("pairwise-fit-scripts"):
        os.makedirs("pairwise-fit-scripts")


    # For the pairwise finish scripts
    if not os.path.exists("pairwise-finish-scripts"):
        os.makedirs("pairwise-finish-scripts")


    pairwise_result_folder = os.path.join("pairwise", "pairwise-results")
    if not os.path.exists(pairwise_result_folder):
        os.makedirs(pairwise_result_folder)





    # make one script for each...

    # all_bootstrap_scripts = set([])

    # all_int_coefs = []
    # all_int_intercepts = []

    # record where the thresholded coefficients are written
    # For integrating these, later.



    try:
        fittimefile = os.path.join("timing", "pairwise_fit_time.csv")
        if not os.path.exists(fittimefile):
            with open(fittimefile, 'w') as csvfile:
                f = csv.writer(csvfile)
                f.writerow(["Name", "Start", "End", "Elapsed"])


        finishtimefile = os.path.join("timing", "pairwise_finish_time.csv")
        if not os.path.exists(finishtimefile):
            with open(finishtimefile, 'w') as csvfile:
                f = csv.writer(csvfile)
                f.writerow(["Name", "Start", "End", "Elapsed"])

        # resulttimefile = os.path.join("timing", "bootstrap_result_time.csv")
        # if not os.path.exists(resulttimefile):
        #     with open(resulttimefile, 'w') as csvfile:
        #         f = csv.writer(csvfile)
        #         f.writerow(["Name", "Start", "End", "Elapsed"])

        with open(os.path.join("timing/timing_list.txt"), 'a') as f:
            f.write(fittimefile + "\n")
            f.write(finishtimefile + "\n")
            # f.write(resulttimefile + "\n")


    except IOError:
        raise IOError("the timing folder does not exist. Please run ./prep_jobs_rand_cv.sh first.")

    pairwise_outmost_name = args.output_name + "-pairwise"
    pairwise_outmost_prefix = os.path.join("pairwise", pairwise_outmost_name)


    # create scripts for pairwise
    pairwise_scripts = [os.path.join("pairwise-fit-scripts", pairwise_outmost_name + "-row-" + str(i) + ".sh")
                         for i in range(len(partition_rows))]
    pairwise_row_prefixes = [pairwise_outmost_prefix + "-row-" + str(i) for i in range(len(partition_rows))]

    command_template = "time python3 fit_pairwise.py -d " + data_file + " -lr " + str(args.load_reps) + \
                         " -o " + "pairwise_row_prefixes[i]" +  " -l " + str(args.lag) + " -rl " + \
                         "row_filename"

    for i, row_filename in zip(list(range(len(partition_rows))), row_filenames):

        # writing results to the pairwise prefix

        command_string = command_template.replace("pairwise_row_prefixes[i]", pairwise_row_prefixes[i]).replace("row_filename", row_filename)

        with open(pairwise_scripts[i], 'w') as outputfile:
                outputfile.write("#!/bin/bash\n")
                outputfile.write("START=$(date)\n")
                outputfile.write("module load python/2.7\n")
                # outputfile.write("module load python/2.7/scipy-mkl\n")
                # outputfile.write("module load python/2.7/numpy-mkl\n")
                outputfile.write("module load anaconda\n")
                outputfile.write(command_string)
                outputfile.write("\n")
                outputfile.write("END=$(date)\n")
                outputfile.write("echo " + pairwise_scripts[i] + ",$START,$END,$SECONDS >> " + fittimefile + "\n")
        os.chmod(pairwise_scripts[i], 0o777)


        print("Scripts made")

        # all_pairwise_scripts = all_pairwise_scripts.union(set(pairwise_scripts))

        # Note the output files

    pairwise_coefs = [pairwise_row_prefix + "_coefs.p" for pairwise_row_prefix in pairwise_row_prefixes]

    pairwise_output_dict = collections.OrderedDict()
    pairwise_output_dict["coef"] = pairwise_coefs

    output_matr_df = pd.DataFrame(pairwise_output_dict)
    output_matr_file = os.path.join("pairwise", pairwise_outmost_name + "_output_matr_list.txt")
    output_matr_df.to_csv(output_matr_file, sep="\t", index=False)
    print("Raw parallelilized output matrices, before integration, written to", output_matr_file)




    int_matr_dict = collections.OrderedDict()
    int_matr_dict["coef"] = os.path.join(pairwise_result_folder, pairwise_outmost_name + "_coefs.p")

    # # append these to the list of final bootstrapped coefficients
    # all_int_coefs.append(int_matr_dict["coef"])
    # all_int_intercepts.append(int_matr_dict["intercept"])

    int_matr_file = pairwise_outmost_prefix +  "_int_matr_list.txt"
    int_matr_df = pd.DataFrame(int_matr_dict, index=[0])
    int_matr_df.to_csv(int_matr_file, sep="\t", index=False)
    print("integrated matrices written to " + int_matr_file)



    # just need to put all of this into the outmost name

    all_pairwise_scripts = [os.path.join("pairwise-fit-scripts", pairwise_outmost_name + "-row-" + str(i) + ".sh")
                             for i in range(len(partition_rows))]


    print("SCRIPTS")

    with open("pairwise_script_list.txt", 'w') as outfile:
        for pairwise_script in all_pairwise_scripts:
            outfile.write("./" + pairwise_script + "\n")
        print("pairwise scripts written to pairwise_script_list.txt")

        if args.parallel_num > 0:
            print("Parallel Number (# processes per job): " + str(args.parallel_num))

            script_groups = pj.partition_inputs(all_pairwise_scripts, number=int(math.ceil(len(all_pairwise_scripts) * 1.0/args.parallel_num)))

            print("Number of script groups ", len(script_groups))

            parallel_scripts = []
            for i, script_group in zip(list(range(len(script_groups))), script_groups):
                appended_script_filenames = ["./" + script_filename for script_filename in script_group]
                parallel_script = " & ".join(appended_script_filenames)
                parallel_scripts.append(parallel_script)

            with open("pairwise_parallel_script_list.txt", 'w') as scriptfile:
                for parallel_script in parallel_scripts:
                    scriptfile.write(parallel_script + "\n")
                print("Parallel script list written to pairwise_parallel_script_list.txt")


    finish_script = os.path.join("pairwise-finish-scripts", "finish.sh")
    with open(finish_script, 'w') as ifile:
        ifile.write("set -e\n")
        ifile.write("START=$(date)\n")
        ifile.write("time python3 integrate_outputs_rand_row.py -i " + output_matr_file + " -o " + int_matr_file +  " -t a \n")
        ifile.write("END=$(date)\n")
        ifile.write("echo " + finish_script + ",$START,$END,$SECONDS >> " + finishtimefile + "\n")
        print("Finish script, written to", finish_script)
        os.chmod(finish_script, 0o777)


def main():
    run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()
