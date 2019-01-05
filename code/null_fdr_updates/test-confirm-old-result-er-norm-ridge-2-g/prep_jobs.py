__author__ = 'jlu96'

import sys
import geneTSmunging as gtm
import itertools
import pickle
import os
import pandas as pd
import numpy as np


def partition_inputs(input, number):
    num_inputs = len(input)
    return [input[num_inputs * i/number:num_inputs * (i+1)/number] for i in range(number)]

def lazy_partition_indices(length, number):
    return [(i * length/ number, (i + 1) * length/ number - 1) for i in range(number)]

def convert_index_to_pair(index, row_length):
    return (index / row_length, index % row_length)

def convert_pair_to_index(pair, row_length):
    return pair[0] * row_length + pair[1]

def lazy_partition_pairs(n_rows, n_cols, number):
    return [(convert_index_to_pair(index1, row_length=n_cols), convert_index_to_pair(index2, row_length=n_cols))
        for index1, index2 in lazy_partition_indices(n_rows * n_cols, number)]



def partition_pair_to_pairlist(partition_pair, row_length):
    pair1, pair2 = partition_pair

    index1 = convert_pair_to_index(pair1, row_length)
    index2 = convert_pair_to_index(pair2, row_length)

    indexlist = range(index1, index2 + 1)
    pairlist = [tuple(convert_index_to_pair(index, row_length)) for index in indexlist]

    return pairlist


def concatenate_matrices(filenames, concatenated_filename=None, axis=-1):
    """
    Concatenate matrices and write out.
    :param filenames: filenames of numpy matrices to load
    :param concatenated_filename: filename of output
    :param axis: axis along which to concatenate

    :return: concatenated_matr
    """
    if axis == -1:
            raise ValueError("You must set the axis of concatenation for the matrices.")

    all_matrs = []
    for filename in filenames:
        if os.path.exists(filename):
            all_matrs.append(pickle.load(open(filename, 'rU')))
        else:
            print "Error", filename, "missing"

    concatenated_matr = np.concatenate(tuple(all_matrs), axis=axis)

    print "Final matrix shape", concatenated_matr.shape

    if concatenated_filename != None:
        with open(concatenated_filename, 'w') as outfile:
            pickle.dump(concatenated_matr, outfile)

        print "Concatenated matr written to ", concatenated_filename

    return concatenated_matr


def add_matrices(filenames, added_filename=None, axis=0):
    """
    :param filenames: filenames of numpy matrices
    :param added_filename: matrix to write filename to
    :param in the list of matrices, the axis to add along. Should be 0.
    :return: added_matr
    """

    all_matrs = []
    for filename in filenames:

        if os.path.exists(filename):
            all_matrs.append(pickle.load(open(filename, 'rU')))
        else:
            print "Error, missing: ", filename

    added_matr = np.sum(all_matrs, axis=axis)

    print "Final matrix shape", added_matr.shape

    if added_filename != None:
        with open(added_filename, 'w') as outfile:
            pickle.dump(added_matr, outfile)

        print "Added matr written to ", added_filename

    return added_matr



def integrate_dfs(filenames, integrated_filename=None,
                  require_full=True):
    all_dfs = []
    for filename in filenames:
        if os.path.exists(filename):
            all_dfs.append(pd.read_csv(filename, sep="\t"))
        else:
            print "Error: ", filename, " does not exist"

            if require_full:
                print "Stopping since missing a file"
                return

    if all_dfs:
        integrated_df = pd.concat(all_dfs)
    else:
        # Submit and finish NOT DONE 12/19/16
        integrated_df =  None
        print "No DF found, nothing written"
        return

    print "Final df shape:", integrated_df.shape

    if integrated_filename != None:
        integrated_df.to_csv(integrated_filename, sep="\t", index=False)

        print "Integrated df written to ", integrated_filename

    return integrated_df


def get_parser():
    # Parse arguments
    import argparse

    description = 'Apply a pre-specified causal test to an input dataset where each row is a gene' \
                  'and its time points'
    parser = argparse.ArgumentParser(description=description)


    parser.add_argument('-d', '--data_file', required=True)

    parser.add_argument('-a', '--args_file', required=True)

    parser.add_argument('-t', '--test', required=True)

    parser.add_argument('-o', '--output_name', required=True)

    parser.add_argument('-n', '--job_num', type=int, default=3)

    parser.add_argument('-p', '--parallel_num', type=int, default=0)

    parser.add_argument('-d2', '--data_file2', default=None, help="The effect genes")

    return parser

def run(args):

    df = gtm.load_file_and_avg(args.data_file)

    genes = df['gene'].values

    n = len(genes)


    script_filenames = []
    output_filenames = []


    partition_pairs = lazy_partition_pairs(n, n, args.job_num)

    for partition_pair, i in zip(partition_pairs, range(len(partition_pairs))):

        script_filename = args.output_name + "-script-" + str(i) + ".sh"
        script_filenames.append(script_filename)


        output_filename = args.output_name + "-" + str(i) + ".p"
        output_filenames.append(output_filename)
        # prepare the job associated with this

        pair_filename = args.output_name + "-pair-" + str(i) + ".txt"

        command_string = "python run_causal.py -d " + args.data_file.split('/')[-1] + " -a " + args.args_file.split('/')[-1] + " -t " + args.test + " -pp " + \
                         str(pair_filename) + " -o " + output_filename


        if args.test == "gp":
            command_string += " -d2 " + args.data_file2.split('/')[-1]

        with open(pair_filename, 'w') as pairfile:
            pairfile.write(str(partition_pair) + "\n")

        print "Partition pair written to ", pair_filename


        with open(script_filename, 'w') as outputfile:
            outputfile.write("#!/bin/bash\n")
            outputfile.write("module load python/2.7\n")
            outputfile.write("module load python/2.7/scipy-mkl\n")
            outputfile.write("module load python/2.7/numpy-mkl\n")
            outputfile.write("module load anaconda\n")
            outputfile.write(command_string)
            outputfile.write("\n")
        os.chmod(script_filename, 0777)

        print "Script written to ", script_filename

    # submit the jobs soon


    with open("script_list.txt", 'w') as scriptfile:
        for script_filename in script_filenames:
            scriptfile.write(script_filename + "\n")
        print "Script list written to script_list.txt"

    with open("output_list.txt", 'w') as outputfile:
        for output_filename in output_filenames:
            outputfile.write(output_filename + "\n")
        print "Output list written to output_list.txt"

    with open("integrate_outputs.sh", 'w') as ifile:
        integrated_filename = args.output_name + ".p"
        ifile.write("python integrate_outputs.py -i output_list.txt -o " + integrated_filename + " -n " + str(n) + "\n")
        print "Integration script written to integrate_outputs.sh"
        os.chmod("integrate_outputs.sh", 0777)


    if args.parallel_num > 0:
        print "Parallel Number (# processes per job): " + str(args.parallel_num)

        script_groups = partition_inputs(script_filenames, number=len(script_filenames)/args.parallel_num)

        print "Number of script groups ", len(script_groups)


        parallel_scripts = []
        for i, script_group in zip(range(len(script_groups)), script_groups):
            appended_script_filenames = ["./" + script_filename for script_filename in script_group]
            parallel_script = " & ".join(appended_script_filenames)
            print "Parallel Script ", i, ":", parallel_script
            parallel_scripts.append(parallel_script)

        with open("parallel_script_list.txt", 'w') as scriptfile:
            for parallel_script in parallel_scripts:
                scriptfile.write(parallel_script + "\n")
            print "Parallel script list written to parallel_script_list.txt"


def main():
    run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()


# Old version
#
# import sys
# import geneTSmunging as gtm
# import itertools
# import pickle
# import os
#
#
# def partition_inputs(input, number):
#     num_inputs = len(input)
#     return [input[num_inputs * i/number:num_inputs * (i+1)/number] for i in range(number)]
#
# def get_parser():
#     # Parse arguments
#     import argparse
#
#     description = 'Apply a pre-specified causal test to an input dataset where each row is a gene' \
#                   'and its time points'
#     parser = argparse.ArgumentParser(description=description)
#
#
#     parser.add_argument('-d', '--data_file', required=True)
#
#     parser.add_argument('-a', '--args_file', required=True)
#
#     parser.add_argument('-t', '--test', required=True)
#
#     parser.add_argument('-plp', '--pairlist_fileprefix', default=None)
#
#     parser.add_argument('-o', '--output_name', required=True)
#
#     parser.add_argument('-n', '--job_num', type=int, default=3)
#
#     return parser
#
# def run(args):
#
#     df = gtm.load_file_and_avg(args.data_file)
#
#     genes = df['gene'].values
#
#     n = len(genes)
#
#     pairs = list(itertools.permutations(range(n), 2))
#
#     pairlists = partition_inputs(pairs, args.job_num)
#
#     if args.pairlist_fileprefix == None:
#         args.pairlist_fileprefix = args.output_name + "-pairs-"
#
#
#     script_filenames = []
#     output_filenames = []
#
#     for pairlist, i in zip(pairlists, range(len(pairlists))):
#         pairlist_filename = args.pairlist_fileprefix + str(i) + ".p"
#
#         pickle.dump(pairlist, open(pairlist_filename, 'w'))
#
#         print "Pair lists written to", pairlist_filename
#
#
#         script_filename = args.output_name + "-script-" + str(i) + ".sh"
#         script_filenames.append(script_filename)
#
#         output_filename = args.output_name + "-" + str(i) + ".p"
#         output_filenames.append(output_filename)
#         # prepare the job associated with this
#
#         command_string = "python run_causal.py -d " + args.data_file.split('/')[-1] + " -a " + args.args_file.split('/')[-1] + " -t " + args.test + " -pl " + \
#                          pairlist_filename + " -o " + output_filename
#
#         with open(script_filename, 'w') as outputfile:
#             outputfile.write("#!/bin/bash\n")
#             outputfile.write("module load python/2.7\n")
#             outputfile.write("module load python/2.7/scipy-mkl\n")
#             outputfile.write("module load python/2.7/numpy-mkl\n")
#             outputfile.write("module load anaconda\n")
#             outputfile.write(command_string)
#             outputfile.write("\n")
#         os.chmod(script_filename, 0777)
#
#         print "Script written to ", script_filename
#
#     # submit the jobs soon
#
#
#     with open("script_list.txt", 'w') as scriptfile:
#         for script_filename in script_filenames:
#             scriptfile.write(script_filename + "\n")
#         print "Script list written to script_list.txt"
#
#     with open("output_list.txt", 'w') as outputfile:
#         for output_filename in output_filenames:
#             outputfile.write(output_filename + "\n")
#         print "Output list written to output_list.txt"
#
#     with open("integrate_outputs.sh", 'w') as ifile:
#         integrated_filename = args.output_name + ".p"
#         ifile.write("python integrate_outputs.py -i output_list.txt -o " + integrated_filename + "\n")
#         print "Integration script written to integrate_outputs.sh"
#         os.chmod("integrate_outputs.sh", 0777)
#
# def main():
#     run(get_parser().parse_args(sys.argv[1:]))
#
# if __name__ == '__main__':
#     main()
