__author__ = 'jlu96'

import sys
import geneTSmunging as gtm
import itertools
import pickle
import os
import pandas as pd
import numpy as np
import gc


def partition_inputs(input, number):
    num_inputs = len(input)
    return [input[int(num_inputs * i/number):int(num_inputs * (i+1)/number)] for i in range(number)]

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

    indexlist = list(range(index1, index2 + 1))
    pairlist = [tuple(convert_index_to_pair(index, row_length)) for index in indexlist]

    return pairlist


def concatenate_matrices(filenames, concatenated_filename=None, axis=-1,
                  require_full=True):
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
            all_matrs.append(pickle.load(open(filename, 'rb')))
        else:
            print("Error", filename, "missing")
            if require_full:
                print("Stopping since missing a file")
                return

    concatenated_matr = np.concatenate(tuple(all_matrs), axis=axis)

    print("Final matrix shape", concatenated_matr.shape)

    if concatenated_filename != None:
        with open(concatenated_filename, 'wb') as outfile:
            pickle.dump(concatenated_matr, outfile)

        print("Concatenated matr written to ", concatenated_filename)

    return concatenated_matr


def add_matrices(filenames, added_filename=None, axis=0,
                 iter_free=True):
    """
    :param filenames: filenames of numpy matrices
    :param added_filename: matrix to write filename to
    :param in the list of matrices, the axis to add along. Should be 0.
    :param iter_free: iteratively add and free the filenames instead of adding all at once.
    :return: added_matr
    """


    if iter_free:
        added_matr = pickle.load(open(filenames[0], 'rb'))

        print("Integrating by iterative adding + freeing")
        for i in range(1, len(filenames)):
            if i % 10 == 0:
                print(i)
            next_matr = pickle.load(open(filenames[i], 'rb'))
            added_matr += next_matr

            del next_matr


    else:
        all_matrs = []
        for i, filename in enumerate(filenames):
            if i % 10 == 0:
                print(i)
            if os.path.exists(filename):
                all_matrs.append(pickle.load(open(filename, 'rb')))
            else:
                print("Error, missing: ", filename)

        added_matr = np.sum(all_matrs, axis=axis)

    print("Final matrix shape", added_matr.shape)

    if added_filename != None:
        with open(added_filename, 'wb') as outfile:
            pickle.dump(added_matr, outfile)

        print("Added matr written to ", added_filename)

    return added_matr



def integrate_dfs(filenames, integrated_filename=None,
                  require_full=True):
    all_dfs = []
    for filename in filenames:
        if os.path.exists(filename):
            print(filename)
            all_dfs.append(pd.read_csv(filename, sep="\t"))
        else:
            print("Error: ", filename, " does not exist")

            if require_full:
                # print
                raise ValueError("Stopping since missing a file")

    if all_dfs:
        integrated_df = pd.concat(all_dfs)
    else:
        # Submit and finish NOT DONE 12/19/16
        integrated_df =  None
        print("No DF found, nothing written")
        return

    print("Final df shape:", integrated_df.shape)

    if integrated_filename != None:
        integrated_df.to_csv(integrated_filename, sep="\t", index=False)

        print("Integrated df written to ", integrated_filename)

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
    raise ValueError("This main function is deprecated, do not run!")

def main():
    run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()
