import os
import pickle
import numpy as np
import sys


def get_parser():
    # Parse arguments
    import argparse

    description = 'Get mean and std of intercepts from bootstraped fits'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-b', '--bootstrap_file_with_names', required=True,
                        help="Name of file with each row as a bootstrap file " +
                             "of intercepts")

    parser.add_argument('-o', '--save_prefix', required=True)

    parser.add_argument('-rsf', '--result_save_folder', required=True)

    return parser


def mean_and_sd_matrices(filenames, axis=0):
    """
    :param filenames: filenames of numpy matrices
    :param added_filename: matrix to write filename to
    :param in the list of matrices, the axis to add along. Should be 0.
    :param iter_free: iteratively add and free the filenames instead of adding all at once.
    :return: added_matr
    """


    all_matrs = []
    prev_matr = []

    for i, filename in enumerate(filenames):
        if i % 10 == 0:
            print(i)
        if os.path.exists(filename):
            new_matr = pickle.load(open(filename, 'rb'))
            if prev_matr != []:
                if not new_matr.shape == prev_matr.shape:
                    raise ValueError("Inconsistent matrix shape: previous was " + str(prev_matr.shape) + " and new is " + str(new_matr.shape))

            prev_matr = new_matr
            all_matrs.append(new_matr)
        else:
            print("Error, missing: ", filename)

    mean_matr = np.mean(all_matrs, axis=axis)
    sd_matr = np.std(all_matrs, axis=axis)

    print("Final matrix shape", mean_matr.shape)

    return mean_matr, sd_matr

def load_and_run(args):

    save_prefix = args.save_prefix

    full_save_prefix = os.path.join(args.result_save_folder, save_prefix)

    filenames = []
    with open(args.bootstrap_file_with_names, 'r') as f:
        for l in f.readlines():
            filenames.append(l.split("\n")[0])

    mean_matr, sd_matr = mean_and_sd_matrices(filenames, axis=0)

    mean_filename = full_save_prefix + "_mean_intercepts.p"
    sd_filename = full_save_prefix + "_sd_intercepts.p"

    pickle.dump(mean_matr, open(mean_filename, 'wb'))
    pickle.dump(sd_matr, open(sd_filename, 'wb'))

    print("Mean written to ", mean_filename)
    print("STD written to ", sd_filename)

def main():
    load_and_run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()