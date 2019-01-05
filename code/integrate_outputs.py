__author__ = 'jlu96'


__author__ = 'jlu96'

import sys
import csv
import pickle
import numpy as np


def array_to_matr(array):
    # Assume the list is of shape (nrows * nrows, nobs)


    array = np.array(array)
    length = array.shape[0]
    nrow = np.sqrt(length)


    print "array shape is ", array.shape
    print "nrow is ", nrow

    if len(array.shape) > 1:
        final_shape = (nrow, nrow) + (array.shape[1:])
    else:
        final_shape = (nrow, nrow)

    print "Final shape is ", final_shape

    matr = array.reshape(final_shape)

    return matr




def get_parser():
    # Parse arguments
    import argparse

    description = 'Integrate output to a pickle file.'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-n', '--num_genes', required=True, type=int)

    parser.add_argument('-i', '--output_listname', default="output_list.txt")

    parser.add_argument('-o', '--integrated_filename', required=True)

    parser.add_argument('-p_pos', type=int, default=0)

    parser.add_argument('-p_add', type=int, default=0)

    return parser




def run(args):

    # load the data

    filenames = []
    with open(args.output_listname, 'rU') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            filenames.append(row[0])

    print "Files to integrate are ", filenames

    output_tuples = []

    for filename in filenames:
        output_tuple = pickle.load(open(filename, 'rU'))
        output_tuples.append(output_tuple)

    output_lists = zip(*output_tuples)

    pre_integrated_list = [np.concatenate(output_entry) for output_entry in output_lists]

    integrated_list = []
    for entry in pre_integrated_list:
        integrated_list.append(array_to_matr(entry))

    for entry, i in zip(integrated_list, range(len(integrated_list))):
        print "Shape of entry ", i, " is ", entry.shape

    print "Args.p_add is ", args.p_add
    print "Args.p_add is ", args.p_add

    if args.p_add:
        print "Adding ones to p-matr at location ", args.p_pos
        p_matr = integrated_list[args.p_pos]
        p_matr += np.diagflat(np.ones(p_matr.shape[0]))
        integrated_list = integrated_list[0:args.p_pos] + [p_matr,] + integrated_list[args.p_pos + 1:]

    with open(args.integrated_filename, 'w') as outfile:

        pickle.dump(integrated_list, outfile)

    print "Integrated written to ", args.integrated_filename

def main():
    run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()



# __author__ = 'jlu96'
#
#
# __author__ = 'jlu96'
#
# import sys
# import csv
# import pickle
# import numpy as np
#
#
#
# def get_parser():
#     # Parse arguments
#     import argparse
#
#     description = 'Integrate output to a pickle file.'
#     parser = argparse.ArgumentParser(description=description)
#
#     parser.add_argument('-i', '--output_listname', default="output_list.txt")
#
#     parser.add_argument('-o', '--integrated_filename', required=True)
#
#     parser.add_argument('-p_pos', type=int, default=0)
#
#     parser.add_argument('-p_add', type=bool, default=True)
#
#     return parser
#
#
#
#
# def run(args):
#
#     # load the data
#
#     filenames = []
#     with open(args.output_listname, 'rU') as csvfile:
#         reader = csv.reader(csvfile, delimiter='\t')
#         for row in reader:
#             filenames.append(row[0])
#
#     print "Files to integrate are ", filenames
#
#     output_tuples = []
#
#     for filename in filenames:
#         output_tuple = pickle.load(open(filename, 'rU'))
#         output_tuples.append(output_tuple)
#
#     output_lists = zip(*output_tuples)
#
#     integrated_tuple = tuple(np.sum(output_list, axis=0) for output_list in output_lists)
#
#     if args.p_add:
#         print "Adding ones to p-matr at location ", args.p_pos
#         p_matr = integrated_tuple[args.p_pos]
#         p_matr += np.diagflat(np.ones(p_matr.shape[0]))
#         integrated_tuple = integrated_tuple[0:args.p_pos] + (p_matr,) + integrated_tuple[args.p_pos + 1:]
#
#     with open(args.integrated_filename, 'w') as outfile:
#
#         pickle.dump(integrated_tuple, outfile)
#
#     print "Integrated written to ", args.integrated_filename
#
# def main():
#     run(get_parser().parse_args(sys.argv[1:]))
#
# if __name__ == '__main__':
#     main()
