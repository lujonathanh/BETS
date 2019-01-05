__author__ = 'jlu96'


import numpy as np
import pandas as pd
import pickle
import sys
import collections
import geneTSmunging as gtm
import network_helpers as nh
import causal_pipeline as cp
import lag_conversion as lc
import os
import math
import time



def transpose_bootstrap_matrices(filenames, length_before_dump=None, save_prefix="",
                                 dump_prefix="",
                                axesnames=None):
    """
    filenames: files of pickled matrices, all with the same shape. some are cause, some are effect

    length_before_dump: number to load before dumping and freeing the matrices

    save_prefix: prefix to save these in

    axesnames: structure such that axesnames[axis][index] gives the name you want associated with that

    return: final_filenames_matr, where[i][j][k] gets you bootstrap filename for that coefficient type in the originalmatrix
    Note this is the same shape as each of the bootstrapped matrices. That coefficient type will just be a raw array.

    A matrix would probably be the best...

    """

    # the number of repeats
    b = len(filenames)

    if length_before_dump == None:
        length_before_dump = b

    num_dump = int(math.ceil(b * 1.0 / length_before_dump))

    # for dumping the arrays


    first_matr = pickle.load(open(filenames[0], 'rU'))


    dump_filenames_matr = np.empty(first_matr.shape + (num_dump,), dtype=object)
    for index, _ in np.ndenumerate(dump_filenames_matr):
        dump_filenames_matr[index] = os.path.join(dump_prefix + "_-dump-index." + "-".join([str(x) for x in index]) + "_" + ".p")


    final_filenames_matr = np.empty(first_matr.shape, dtype=object)
    if axesnames == None:
        for index, _ in np.ndenumerate(final_filenames_matr):
            final_filenames_matr[index] = save_prefix + "_index." + "-".join([str(x) for x in index]) + "_" + ".p"
    else:
        for index, _ in np.ndenumerate(final_filenames_matr):
            final_filenames_matr[index] = save_prefix + "_index." + "-".join([str(x) for x in index]) + "_" + \
            "-".join([axesnames[i][x] for i, x in enumerate(index)]) + "_.p"



    for d in range(num_dump):
        print "Num dump: ", d
        matrs = []
        for i in range(length_before_dump):

            num = d * length_before_dump + i

            if num >= b:
                break


            matrs.append(pickle.load(open(filenames[num], 'rU')))

            if num % 10 == 0:
                print num

        matr = np.stack(matrs, axis=-1)

        # iterate outside it
        for index, _ in np.ndenumerate(first_matr):
            pickle.dump(matr[index], open(dump_filenames_matr[index + (d,)], 'w'))

        del matrs
        del matr


    # we've dumped all the coefs

    # now just merge each of these

    print "merging"

    t = time.time()
    for index, _ in np.ndenumerate(first_matr):
        if index[0] % 10 == 0 and index[1] % 100 == 0 and index[-1] % 100 == 0:
            print "Merging index ", index
            print time.time() - t
            t = time.time()
        arrs = []

        for d in range(num_dump):
            arrs.append(pickle.load(open(dump_filenames_matr[index + (d,)], 'rB')))

        arr = np.concatenate(arrs)

        pickle.dump(arr, open(final_filenames_matr[index], 'w'))

        del arrs
        del arr


    return final_filenames_matr



def compute_bootstrap_stats_matr(bootstrap_coef_file_matr):
    """
    bootstrap_coef_file_matr: computes the coefficients for each bootstrap matrix

    Iteratively compute statistics across each coefficient distribution


    return: stats_matr_dict
    A dictionary where entries are matrices of statistics on the bootstrap coefficient
    """



    stats_matr_dict = collections.OrderedDict()
    stats_matr_dict["mean"] = np.zeros(bootstrap_coef_file_matr.shape)
    stats_matr_dict["std"] = np.zeros(bootstrap_coef_file_matr.shape)
    stats_matr_dict["freq"] = np.zeros(bootstrap_coef_file_matr.shape)
    stats_matr_dict["median"] = np.zeros(bootstrap_coef_file_matr.shape)
    stats_matr_dict["1%"] = np.zeros(bootstrap_coef_file_matr.shape)
    stats_matr_dict["2.5%"] = np.zeros(bootstrap_coef_file_matr.shape)
    stats_matr_dict["5%"] = np.zeros(bootstrap_coef_file_matr.shape)
    stats_matr_dict["95%"] = np.zeros(bootstrap_coef_file_matr.shape)
    stats_matr_dict["97.5%"] = np.zeros(bootstrap_coef_file_matr.shape)
    stats_matr_dict["99%"] = np.zeros(bootstrap_coef_file_matr.shape)

    t = time.time()
    for index, filename in np.ndenumerate(bootstrap_coef_file_matr):
        if index[0] % 10 ==0 and index[1] % 100 == 0 and index[-1] % 100 == 0:
            print "Merging index ", index
            print time.time() - t
            t = time.time()
        coefs = pickle.load(open(filename, 'rU'))

        stats_matr_dict["mean"][index] = np.mean(coefs)
        stats_matr_dict["std"][index] = np.std(coefs, ddof=1)
        stats_matr_dict["freq"][index] = np.mean(coefs != 0)
        stats_matr_dict["median"][index] = np.median(coefs)
        stats_matr_dict["1%"][index] = np.percentile(coefs, 1, interpolation="nearest")
        stats_matr_dict["2.5%"][index] = np.percentile(coefs, 2.5, interpolation="nearest")
        stats_matr_dict["5%"][index] = np.percentile(coefs, 5, interpolation="nearest")
        stats_matr_dict["95%"][index] = np.percentile(coefs, 95, interpolation="nearest")
        stats_matr_dict["97.5%"][index] = np.percentile(coefs, 97.5, interpolation="nearest")
        stats_matr_dict["99%"][index] = np.percentile(coefs, 99, interpolation="nearest")

        del coefs

    return stats_matr_dict




def get_parser():
    # Parse arguments
    import argparse

    description = 'Get results from bootstraped coefficient matrices'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-b', '--bootstrap_file_with_names', required=True,
                        help="Name of file with each row as a bootstrap file " +
                             "of coefficients")

    parser.add_argument('-df', '--data_file', required=True)
    # only aded for genes
    #
    parser.add_argument('-lr', '--load_reps', required=True, type=int)
    #
    parser.add_argument('-o', '--save_prefix', required=True)

    parser.add_argument('-osf', '--outer_save_folder', required=True)

    parser.add_argument('-rsf', '--result_save_folder', required=True)
    #

    parser.add_argument('-l', '--lag',  type=int, required=True)
    #
    parser.add_argument('-tn', '--test_name', required=True)

    parser.add_argument('-da', '--do_align', default=1, type=int)

    parser.add_argument('-dl', '--do_lite', default=0, type=int)

    parser.add_argument('-dr', '--dump_raw', default=1, type=int,
                        help="Save the raw coefficients-- BEFORE alignment.")

    parser.add_argument('-lbd', '--length_before_dump',
                        default=None, type=int)

    parser.add_argument('-tbf', '--transpose_bootstrap_folder', default=None,
                        type=str)

    parser.add_argument('-uabrd', '--unalign_before_raw_dump', type=int,
                        default=0)



    return parser

def load_and_run(args):


    lag = args.lag
    save_prefix = args.save_prefix
    full_save_prefix = os.path.join(args.result_save_folder, save_prefix)



     # Load data file and prepare a file to pass to plotters
    if args.load_reps:
        # load
        genes, _ = gtm.load_basic_rep_file_list(args.data_file)
        # _, genes, _, _, _, _  = gtm.load_rep_file_list(args.data_file)

        # dfs, genes, geneTS, df, timekeys, num_per_keys  = gtm.load_rep_file_list(args.data_file)

        # print "Timekeys: ", timekeys
        # print "Num per key: ", num_per_keys


    else:
        df = pd.read_csv(args.data_file, sep="\t")
        genes, _ = gtm.get_gene_TS(df)
        # dfr = pd.read_csv(args.rand_data_file, sep="\t")
        # genesr, geneTSr = gtm.get_gene_TS(dfr)
        #
        # timekeys = df.columns.values[1:]
        # print "Timekeys: ", timekeys
        #
        # # Num. replicates per key
        # num_per_keys = None

    # load the coef matrices

    # bootstrap_lag_to_matrs = dict([(i,[]) for i in range(1, lag+1)])
    # matr_shape = None
    # with open(args.bootstrap_file_with_names, 'rU') as f:
    #     # check the same shape
    #
    #     # the alignment causes them to be lag x n x n
    #     for line in f.readlines():
    #
    #         print "Loading ", line
    #
    #         # if args.bootstrap_type == "p":
    #         bootstrap_matr = pickle.load(open(line.splitlines()[0], 'rB'))\
    #
    #         # else:
    #         #     raise ValueError("Wrong bootstrap file type specified")
    #
    #
    #         if matr_shape == None:
    #             matr_shape = bootstrap_matr.shape
    #         else:
    #             assert matr_shape == bootstrap_matr.shape
    #
    #
    #         aligned_bootstrap_matr = lc.align_coefs(bootstrap_matr, lag)
    #         for i in range(1, lag+1):
    #             bootstrap_lag_to_matrs[i].append(aligned_bootstrap_matr[i-1])
    #
    # print "Number of bootstrap matrices: ", len(bootstrap_lag_to_matrs[1])





    with open(args.bootstrap_file_with_names, 'rU') as f:
        filenames = [line.split("\n")[0] for line in f.readlines()]




        if args.do_lite:

            stats_matr_dict = cp.bootstrap_matrices_iter_free(filenames)


        else:

            if args.transpose_bootstrap_folder == None:
                raise ValueError("If doing bootstrap calculation, transpose is required")

        # allow the other problem




            transpose_bootstrap_folder = os.path.join(args.outer_save_folder, args.transpose_bootstrap_folder)

            if not os.path.exists(transpose_bootstrap_folder):
                os.makedirs(transpose_bootstrap_folder)
            if not os.path.exists(os.path.join(transpose_bootstrap_folder, "dump-" + str(args.length_before_dump))):
                os.makedirs(os.path.join(transpose_bootstrap_folder, "dump-" + str(args.length_before_dump)))


            transpose_prefix = os.path.join(transpose_bootstrap_folder,
                                            save_prefix)
            dump_prefix = os.path.join(transpose_bootstrap_folder, "dump-" + str(args.length_before_dump), save_prefix)


            t = time.time()
            bootstrap_coef_file_matr = transpose_bootstrap_matrices(filenames,
                                                                    length_before_dump=args.length_before_dump,
                                                                    save_prefix=transpose_prefix,
                                                                    dump_prefix=dump_prefix
                                                                    )
            print "Time to transpose: ", time.time() - t

            bootstrap_coef_filename = dump_prefix + "-NAMES.p"

            pickle.dump(bootstrap_coef_file_matr, open(bootstrap_coef_filename, 'w'))

            print "Bootstrap coef matrix dumped to ", bootstrap_coef_filename

            t = time.time()
            stats_matr_dict = compute_bootstrap_stats_matr(bootstrap_coef_file_matr)
            print "Time to get stats: ", time.time() - t



        # align results

    if args.dump_raw:
        dump_stats_matr_dict = stats_matr_dict.copy()

        if args.unalign_before_raw_dump:
            for k in dump_stats_matr_dict:
                dump_stats_matr_dict[k] = lc.unalign_coefs(dump_stats_matr_dict[k],
                                                           lag)


        for k in dump_stats_matr_dict:
            outfile = full_save_prefix + "_raw_" + k + "_coefs.p"
            with open(outfile, 'w') as f:
                pickle.dump(dump_stats_matr_dict[k], f)

            print "For ", k , "Saved to ", outfile





    if args.do_align:
        for k in stats_matr_dict:
            stats_matr_dict[k] = lc.align_coefs(stats_matr_dict[k], lag)




    # Save the gene matrices

    # Note bootstrap_matr is of form lag x n x n

    full_nets = []
    for i in range(1, lag + 1):
        print "Lag: ", i

        print "Aggregating results"
        #bootstrap_mean, bootstrap_std, bootstrap_freq = cp.get_bootstrap_results(bootstrap_lag_to_matrs[i])


        extra_dict = collections.OrderedDict()
        extra_dict["Test"] = args.test_name
        extra_dict["Lag"] = lag
        extra_dict["Coef"] = i



        nets = []

        for k in stats_matr_dict:
            raw_matr = stats_matr_dict[k][i-1]
            matr_filename = full_save_prefix + "-" + str(i) + "-bootstrap-" + k + "-matrix.txt"

            matr = gtm.save_gene_matrix(matr_filename, matrix=raw_matr, genes=genes)

            print "Saved ", k, " to ", matr_filename

            if k == "mean":
                net = nh.matr_to_net(matr, make_type=False, edge_name="Bootstrap:" + k.capitalize(),
                                      abs_name="AbsBootstrap:" + k.capitalize(),
                                     do_sort=False, extra_dict=extra_dict)
            else:
                net = nh.matr_to_net(matr, make_type=False, edge_name="Bootstrap:" + k.capitalize(),
                                      no_abs=True,
                                     do_sort=False, extra_dict=extra_dict)

            nets.append(net)

        full_net = nets[0]

        for j in range(1, len(nets)):
            full_net = full_net.merge(nets[j], how='outer')



        print "Final net: ", full_net.shape[0]

        sortby = "Bootstrap:Freq"
        print "Sorting by :", sortby
        full_net.sort_values(sortby, inplace=True, ascending=False)

        full_net_filename = full_save_prefix +"-" + str(i) + "-bootstrap-network.txt"
        full_net.to_csv(full_net_filename, sep="\t", index=False)
        print "Written to ", full_net_filename

        full_nets.append(full_net)

    union_net_filename = full_save_prefix + "-union-bootstrap-network.txt"

    if lag > 1:

        m_net = cp.get_max_network(full_nets, max_col="AbsBootstrap:Mean", index_col="Cause-Effect")
        union_net = cp.get_union_network(full_nets + [m_net], suffixes=[str(i) for i in range(1, lag + 1)] + [""])
        print "Max network edges: ", m_net.shape
        print "Union network edges: ", union_net.shape
    else:
        union_net = full_nets[0]

    sortby = "Bootstrap:Freq"
    print "Sorting by :", sortby
    union_net.sort_values(sortby, inplace=True, ascending=False)

    union_net.to_csv(union_net_filename, sep="\t", index=False)
    print "Unioned bootstrap network written to ", union_net_filename





def main():
    load_and_run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()