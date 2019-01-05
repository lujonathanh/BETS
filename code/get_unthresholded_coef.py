__author__ = 'jlu96'


import numpy as np
import pandas as pd
import pickle
import sys
import os
import collections
import geneTSmunging as gtm
import fdr_control as fc
import network_helpers as nh
import causal_pipeline as cp
import lag_conversion as lc




def get_parser():
    # Parse arguments
    import argparse

    description = 'Get a network file with all coefficients, unthresholded'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-df', '--data_file', required=True)

    parser.add_argument('-lr', '--load_reps', required=True, type=int)

    parser.add_argument('-o', '--save_prefix', required=True)

    parser.add_argument('-cf', '--coef_file', required=True)

    parser.add_argument('-l', '--lag',  type=int, required=True)

    parser.add_argument('-tn', '--test_name', required=True)

    return parser

def load_and_run(args):

    lag = args.lag
    save_prefix = args.save_prefix


     # Load data file and prepare a file to pass to plotters
    if args.load_reps:
        # load
        dfs, genes, geneTS, df, timekeys, num_per_keys  = gtm.load_rep_file_list(args.data_file)
        dfsr, genesr, geneTSr, dfr, __, __  = gtm.load_rep_file_list(args.rand_data_file)

        # get shared prefix timekeys

        print "Timekeys: ", timekeys
        print "Num per key: ", num_per_keys


    else:
        df = pd.read_csv(args.data_file, sep="\t")
        genes, geneTS = gtm.get_gene_TS(df)
        dfr = pd.read_csv(args.rand_data_file, sep="\t")
        genesr, geneTSr = gtm.get_gene_TS(dfr)

        timekeys = df.columns.values[1:]
        print "Timekeys: ", timekeys

        # Num. replicates per key
        num_per_keys = None


    assert (geneTS.shape == geneTSr.shape)
    assert (genes == genesr).all()




    coefs = pickle.load(open(args.coef_file, 'rB'))



    # Align the coefs


    # print "Aligning coefficients"
    acoefs = lc.align_coefs(coefs, lag)


    print "Removing alphas (gene-on-self effects) "

    acoefs = lc.remove_alphas(acoefs, lag)




    coef_nets = []
    coefr_nets = []


    # Save the gene matrices
    for i in range(acoefs.shape[0]):
        coef_matr_filename = save_prefix + "-" + str(i+1) + "-matrix.txt"
        coefr_matr_filename = save_prefix + "-" + str(i+1) + "-r-matrix.txt"

        coef_net_filename = save_prefix + "-" + str(i+1) + "-network.txt"
        coefr_net_filename = save_prefix + "-" + str(i+1) + "-r-network.txt"



        coef_matr = gtm.save_gene_matrix(filename= coef_matr_filename,matrix=acoefs[i], genes=genes)


        extra_dict = collections.OrderedDict()
        extra_dict["Test"] = args.test_name
        extra_dict["Lag"] = acoefs.shape[0]
        extra_dict["Coef"] = i + 1


        coef_net = nh.matr_to_net(coef_matr, extra_dict=extra_dict, make_type=False)

        coef_net.to_csv(coef_net_filename, sep="\t", index=False)

        coef_nets.append(coef_net)


        print "Coef ", i+1
        print "Networks written to "
        print coef_net_filename
        print coefr_net_filename


    # max_net_filename = save_prefix + "-max-network.txt"
    # max_r_net_filename = save_prefix + "-max-r-network.txt"
    union_net_filename = save_prefix + "-union-network.txt"
    union_r_net_filename = save_prefix + "-union-r-network.txt"

    if acoefs.shape[0] > 1:
        m_net = cp.get_max_network(coef_nets, max_col="AbsWeight", index_col="Cause-Effect")
        union_net = cp.get_union_network(coef_nets + [m_net], suffixes=[str(i) for i in range(1, acoefs.shape[0] + 1)] + [""])
        print "Max network edges: ", m_net.shape
        print "Union network edges: ", union_net.shape
    else:
        union_net = coef_nets[0]
    union_net.to_csv(union_net_filename, sep="\t", index=False)


    # print "Max networks written to "
    # print max_net_filename
    # print max_r_net_filename
    print "Unioned networks written to "
    print union_net_filename
    print union_r_net_filename







    if not os.path.exists("plots"):
        os.makedirs("plots")
    if not os.path.exists("plots" + os.sep + "betas"):
        os.makedirs("plots" + os.sep + "betas")

    # Plot the betas
    for i in range(acoefs.shape[0]):

        if len(np.nonzero(acoefs[i])[0]) > 0 and len(np.nonzero(acoefsr[i])[0]) > 0:


            fc.plot_betas(acoefs[i][np.nonzero(acoefs[i])].flatten(), acoefsr[i][np.nonzero(acoefsr[i])].flatten(), filename="plots" + os.sep + "betas" + os.sep + "beta_nonzero_coef-" + str(i+1), title="Causal coefs, Coef " + str(i+1), xlabel="Causal Coefficient")
            fc.plot_betas(acoefs[i][np.nonzero(acoefs[i])].flatten(), acoefsr[i][np.nonzero(acoefsr[i])].flatten(), filename="plots" + os.sep + "betas" + os.sep + "beta_nonzero_coef-" + str(i+1) + "_zoom-in-90", zoom_in_top_percentile=95, zoom_in_bottom_percentile=5, title="Causal coefs, Coef " + str(i+1), xlabel="Causal Coefficient")


            fc.plot_betas(np.absolute(acoefs[i][np.nonzero(acoefs[i])].flatten()), np.absolute(acoefsr[i][np.nonzero(acoefsr[i])].flatten()), filename="plots" + os.sep + "betas" + os.sep + "beta_abs_coef-" + str(i+1), title="Absolute causal coefs, Coef " + str(i+1), xlabel="Absolute Causal Coefficient")
            fc.plot_betas(np.absolute(acoefs[i][np.nonzero(acoefs[i])].flatten()), np.absolute(acoefsr[i][np.nonzero(acoefsr[i])].flatten()), filename="plots" + os.sep + "betas" + os.sep + "beta_abs_coef-" + str(i+1) + "_zoom-in-bottom-95", zoom_in_top_percentile=95, title="Absolute causal coefs, Coef " + str(i+1), xlabel="Absolute Causal Coefficient")
            fc.plot_betas(np.absolute(acoefs[i][np.nonzero(acoefs[i])].flatten()), np.absolute(acoefsr[i][np.nonzero(acoefsr[i])].flatten()), filename="plots" + os.sep + "betas" + os.sep + "beta_abs_coef-" + str(i+1) + "_zoom-in-top-5", zoom_in_bottom_percentile=95, title="Absolute causal coefs, Coef " + str(i+1), xlabel="Absolute Causal Coefficient")


        print "Coef ", i+1
        print "Plots of betas written to: plots" + os.sep + "betas"





    # get FDRS
    fdrs = [0.01, 0.05, 0.1, 0.2]


    acoefs_fdrs = []
    sf_dfs = []

    for fdr in fdrs:

        fdr_dir = "fdr-" + str(fdr) + "-" + stratify_by
        if not os.path.exists(fdr_dir):
            os.makedirs(fdr_dir)

        fdr_prefix = fdr_dir + os.sep + save_prefix

        acoefs_fdr = np.zeros(acoefs.shape)



        fdr_nets = []

        print "*************"
        for i in range(acoefs.shape[0]):
            print "-----"
            print "FDR = ", fdr
            print "Lag ", lag
            print "Coef ", i + 1
            print "Stratify ", stratify_by
            acoefs_fdr[i], threshes = fc.get_abs_thresh(acoefs[i], acoefsr[i], fdr, stratify_by=stratify_by)
            # print "Threshes", threshes



            fdr_matr_filename = fdr_prefix + "-" + str(i+1) + "-fdr-" + str(fdr) + "-" + stratify_by + "-matrix.txt"
            fdr_net_filename = fdr_prefix + "-" + str(i+1) + "-fdr-" + str(fdr) + "-" + stratify_by + "-network.txt"

            fdr_matr = gtm.save_gene_matrix(fdr_matr_filename,  matrix=acoefs_fdr[i], genes=genes)
            pickle.dump(threshes, open(fdr_prefix + "-" + str(i+1) + "-fdr-" + str(fdr) + "-" + stratify_by  + "-threshes.p", 'wB'))


            extra_dict = collections.OrderedDict()
            extra_dict["Test"] = args.test_name
            extra_dict["Lag"] = acoefs.shape[0]
            extra_dict["Coef"] = i + 1


            fdr_net = nh.matr_to_net(fdr_matr, extra_dict=extra_dict, make_type=False)
            fdr_net.to_csv(fdr_net_filename, sep="\t", index=False)
            fdr_nets.append(fdr_net)

            # write summary readme
            sf_df = fc.summarize_fdr(matr=acoefs_fdr[i],  test=args.test_name,
                                     fdr=fdr, lag=lag, coef=i+1, hyper=best_hyper, thresh=threshes,
                                    readme_name=fdr_prefix + "-" + str(i+1) + "-fdr-" + str(fdr) + "-" + stratify_by + "-README.txt",
                                    matrixname=fdr_matr_filename,
                                    filename=fdr_net_filename)

            sf_dfs.append(sf_df)

            print "Network edges: ", fdr_net.shape[0]


        if acoefs_fdr.shape[0] > 1:
            m_net = cp.get_max_network(fdr_nets, max_col="AbsWeight", index_col="Cause-Effect")
            union_net = cp.get_union_network(fdr_nets + [m_net], suffixes=[str(i) for i in range(1, acoefs_fdr.shape[0] + 1)] + [""])

        else:
            union_net = fdr_nets[0]

        union_net_filename = fdr_prefix + "-union-fdr-" + str(fdr) + "-" + stratify_by +  "-network.txt"
        union_net.to_csv(union_net_filename, sep="\t", index=False)

        print "Union network edges", union_net.shape[0]
        print "Written to ", union_net_filename


        acoefs_fdrs.append(acoefs_fdr.copy())


    all_sf_dfs = pd.concat(sf_dfs)

    all_sf_dfs.to_csv("fit_all_summary_fdr-" + stratify_by + ".txt", sep="\t", index=False)
    print "********"
    print "Summaries of all fdrs written to fit_all_summary_fdr-" + stratify_by + ".txt"
    print "Matrices done."



    with open("matrices_done.txt", 'w') as donefile:
        donefile.write("done\n")



def main():
    load_and_run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()