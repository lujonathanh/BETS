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

    description = 'Run causal pipeline on input datasets'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-df', '--data_file', required=True)

    parser.add_argument('-rdf', '--rand_data_file', required=True)

    parser.add_argument('-lr', '--load_reps', required=True, type=int)

    parser.add_argument('-o', '--save_prefix', required=True)

    parser.add_argument('-of', '--output_folder', default=None)

    parser.add_argument('-bh', '--best_hyper_file', default=None)

    parser.add_argument('-cf', '--coef_file', required=True)

    parser.add_argument('-if', '--intercept_file', required=True)

    parser.add_argument('-cfr', '--coef_rand_file', required=True)

    # parser.add_argument('-ifr', '--intercept_rand_file', required=True)

    parser.add_argument('-fr', '--fit_result_file', required=True)

    parser.add_argument('-frr', '--fit_result_rand_file', required=True)

    parser.add_argument('-l', '--lag',  type=int, required=True)

    parser.add_argument('-tn', '--test_name', required=True)

    parser.add_argument('-sb', '--stratify_by', type=str, required=True)

    # parser.add_argument('-pcf', '--plot_coef_fdr', type=int, default=0)

    parser.add_argument('-pa', '--plot_all', type=int, default=0,
                        help="Deprecated, don't worry about this. It's not used.")

    parser.add_argument('-pc', '--plot_coef', type=int, default=1,
                        help="Plot the coefficient distribution")


    return parser

def load_and_run(args):

    lag = args.lag
    save_prefix = args.save_prefix

    assert args.stratify_by in {"e", "n"}

    stratify_by = cp.args2stratify_by[args.stratify_by]

    if args.output_folder == None:
        args.output_folder = "."


     # Load data file and prepare a file to pass to plotters
    if args.load_reps:
        # load

        genes, geneTS = gtm.load_basic_rep_file_list(args.data_file)
        genesr, geneTSr = gtm.load_basic_rep_file_list(args.rand_data_file)



        # dfs, genes, geneTS, df, timekeys, num_per_keys  = gtm.load_rep_file_list(args.data_file)
        # dfsr, genesr, geneTSr, dfr, __, __  = gtm.load_rep_file_list(args.rand_data_file)

        # get shared prefix timekeys

        # print "Timekeys: ", timekeys
        # print "Num per key: ", num_per_keys


    else:
        df = pd.read_csv(args.data_file, sep="\t")
        genes, geneTS = gtm.get_gene_TS(df)
        dfr = pd.read_csv(args.rand_data_file, sep="\t")
        genesr, geneTSr = gtm.get_gene_TS(dfr)

        timekeys = df.columns.values[1:]
        print("Timekeys: ", timekeys)

        # Num. replicates per key
        num_per_keys = None


    assert (geneTS.shape == geneTSr.shape)
    assert (genes == genesr).all()




    coefs = pickle.load(open(args.coef_file, 'rb'))
    intercepts = pickle.load(open(args.intercept_file, 'rb'))
    fit_result_df = pd.read_csv(args.fit_result_file, sep="\t")

    coefsr = pickle.load(open(args.coef_rand_file, 'rb'))
    # interceptsr = pickle.load(open(args.intercept_rand_file, 'rb'))
    fit_result_dfr = pd.read_csv(args.fit_result_rand_file, sep="\t")

    if args.best_hyper_file != None:
        best_hyper = pickle.load(open(args.best_hyper_file, 'rb'))
    else:
        best_hyper = None



    print("RESULTS")
    
    
    print("*************************")
    print("RESIDUALS: ")
    
    print("*************************")
    print("NORMAL: ")
    cp.summarize_fit(coefs, intercepts, fit_result_df, filename= os.path.join(args.output_folder, "fit_all_summary_normal.txt"), hyper=best_hyper,
                     test_name=args.test_name, lag=lag)



    # Align the coefs


    # print "Aligning coefficients"
    acoefs = lc.align_coefs(coefs, lag)
    acoefsr = lc.align_coefs(coefsr, lag)


    print("Removing alphas (gene-on-self effects) ")

    acoefs = lc.remove_alphas(acoefs, lag)
    acoefsr = lc.remove_alphas(acoefsr, lag)




    coef_nets = []
    coefr_nets = []


    # Save the gene matrices
    for i in range(acoefs.shape[0]):
        coef_matr_filename = os.path.join(args.output_folder, save_prefix + "-" + str(i+1) + "-matrix.txt")
        coefr_matr_filename = os.path.join(args.output_folder, save_prefix + "-" + str(i+1) + "-r-matrix.txt")

        coef_net_filename = os.path.join(args.output_folder, save_prefix + "-" + str(i+1) + "-network.txt")
        coefr_net_filename = os.path.join(args.output_folder, save_prefix + "-" + str(i+1) + "-r-network.txt")



        coef_matr = gtm.save_gene_matrix(filename= coef_matr_filename,matrix=acoefs[i], genes=genes)
        coefr_matr = gtm.save_gene_matrix(filename=coefr_matr_filename, matrix=acoefsr[i], genes=genes)


        extra_dict = collections.OrderedDict()
        extra_dict["Test"] = args.test_name
        extra_dict["Lag"] = acoefs.shape[0]
        extra_dict["Coef"] = i + 1


        coef_net = nh.matr_to_net(coef_matr, extra_dict=extra_dict, make_type=False)
        coefr_net = nh.matr_to_net(coefr_matr, extra_dict=extra_dict, make_type=False)

        coef_net.to_csv(coef_net_filename, sep="\t", index=False)
        coefr_net.to_csv(coefr_net_filename, sep="\t", index=False)

        coef_nets.append(coef_net)
        coefr_nets.append(coefr_net)


        print("Coef ", i+1)
        print("Networks written to ")
        print(coef_net_filename)
        print(coefr_net_filename)


    # max_net_filename = save_prefix + "-max-network.txt"
    # max_r_net_filename = save_prefix + "-max-r-network.txt"
    union_net_filename = os.path.join(args.output_folder, save_prefix + "-union-network.txt")
    union_r_net_filename = os.path.join(args.output_folder, save_prefix + "-union-r-network.txt")

    if acoefs.shape[0] > 1:
        m_net = cp.get_max_network(coef_nets, max_col="AbsWeight", index_col="Cause-Effect")
        union_net = cp.get_union_network(coef_nets + [m_net], suffixes=[str(i) for i in range(1, acoefs.shape[0] + 1)] + [""])
        print("Max network edges: ", m_net.shape)
        print("Union network edges: ", union_net.shape)
    else:
        union_net = coef_nets[0]
    union_net.to_csv(union_net_filename, sep="\t", index=False)

    if acoefsr.shape[0] > 1:
        m_net = cp.get_max_network(coefr_nets, max_col="AbsWeight", index_col="Cause-Effect")
        union_r_net = cp.get_union_network(coefr_nets + [m_net], suffixes=[str(i) for i in range(1, acoefs.shape[0] + 1)] + [""])
    else:
        union_r_net = coefr_nets[0]
    union_r_net.to_csv(union_r_net_filename, sep="\t", index=False)


    # print "Max networks written to "
    # print max_net_filename
    # print max_r_net_filename
    print("Unioned networks written to ")
    print(union_net_filename)
    print(union_r_net_filename)







    if not os.path.exists(os.path.join(args.output_folder,"plots")):
        os.makedirs(os.path.join(args.output_folder,"plots"))

    if args.plot_coef:
        if not os.path.exists(os.path.join(args.output_folder, "plots", "betas")):
            os.makedirs(os.path.join(args.output_folder, "plots", "betas"))

        # Plot the betas
        for i in range(acoefs.shape[0]):

            if len(np.nonzero(acoefs[i])[0]) > 0 and len(np.nonzero(acoefsr[i])[0]) > 0:


                fc.plot_betas(acoefs[i][np.nonzero(acoefs[i])].flatten(), acoefsr[i][np.nonzero(acoefsr[i])].flatten(), filename= os.path.join(args.output_folder, "plots", "betas", "beta_nonzero_coef-" + str(i+1)), title="Causal coefs, Coef " + str(i+1), xlabel="Causal Coefficient")
                fc.plot_betas(acoefs[i][np.nonzero(acoefs[i])].flatten(), acoefsr[i][np.nonzero(acoefsr[i])].flatten(), filename= os.path.join(args.output_folder, "plots", "betas", "beta_nonzero_coef-" + str(i+1) + "_zoom-in-90"), zoom_in_top_percentile=95, zoom_in_bottom_percentile=5, title="Causal coefs, Coef " + str(i+1), xlabel="Causal Coefficient")


                fc.plot_betas(np.absolute(acoefs[i][np.nonzero(acoefs[i])].flatten()), np.absolute(acoefsr[i][np.nonzero(acoefsr[i])].flatten()), filename= os.path.join(args.output_folder, "plots", "betas", "beta_abs_coef-" + str(i+1)), title="Absolute causal coefs, Coef " + str(i+1), xlabel="Absolute Causal Coefficient")
                fc.plot_betas(np.absolute(acoefs[i][np.nonzero(acoefs[i])].flatten()), np.absolute(acoefsr[i][np.nonzero(acoefsr[i])].flatten()), filename=  os.path.join(args.output_folder, "plots", "betas", "beta_abs_coef-" + str(i+1) + "_zoom-in-bottom-95"), zoom_in_top_percentile=95, title="Absolute causal coefs, Coef " + str(i+1), xlabel="Absolute Causal Coefficient")
                fc.plot_betas(np.absolute(acoefs[i][np.nonzero(acoefs[i])].flatten()), np.absolute(acoefsr[i][np.nonzero(acoefsr[i])].flatten()), filename= os.path.join(args.output_folder, "plots", "betas", "beta_abs_coef-" + str(i+1) + "_zoom-in-top-5"), zoom_in_bottom_percentile=95, title="Absolute causal coefs, Coef " + str(i+1), xlabel="Absolute Causal Coefficient")


            print("Coef ", i+1)
            print("Plots of betas written to: ", os.path.join(args.output_folder, "plots", "betas"))





    # get FDRS
    fdrs = [0.01, 0.05, 0.1, 0.2]


    acoefs_fdrs = []
    sf_dfs = []

    for fdr in fdrs:

        fdr_dir = os.path.join(args.output_folder, "fdr-" + str(fdr) + "-" + stratify_by)
        if not os.path.exists(fdr_dir):
            os.makedirs(fdr_dir)

        fdr_prefix = fdr_dir + os.sep + save_prefix

        # in case we want there to be an intermediate directory for fdr, like the bootstrap case.
        # if not os.path.exists(os.path.dirname(fdr_prefix)):
        #     os.makedirs(os.path.dirname(fdr_prefix))

        acoefs_fdr = np.zeros(acoefs.shape)



        fdr_nets = []

        print("*************")
        for i in range(acoefs.shape[0]):
            print("-----")
            print("FDR = ", fdr)
            print("Lag ", lag)
            print("Coef ", i + 1)
            print("Stratify ", stratify_by)
            acoefs_fdr[i], threshes = fc.get_abs_thresh(acoefs[i], acoefsr[i], fdr, stratify_by=stratify_by)
            # print "Threshes", threshes



            fdr_matr_filename = fdr_prefix + "-" + str(i+1) + "-fdr-" + str(fdr) + "-" + stratify_by + "-matrix.txt"
            fdr_net_filename = fdr_prefix + "-" + str(i+1) + "-fdr-" + str(fdr) + "-" + stratify_by + "-network.txt"

            fdr_matr = gtm.save_gene_matrix(fdr_matr_filename,  matrix=acoefs_fdr[i], genes=genes)
            pickle.dump(threshes, open(fdr_prefix + "-" + str(i+1) + "-fdr-" + str(fdr) + "-" + stratify_by  + "-threshes.p", 'wb'))


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

            print("Network edges: ", fdr_net.shape[0])


        if acoefs_fdr.shape[0] > 1:
            m_net = cp.get_max_network(fdr_nets, max_col="AbsWeight", index_col="Cause-Effect")
            union_net = cp.get_union_network(fdr_nets + [m_net], suffixes=[str(i) for i in range(1, acoefs_fdr.shape[0] + 1)] + [""])

        else:
            union_net = fdr_nets[0]

        union_net_filename = fdr_prefix + "-union-fdr-" + str(fdr) + "-" + stratify_by +  "-network.txt"
        union_net.to_csv(union_net_filename, sep="\t", index=False)

        print("Union network edges", union_net.shape[0])
        print("Written to ", union_net_filename)



        fdr_agg_matr_filename = fdr_prefix + "-union-fdr-" + str(fdr) + "-" + stratify_by +  "-coefs.p"
        pickle.dump(acoefs_fdr, open(fdr_agg_matr_filename, 'wb'))

        print("Thresholded matrix written as pickle file: ", fdr_agg_matr_filename)

        acoefs_fdrs.append(acoefs_fdr.copy())


    all_sf_dfs = pd.concat(sf_dfs)


    # Hack to allow the base to still be fit_all_summary_fdr-stratby.txt
    # While the bootstrap will write to its own file, in its own corresponding folder
    # bullshit. just sent the output folder


    save_file = os.path.join(args.output_folder ,"fit_all_summary_fdr-" + stratify_by + ".txt")
    all_sf_dfs.to_csv(save_file, sep="\t", index=False)
    print("********")
    print("Summaries of all fdrs written to ", save_file)
    print("Matrices done.")



    with open(os.path.join(args.output_folder, "matrices_done.txt"), 'w') as donefile:
        donefile.write("done\n")




    #
    #
    # if args.plot_coef_fdr:
    #     print "*******"
    #     print "PLOTS"
    #     for i, fdr in zip(range(len(fdrs)), fdrs):
    #         acoefs_fdr = acoefs_fdrs[i]
    #
    #         plot_fdr_folder = os.path.join(args.output_folder, "plots", "fdr-" + str(fdr))
    #
    #         if not os.path.exists(plot_fdr_folder):
    #             os.makedirs(plot_fdr_folder)
    #
    #         # Only plot the bar if replicates were loaded
    #         cp.plot_all_coef(acoefs_fdr, df, genes, lag, file_prefix=os.path.join(plot_fdr_folder, save_prefix+ "-" ),
    #                          plot_bar=args.load_reps, keys=timekeys, num_per_keys=num_per_keys,
    #                            linewidth=2,
    #                            capsize=5,
    #                            capwidth=2, verbose=True)
    #
    #         # Plot them without error bars just to check
    #         if args.load_reps:
    #             cp.plot_all_coef(acoefs_fdr, df, genes, lag, file_prefix=os.path.join(plot_fdr_folder, save_prefix+ "-nobar-"),
    #                             plot_bar=False, keys=timekeys, num_per_keys=num_per_keys,
    #                                linewidth=2,
    #                                capsize=5,
    #                                capwidth=2)
    #
    #
    #
    #         print "FDR plots written to: ", plot_fdr_folder



    # Plot all the coefs
    # NOTE: this will take a long time!
    # I'm declaring the below deprecated, the user will never want to do this -- JL
    # if args.plot_all:
    #
    #     raise ValueError("Fix all the below first before trying to do plot all")
    #
    #     if not os.path.exists("plots" + os.sep + "original"):
    #         os.makedirs("plots" + os.sep + "original")
    #     cp.plot_all_coef(acoefs, df, genes, lag, file_prefix="plots" + os.sep + "original" + os.sep + save_prefix + "-",
    #                         plot_bar=args.load_reps, keys=timekeys, num_per_keys=num_per_keys,
    #                        linewidth=2,
    #                        capsize=5,
    #                        capwidth=2)
    #     print "Original plots written to: ", "plots" + os.sep + "original"
    #
    #     if not os.path.exists("plots" + os.sep + "randomized"):
    #         os.makedirs("plots" + os.sep + "randomized")
    #     cp.plot_all_coef(acoefsr, dfr, genes, lag, file_prefix="plots" + os.sep + "randomized" + os.sep + save_prefix+ "-",
    #                         plot_bar=args.load_reps, keys=timekeys, num_per_keys=num_per_keys,
    #                        linewidth=2,
    #                        capsize=5,
    #                        capwidth=2)
    #
    #     print "Randomized plots written to: ", "plots" + os.sep + "randomized"




def main():
    load_and_run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()