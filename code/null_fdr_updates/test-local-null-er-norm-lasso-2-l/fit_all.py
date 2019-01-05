__author__ = 'jlu96'


import causal_pipeline as cp
import sys
import pickle
import pandas as pd
import geneTSmunging as gtm


def get_parser():
    # Parse arguments
    import argparse

    description = 'Perform all the fits.'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-d', '--data_file', required=True)

    parser.add_argument('-rd', '--rand_data_file', required=True)

    parser.add_argument('-lr', '--load_reps', required=True, type=int)

    parser.add_argument('-o', '--out_prefix', required=True, help="Prefix for saving results")

    parser.add_argument('-bh', '--best_hyper_file', required=True)

    parser.add_argument('-t', '--test', required=True)

    parser.add_argument('-l', '--lag',  type=int, required=True)

    parser.add_argument('-rl', '--row_file', default=None)

    parser.add_argument('-n', '--null', required=True)


    return parser




def load_and_run(args):

    data_file = args.data_file
    rand_data_file = args.rand_data_file
    save_prefix = args.out_prefix


    assert args.test in {'e', 'l', 'r'}
    fit_method = cp.test2fit_method[args.test]

    lag = args.lag
    best_hyper = pickle.load(open(args.best_hyper_file, 'rB'))

    if args.row_file != None:
        rows = pickle.load(open(args.row_file, 'rB'))
    else:
        rows = None

    assert args.null in {"l", "g"}


    # Load data file    # Load data file
    if args.load_reps:
        dfs, genes, geneTS, df, __, __  = gtm.load_rep_file_list(data_file)
        dfsr, genesr, geneTSr, dfr, __, __  = gtm.load_rep_file_list(rand_data_file)

    else:
        df = pd.read_csv(data_file, sep="\t")
        genes, geneTS = gtm.get_gene_TS(df)
        dfr = pd.read_csv(rand_data_file, sep="\t")
        genesr, geneTSr = gtm.get_gene_TS(dfr)


    assert (geneTS.shape == geneTSr.shape)
    assert (genes == genesr).all()


    coefs, intercepts, fit_result_df, coefsr, fit_result_dfr = cp.run(geneTS,
                                                                           geneTSr,
                                                                            hyper=best_hyper,
                                                                            fit_method=fit_method,
                                                                            lag=lag,
                                                                            rows=rows,
                                                                            save_prefix=save_prefix,
                                                                            has_reps=args.load_reps,
                                                                            null=args.null
                                                                           )


    print "RESULTS of causal fit"
    print "*************************"
    print "NORMAL: "
    cp.summarize_fit(coefs, intercepts, fit_result_df)


    # print "*************************"
    # print "RANDOM:"
    # cp.summarize_fit(coefsr, interceptsr, fit_result_dfr)







def main():
    load_and_run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()