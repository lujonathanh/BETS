__author__ = 'jlu96'

import causal_pipeline as cp
import sys
import pickle
import pandas as pd
import geneTSmunging as gtm


def get_parser():
    # Parse arguments
    import argparse

    description = 'Run cross-validate using input datasets'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-d', '--data_file', required=True)

    parser.add_argument('-lr', '--load_reps', required=True, type=int)

    parser.add_argument('-o', '--output_names', nargs='+', required=True, help="Place where CV results are saved")

    parser.add_argument('-hl', '--hyper_file', required=True)

    parser.add_argument('-t', '--test', required=True)

    parser.add_argument('-l', '--lag',  type=int, required=True)

    parser.add_argument('-rl', '--row_file', default=None)


    return parser



def load_and_run(args):

    data_file = args.data_file
    output_names = args.output_names


    assert args.test in {'e', 'l', 'r'}
    fit_method = cp.test2fit_method[args.test]

    lag = args.lag
    hyperlist = pickle.load(open(args.hyper_file, 'rB'))

    if args.row_file != None:
        rows = pickle.load(open(args.row_file, 'rB'))
    else:
        rows = None


    # Load data file
    if args.load_reps:

        genes, geneTS = gtm.load_basic_rep_file_list(data_file)
        #dfs, genes, geneTS, df, __, __  = gtm.load_rep_file_list(data_file)
    else:
        df = pd.read_csv(data_file, sep="\t")
        genes, geneTS = gtm.get_gene_TS(df)





    best_hyper, best, hyper_df, hyper_fit_dfs = cp.run_cross_validate(geneTS, fit_method=fit_method,
                                                                      hyperlist=hyperlist, lag=lag, rows=rows,
                                                                      has_reps=args.load_reps)



    print "Best hyper is : ", best_hyper
    print "Best result : ", best

    print "Hyper df: "
    print hyper_df

    for output_name, hyper_fit_df, hyper in zip(output_names, hyper_fit_dfs, hyperlist):
        hyper_fit_df.to_csv(output_name, sep="\t", index=0)

        print "Result for ", hyper, " written to ", output_name



def main():
    load_and_run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()