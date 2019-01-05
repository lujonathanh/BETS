__author__ = 'jlu96'


import sys
import CausalTests as ct
import geneTSmunging as gtm
import pickle
import pandas as pd




def get_parser():
    # Parse arguments
    import argparse

    description = 'Apply a pre-specified causal test to an input dataset and randomized dataset where each row is a geene' \
                  'and its tim points, specifying which rows to test as effect,'\
                    'Save the results (and parameters if needed), write output coefficients to a pickle file.' \

    parser = argparse.ArgumentParser(description=description)


    parser.add_argument('-d', '--data_file', required=True)

    parser.add_argument('-rd', '--rand_data_file', required=True)

    parser.add_argument('-a', '--args_file', required=True)

    parser.add_argument('-t', '--test', required=True)

    parser.add_argument('-rl', '--rowlist_file', required=True)

    parser.add_argument('-o', '--output_name', required=True)

    parser.add_argument('-oa', '--output_all_name', default=None)

    parser.add_argument('-ou', '--output_use_name', default=None)


    parser.add_argument('-or', '--output_rand_name', required=True)
    parser.add_argument('-ora', '--output_rand_all_name', default=None)

    parser.add_argument('-oru', '--output_rand_use_name', default=None)

    return parser



def run(args):

    if args.test not in {"p", "e"}:
        raise ValueError("args.test must be p (pairwise granger) or e (elastic net)")

    # load the data

    df = pd.read_csv(args.data_file, sep="\t")

    genes = df['gene'].values

    found_genes, geneTS = gtm.get_gene_TS(df, genes)


    dfr = pd.read_csv(args.rand_data_file, sep="\t")

    genesr = dfr['gene'].values

    found_genesr, geneTSr = gtm.get_gene_TS(dfr, genesr)

    n = geneTSr.shape[0]

    args_dict = ct.load_kwargs_file(argsfile=args.args_file)

    print "Arguments: "
    print args_dict




    if args.rowlist_file != None:
        with open(args.rowlist_file, 'rU') as f:
            rowlist = eval(f.readline())
    else:
        rowlist = range(n)


    if args.test == "e":
        beta_tuple, all_res_df, use_df = ct.enet_granger_causality_row_cv(geneTS, geneTS, rowlist, **args_dict)
        with open(args.output_name, 'w') as outfile:
            pickle.dump(beta_tuple, outfile)
        all_res_df.to_csv(args.output_all_name, sep="\t", index=False)
        use_df.to_csv(args.output_use_name, sep="\t", index=False)


        param_df = use_df[["alpha", "lambda.min", "Row"]]

        print "Orig shape", geneTS.shape
        print "Rand shape ", geneTSr.shape

        rand_beta_tuple, rand_all_res_df, rand_use_df = ct.enet_granger_causality_row_load(geneTSr, geneTS, rowlist, param_df, **args_dict)

        with open(args.output_rand_name, 'w') as outfile:
            pickle.dump(rand_beta_tuple, outfile)

        rand_all_res_df.to_csv(args.output_rand_all_name, sep="\t", index=False)
        rand_use_df.to_csv(args.output_rand_use_name, sep="\t", index=False)

        print "HIIIIIII"
        print "Output written to ", args.output_name
        print "All results written to ", args.output_all_name
        print "Used params written to ", args.output_use_name

        print "Rand output written to ", args.output_rand_name
        print "All rand results written to ", args.output_rand_all_name
        print "Used rand params written to ", args.output_rand_use_name

    elif args.test == "p":
        beta_tuple = ct.pairwise_granger_causality_row(geneTS, geneTS, rowlist, **args_dict)
        with open(args.output_name, 'w') as outfile:
            pickle.dump(beta_tuple, outfile)

        rand_beta_tuple = ct.pairwise_granger_causality_row(geneTSr, geneTS, rowlist, **args_dict)
        with open(args.output_rand_name, 'w') as outfile:
            pickle.dump(rand_beta_tuple, outfile)



def main():
    run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()
