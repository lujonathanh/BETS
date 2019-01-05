__author__ = 'jlu96'


import sys
import pickle
import pandas as pd
import prep_jobs as pj

def get_parser():
    # Parse arguments
    import argparse

    description = 'Given the baseline, per gene hyperparameter fit results, choose the best hyperparameter'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-hfd', '--hyper_fit_names_to_int', required=True, help="Data file with columns containing the hyperparameters")

    parser.add_argument('-ind', '--int_name_dfname', required=True)

    parser.add_argument('-hl', '--hyper_file', required=True, help="Pickle of the original list. Make sure this has same order of the row names in integration")


    return parser


def load_and_run(args):


    hyperlist = pickle.load(open(args.hyper_file, 'rB'))

    # Integrate all hyperparameter files into dfs
    output_df = pd.read_csv(args.hyper_fit_names_to_int, sep="\t")
    int_name_df = pd.read_csv(args.int_name_dfname, sep="\t")

    print "Integrated into: "
    print int_name_df.head()

    assert set(output_df.columns.values) == set(int_name_df.columns.values)

    hyper_fit_dfs = []

    for x, hyper in zip(output_df, hyperlist):
        filenames = output_df[x].values

        print
        print x
        print "hyper is ", hyper
        print "Files to integrate are ", filenames
        integrated_filename = int_name_df[x].values[0]
        print "Integrated_file is ", integrated_filename

        hyper_fit_df = pj.integrate_dfs(filenames, integrated_filename=integrated_filename)

        hyper_fit_dfs.append(hyper_fit_df)



def main():
    load_and_run(get_parser().parse_args(sys.argv[1:]))




if __name__ == '__main__':
    main()