__author__ = 'jlu96'


def get_parser():
    # Parse arguments
    import argparse

    description = 'Run causal pipeline on input datasets'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-df', '--data_file', required=True)

    parser.add_argument('-rdf', '--rand_data_file', required=True)

    parser.add_argument('-o', '--out_prefix', required=True)

    parser.add_argument('-h', '--hyper_file', required=True)

    parser.add_argument('-nf', '--num_folds', type=int, required=True)

    parser.add_argument('-t', '--test', required=True)

    parser.add_argument('-l', '--lag',  type=int, required=True)

    parser.add_argument('-sb', '--stratify_by', type=str, required=True)

    return parser