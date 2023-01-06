__author__ = 'jlu96'


__author__ = 'jlu96'

import sys
import pandas as pd
import prep_jobs as pj



def get_parser():
    # Parse arguments
    import argparse

    description = 'Integrate output to a pickle or integrated matrix file.'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-i', '--output_dfname', required=True)

    parser.add_argument('-t', '--type', required=True, help="Args.type must be m (matrix) or d (pandas dataframe)")
    # if "m", then pickle concatenate
    # if "d", then pandas
    # if "a", then pickle addition

    parser.add_argument('-if', '--iter_free', type=int, default=1, help="Iteratively free the arguments while integrating." + \
                                                                        "Only for type free.")

    parser.add_argument('-a', '--axis', type=int, default = -1)

    parser.add_argument('-o', '--int_name_dfname', required=True)
    return parser




def run(args):

    # load the data

    output_df = pd.read_csv(args.output_dfname, sep="\t")
    int_name_df = pd.read_csv(args.int_name_dfname, sep="\t")
    print(int_name_df.head())

    assert set(output_df.columns.values) == set(int_name_df.columns.values)

    for x in output_df:
        filenames = output_df[x].values

        print()
        print(x)
        print("Files to integrate are ", filenames)
        integrated_filename = int_name_df[x].values[0]
        print("Integrated_file is ", integrated_filename)


        if args.type == "m":
            pj.concatenate_matrices(filenames, concatenated_filename=integrated_filename, axis=args.axis)
        elif args.type == "d":
            pj.integrate_dfs(filenames, integrated_filename=integrated_filename)
        elif args.type == "a":
            pj.add_matrices(filenames, added_filename=integrated_filename, iter_free=args.iter_free)

        else:
            raise ValueError("Args.type must be m (matrix concatenation), d (pandas dataframe), or a (matrix addition")



def main():
    run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()
