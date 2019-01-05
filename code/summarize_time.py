__author__ = 'jlu96'

import pandas as pd
import sys
from dateutil.parser import parse



def get_parser():
    # Parse arguments
    import argparse

    description = 'Aggregate timing files'
    parser = argparse.ArgumentParser(description=description)


    parser.add_argument('-i', '--timefile',  required=True)

    parser.add_argument('-o', '--outfile', required=True)

    parser.add_argument('-oo', '--overallfile', required=True)

    return parser



def run(args):
    """


    :param filenames: Names of csv files that have "Name", "Start", "End", and "Elapsed" fields
    :return:
    """

    print "Assumes that the first row is the start, last row is the end."


    filenames = []
    starts = []
    ends = []
    elapseds = []
    nscripts = []

    with open(args.timefile, 'rU') as f:
        for l in f.readlines():
            filename = l.split("\n")[0]
            print filename
            time_df = pd.read_csv(filename)

            filenames.append(filename)
            starts.append(time_df["Start"].values[0])
            ends.append(time_df["End"].values[-1])
            elapseds.append(time_df["Elapsed"].sum())

            if "Nscripts" not in time_df:
                nscripts.append(time_df.shape[0])
            else:
                nscripts.append(time_df["Nscripts"].sum())


    out_df = pd.DataFrame()
    out_df["Name"] = filenames
    out_df["Start"] = starts
    out_df["End"] = ends
    out_df["Real_Seconds"] = [(parse(row["End"]) - parse(row["Start"])).total_seconds() for i, row in out_df.iterrows()]
    out_df["Real_Days"] = [round( (x / (3600 * 24.)), 1) for x in out_df["Real_Seconds"]]
    out_df["Elapsed"] = elapseds
    out_df["Nscripts"] = nscripts

    print "Summary: ", out_df
    print "Writing to ", args.outfile
    out_df.to_csv(args.outfile, index=0)

    overall_df = pd.DataFrame()
    overall_df["Start"] = [starts[0]]
    overall_df["End"] = [ends[-1]]
    overall_df["Real_Seconds"] = [(parse(row["End"]) - parse(row["Start"])).total_seconds() for i, row in overall_df.iterrows()]
    overall_df["Real_Days"] = [round((x / (3600 * 24.)), 1) for x in overall_df["Real_Seconds"]]
    overall_df["Elapsed"] = [sum(elapseds)]
    overall_df["Nscripts"] = [sum(nscripts)]

    print "Overall: ", overall_df
    print "Writing to ", args.overallfile
    overall_df.to_csv(args.overallfile, index=0)




def main():
    run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()
