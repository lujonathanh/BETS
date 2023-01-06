__author__ = 'jlu96'


import os
import sys

def get_parser():
    # Parse arguments
    import argparse

    description = 'Get the new filename of network for della'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-f', '--filename', required=True)

    return parser

def run(args):

    f = args.filename

    names = f.split("-")

    genes = names[0]
    norm = names[1]
    reps = names[2]
    deg = names[3]
    test = names[5]
    lag = names[6]
    null = names[7]
    coef = names[8]
    fdr_thresh = names[10]
    fdr_type = names[11]
    ending = "-".join(names[12:])

    if fdr_type == "none":
        fdr_type = "g"
    elif fdr_type == "effect":
        fdr_type = "l"
    else:
        raise ValueError("FDR type must be none or effect. Given is "  + fdr_type)

    if norm == "0mean1var":
        norm = "0mean-1var"
    elif norm == "0mean":
        norm = "0mean-unnormalized"
    else:
        raise ValueError("Norm must be 0mean1var or 0mean")

    prefix = "-".join([genes, deg, reps])
    folder = "_".join([norm, null + "-null", fdr_type + "-fdr-" + fdr_thresh])
    suffix = "-".join([test, lag, coef, ending])

    newname = "_".join([prefix, folder, suffix])

    print(newname)




def main():
    run(get_parser().parse_args(sys.argv[1:]))




if __name__ == '__main__':
    main()