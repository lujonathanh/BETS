__author__ = 'jlu96'

import sys

def get_parser():
    # Parse arguments
    import argparse

    description = 'Get all variables declared, assuming they are exported. Designed for package_params_cpipeline.sh'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-f', '--filename', required=True)

    parser.add_argument('-o', '--outputname', required=True)

    return parser

def run(args):

    keys = []
    print "Loading script: ", args.filename
    with open(args.filename, 'rU') as f:
        for line in f.readlines():
            if line.startswith("export"):
                newline = line[len("export"):].strip().split("#")[0].rstrip()
                key_value = newline.split("=")
                keys.append(key_value[0])

    print "Writing declared vars: ", args.outputname
    with open(args.outputname, 'w') as o:
        for key in keys:
            o.write(key + "\n")





def main():
    run(get_parser().parse_args(sys.argv[1:]))




if __name__ == '__main__':
    main()