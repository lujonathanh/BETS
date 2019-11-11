__author__ = 'jlu96'

import sys

def get_parser():
    # Parse arguments
    import argparse

    description = 'Get all variables declared, assuming they are exported. Designed for package_params_cpipeline.sh'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-e', '--envfile', required=True, help="List the environment variables resulting from printenv")

    parser.add_argument('-p', '--paramfile', required=True)

    parser.add_argument('-o', '--outputfile', required=True)

    return parser

def run(args):

    keys = []
    print("Loading params: ", args.paramfile)
    with open(args.paramfile, 'r') as f:
        for line in f.readlines():
            keys.append(line.split("\n")[0])

    env_dict = {}
    print("Getting values: ", args.envfile)
    with open(args.envfile, 'r') as e:
        for line in e.readlines():
            key = line.split("=")[0]
            value = line.split("\n")[0].split("=")[1]

            env_dict[key] = value

    print("Writing to: ", args.outputfile)
    with open(args.outputfile, 'w') as o:
        for key in keys:
            o.write(key + "," + env_dict[key] + "\n")








def main():
    run(get_parser().parse_args(sys.argv[1:]))




if __name__ == '__main__':
    main()