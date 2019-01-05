__author__ = 'jlu96'



import stochpy
import numpy as np
import sys
import csv
import os
import collections




def get_sim(filename, filedir=None, outdir=None, argDict=dict(), paramDict=dict(), interval=0.5,
            suffix='_TS_data.txt'):
    """


    :param filename: name of .psc file (specifying diff eq's, values)
    :param filedir: directory of .psc file
    :param outdir: directory of output file

    :param argDict: dict containing keyword args for the smod.Model() call
    :param paramDict: dict containing parameter values. This can be empty;
        otherwise the default values are used
    :param interval: interval to keep in our array. For example, if the simulation
    produces a value for each second (which we can't control), setting interval to 5
    will result in a vector spaced at 5 seconds.
    :param suffix: any suffix to append
    :return: smod, newdictToVector, times

    smod: the original stochpy simulation model
    newdictToVector: a dictionary of form:
        "A" -> vector of "A"'s values
        "time" -> vector of time
    times: the timepoints of our simulation
    """
    smod = stochpy.SSA()
    smod.Model(filename, filedir)
    if paramDict:

        # suffix will note the parameters that were set
        for param in paramDict:
            smod.ChangeParameter(param, paramDict[param])
            val = paramDict[param]
            if val != 0:
                exp = int(np.log10(val))
                new_val = int(val * 1.0 / 10**exp)
                suffix += '_' + param + '_' + str(new_val) + 'e' + str(exp)
            else:
                suffix += '_' + param + '_' + str(val) + 'e0'
        suffix += '.txt'

    smod.DoStochSim(**argDict)

    # species was simulated
    smod.Export2File(analysis='timeseries', datatype='species', directory=outdir)

    # now load it into dictionary vector
    import csv
    with open(os.path.join(outdir, filename + '_species_timeseries1.txt'), 'rU') as csvfile:
        csvfile.readline()
        reader = csv.DictReader(csvfile, delimiter='\t')
        dictToVector = collections.OrderedDict()
        for field in reader.fieldnames:
            dictToVector[field] = []
        for row in reader:
            for field in dictToVector:
                dictToVector[field].append(eval(row[field]))

    for field in dictToVector:
        dictToVector[field] = np.array(dictToVector[field])


    max_time = max(dictToVector['Time'])

    # get the spaced intervals we want
    times = np.arange(0, max_time, interval)
    indices = [np.nonzero(dictToVector['Time'] > t)[0][0] for t in times]

    newdictToVector = collections.OrderedDict()
    for key in dictToVector:
        newdictToVector[key] = dictToVector[key][np.ix_(indices)]


    # Save
    out_file = os.path.join(outdir, filename + suffix)
    with open(out_file, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        writer.writerow([key for key in newdictToVector])
        for i in range(len(indices)):
            writer.writerow([newdictToVector[key][i] for key in newdictToVector])
    print "Data written to ", out_file

    return smod, newdictToVector, times


def write_simulation(dictToVector, times, filename):
    fieldnames = ['Time'] + dictToVector.keys()


    with open(filename, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for t in range(len(times)):
            outrow = {}
            for field in fieldnames:
                if field == 'Time':
                    outrow[field] = times[t]
                else:
                    outrow[field] = dictToVector[field][t]
            writer.writerow(outrow)



def run(args):
    if not args.simulationfile:
        args.simulationfile = args.pscfile + "-simulations.txt"


    with open(args.argfile, 'rU') as csvfile:

        reader = csv.DictReader(csvfile, delimiter='\t')
        argDict = reader.next()
        for field in argDict:
            argDict[field] = eval(argDict[field])


    with open(args.paramfile, 'rU') as csvfile:

        reader = csv.DictReader(csvfile, delimiter='\t')
        paramDict = reader.next()
        for field in paramDict:
            paramDict[field] = eval(paramDict[field])

    smod, newdictToVector, times = get_sim(args.pscfile, filedir=args.filedir, outdir=args.outfiledir,
                                           argDict=argDict, paramDict=paramDict, interval=args.interval)


    write_simulation(newdictToVector, times, args.simulationfile)



def get_parser():
    # Parse arguments
    import argparse
    description = 'Get the gene dynamical simulation of the given file with given parameters.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-psc', '--pscfile', default="mugler.psc")
    parser.add_argument('-fdr', '--filedir', default='/Users/jlu96/v-causal-snps/code/models/')
    parser.add_argument('-ofdr', '--outfiledir', default='/Users/jlu96/v-causal-snps/data/')
    parser.add_argument('-i', '--interval', type=float, default=100.0)
    parser.add_argument('-a', '--argfile', default='arg.txt')
    parser.add_argument('-p', '--paramfile', default='param.txt')
    parser.add_argument('-sf', '--simulationfile', default=None)

    return parser


def main():
    run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()