__author__ = 'jlu96'

import prep_jobs as pj
import sys
import os
import pickle
import csv
import pandas as pd
import math
import collections
import itertools
import geneTSmunging as gtm

def get_parser():
    # Parse arguments
    import argparse

    description = 'Prepare cluster jobs by partitioning tests by rows and hyper-parameters.'
    parser = argparse.ArgumentParser(description=description)


    parser.add_argument('-d', '--data_file', required=True)

    parser.add_argument('-d2', '--rand_data_file', required=True, help="The effect genes")

    parser.add_argument('-lr', '--load_reps', required=True, type=int)

    parser.add_argument('-o', '--output_name', required=True)

    parser.add_argument('-hlf', '--hyper_list_file', required=True)

    parser.add_argument('-t', '--test', required=True)

    parser.add_argument('-tn', '--test_name', required=True)

    parser.add_argument('-sn', '--script_num', type=int, default=3)

    parser.add_argument('-p', '--parallel_num', type=int, default=0)

    parser.add_argument('-l', '--lag', type=int, required=True)

    parser.add_argument('-n', '--null', type=str, required=True)

    parser.add_argument('-cv', '--cv', type=int, default=1, help="Do prep with CV or not. If 0, then skip the CV-making step.")

    parser.add_argument('-oa', '--only_array', type=int, default=0, help="Whehter to only save output coefs as arrays (1) or as whole matrices that are integrated by adding (0)")


    return parser

def run(args):
    if args.test not in {"r", "l", "e"}:
        raise ValueError("args.test must be r (ridge), l (lasso) or e (elastic net)")

    if args.null not in {"l", "g"}:
        raise ValueError("args.null must be l (local) or g (global)")

    # Load files
    data_file = args.data_file
    rand_data_file = args.rand_data_file

    if args.load_reps:
        genes, geneTS = gtm.load_basic_rep_file_list(data_file)
        #dfs, genes, geneTS, df, __, __ = gtm.load_rep_file_list(data_file)
    else:
        df = pd.read_csv(data_file, sep="\t")
        genes, geneTS = gtm.get_gene_TS(df)
    n = len(genes)

    hyperlist = pickle.load(open(args.hyper_list_file, 'rb'))
    # hyper_names = cp.hyperlist_to_namelist(hyperlist)


    # Make hyper files for cross_validate loading.

    hyper_filenames = []

    print("*************")
    print("HYPERS")
    print("*************")

    if not os.path.exists("hyper"):
        os.makedirs("hyper")


    # for hyper, hyper_name in zip(hyperlist, hyper_names):
    for hyper, h in zip(hyperlist, list(range(len(hyperlist)))):
        hyper_filename = "hyper" + os.sep + args.output_name + "-hyper-" + str(h) + ".p"

        hyper_filenames.append(hyper_filename)

        pickle.dump([hyper], open(hyper_filename, 'wb'))

    print("Hypers written in format: ", hyper_filename)


    # Make row files
    # Split up the rows according to number of input scripts
    partition_rows = pj.partition_inputs(list(range(n)), args.script_num)

    row_filenames = []


    print("*************")
    print("ROWS")
    print("*************")


    if not os.path.exists("rows"):
        os.makedirs("rows")

    for partition_row, i in zip(partition_rows, list(range(len(partition_rows)))):

        row_filename = os.path.join("rows", args.output_name + "-row-" + str(i) + ".p")
        row_filenames.append(row_filename)

        pickle.dump(partition_row, open(row_filename, 'wb'))

    print("Row written in format: ", row_filename)


    if not os.path.exists("timing"):
        os.makedirs("timing")
        print("Folder timing created")
    resulttimefile = os.path.join("timing", "result_time.csv")
    if not os.path.exists(resulttimefile):
        with open(resulttimefile, 'w') as csvfile:
            f = csv.writer(csvfile)
            f.writerow(["Name", "Start", "End", "Elapsed"])


    if args.cv != 0:
        print("*************")
        print("CV")
        print("*************")

        # Make CV scripts

        cv_scripts = []

        hyper_output_dict = collections.OrderedDict()
        hyper_int_dict = collections.OrderedDict()

        if not os.path.exists("cv-scripts"):
            os.makedirs("cv-scripts")


        cvtimefile = os.path.join("timing", "hyper_time.csv")
        if not os.path.exists(cvtimefile):
            with open(cvtimefile, 'w') as csvfile:
                f = csv.writer(csvfile)
                f.writerow(["Name", "Start", "End", "Elapsed"])



        for hyper, h, hyper_filename in zip(hyperlist, list(range(len(hyperlist))), hyper_filenames):

            hyper_output_group = []

            for partition_row, i, row_filename in zip(partition_rows, list(range(len(partition_rows))), row_filenames):

                cv_prefix = args.output_name + "-cv-" + str(h) + "-row-" + str(i)

                cv_script = os.path.join("cv-scripts", cv_prefix + ".sh")
                cv_scripts.append(cv_script)

                cv_output  = "hyper" + os.sep + cv_prefix + "-result.txt"
                hyper_output_group.append(cv_output)

                command_string = "time python cross_validate.py -d " + data_file + " -lr " + str(args.load_reps) +  " -o " + cv_output + " -hl " + str(hyper_filename) \
                                 + " -t " + args.test + " -l " + str(args.lag) + " -rl " + str(row_filename)


                with open(cv_script, 'w') as outputfile:
                    outputfile.write("#!/bin/bash\n")
                    outputfile.write("START=$(date)\n")
                    #outputfile.write("module load python/2.7\n")
                    # outputfile.write("module load python/2.7/scipy-mkl\n")
                    # outputfile.write("module load python/2.7/numpy-mkl\n")
                    #outputfile.write("module load anaconda\n")
                    outputfile.write("module load anaconda3\n")
                    outputfile.write(command_string)
                    outputfile.write("\n")
                    outputfile.write("END=$(date)\n")
                    outputfile.write("echo " + cv_script + ",$START,$END,$SECONDS >> " + cvtimefile + "\n")
                os.chmod(cv_script, 0o777)


            # Set the output names, prepare for integration of all the hyper parameter fit results
            hyper_output_dict[str(hyper)] = hyper_output_group
            hyper_int_dict[str(hyper)] = "hyper" + os.sep + args.output_name + "-cv-" + str(h) + "-result.txt"



        hyper_output_df = pd.DataFrame(hyper_output_dict)
        hyper_int_df = pd.DataFrame(hyper_int_dict, index=[0])

        print("Hyper output df is in form", hyper_output_df.head(n=5))

        hyper_output_df.to_csv("cv_outputs.txt", sep="\t", index=0)
        hyper_int_df.to_csv("cv_integrated.txt", sep="\t", index=0)

        print("Partitioned CV fit_result_dfs in cv_outputs.txt", "Integrated CV fit_result_dfs in cv_integrated.txt")

        with open("cv_script_list.txt", 'w') as outfile:
            for cv_script in cv_scripts:
                outfile.write(cv_script + "\n")
            print("CV scripts written to cv_script_list.txt")

        if args.parallel_num > 0:
            print("Parallel Number (# processes per job): " + str(args.parallel_num))

            script_groups = pj.partition_inputs(cv_scripts, number=int(math.ceil(len(cv_scripts) * 1.0/args.parallel_num)))

            print("Number of script groups ", len(script_groups))

            parallel_scripts = []
            for i, script_group in zip(list(range(len(script_groups))), script_groups):
                appended_script_filenames = ["./" + script_filename for script_filename in script_group]
                parallel_script = " & ".join(appended_script_filenames)
                parallel_scripts.append(parallel_script)

            with open("cv_parallel_script_list.txt", 'w') as scriptfile:
                for parallel_script in parallel_scripts:
                    scriptfile.write(parallel_script + "\n")
                print("Parallel script list written to cv_parallel_script_list.txt")






        # Integrate hyperparameters
        # Begin whole normal fit

        hyper_script = "set_hyper.sh"

        with open(hyper_script, 'w') as outputfile:
            outputfile.write("#!/bin/bash\n")
            outputfile.write("START=$(date)\n")
            outputfile.write("set -e\n")
            outputfile.write("time python integrate_hyper.py -hfd cv_outputs.txt -ind cv_integrated.txt -hl " + args.hyper_list_file + "\n")
            outputfile.write("time python set_hyper.py -ind cv_integrated.txt -r " + "hyper" + os.sep + "hyper_df.txt -o " + "hyper" +
                             os.sep + "best_hyper.p -hl " + args.hyper_list_file + " -tn " + args.test_name + " \n")
            outputfile.write("END=$(date)\n")
            outputfile.write("echo " + hyper_script + ",$START,$END,$SECONDS >> " + resulttimefile + "\n")
        os.chmod(hyper_script, 0o777)

        print("set_hyper.sh written")




    print("*************")
    print("FITTING")
    print("*************")


    # Run the actual fit
    if not os.path.exists("fit"):
        os.makedirs("fit")

    if not os.path.exists("fit-scripts"):
        os.makedirs("fit-scripts")

    fittimefile = os.path.join("timing", "fit_time.csv")
    if not os.path.exists(fittimefile):
        with open(fittimefile, 'w') as csvfile:
            f = csv.writer(csvfile)
            f.writerow(["Name", "Start", "End", "Elapsed"])


    fit_scripts = []
    fit_output_prefixes = []
    for partition_row, i, row_filename in zip(partition_rows, list(range(len(partition_rows))), row_filenames):

        fit_prefix = args.output_name + "-fit-row-" + str(i)

        fit_script = os.path.join("fit-scripts", fit_prefix + ".sh")
        fit_scripts.append(fit_script)

        fit_output_prefix  = "fit" + os.sep + fit_prefix
        fit_output_prefixes.append(fit_output_prefix)


        command_string = "time python fit_all.py -d " + data_file + " -rd " + rand_data_file + " -lr " + str(args.load_reps) + \
                         " -o " + fit_output_prefix + " -bh " + \
                        "hyper" + os.sep + "best_hyper.p" + " -t " + args.test + " -l " + str(args.lag) + " -rl " + \
                         str(row_filename) + " -n " + args.null + " -oa " + str(args.only_array)


        with open(fit_script, 'w') as outputfile:
            outputfile.write("#!/bin/bash\n")
            outputfile.write("START=$(date)\n")
            #outputfile.write("module load python/2.7\n")
            # outputfile.write("module load python/2.7/scipy-mkl\n")
            # outputfile.write("module load python/2.7/numpy-mkl\n")
            outputfile.write("module load anaconda3\n")
            outputfile.write(command_string)
            outputfile.write("\n")
            outputfile.write("END=$(date)\n")
            outputfile.write("echo " + fit_script + ",$START,$END,$SECONDS >> " + fittimefile + "\n")
        os.chmod(fit_script, 0o777)


    with open("fit_script_list.txt", 'w') as outfile:
        for fit_script in fit_scripts:
            outfile.write("./" + fit_script + "\n")
        print("Fit scripts written to fit_script_list.txt")

    if args.parallel_num > 0:
        print("Parallel Number (# processes per job): " + str(args.parallel_num))

        script_groups = pj.partition_inputs(fit_scripts, number=int(math.ceil(len(fit_scripts) * 1.0/args.parallel_num)))

        print("Number of script groups ", len(script_groups))

        parallel_scripts = []
        for i, script_group in zip(list(range(len(script_groups))), script_groups):
            appended_script_filenames = ["./" + script_filename for script_filename in script_group]
            parallel_script = " & ".join(appended_script_filenames)
            parallel_scripts.append(parallel_script)

        with open("fit_parallel_script_list.txt", 'w') as scriptfile:
            for parallel_script in parallel_scripts:
                scriptfile.write(parallel_script + "\n")
            print("Parallel script list written to fit_parallel_script_list.txt")





    # Note the output files

    fit_coefs = [fit_output_prefix + "_coefs.p" for fit_output_prefix in fit_output_prefixes]
    fit_intercepts = [fit_output_prefix + "_intercepts.p" for fit_output_prefix in fit_output_prefixes]
    fit_results = [fit_output_prefix + "_fit_result_df.txt" for fit_output_prefix in fit_output_prefixes]
    fit_coefsr = [fit_output_prefix + "_coefsr.p" for fit_output_prefix in fit_output_prefixes]
    # fit_interceptsr = [fit_output_prefix + "_interceptsr.p" for fit_output_prefix in fit_output_prefixes]
    fit_resultsr = [fit_output_prefix + "_fit_result_dfr.txt" for fit_output_prefix in fit_output_prefixes]

    fit_output_dict = collections.OrderedDict()
    fit_output_dict["coef"] = fit_coefs
    fit_output_dict["coefr"] = fit_coefsr
    fit_output_dict["intercept"] = fit_intercepts
    # fit_output_dict["interceptr"] = fit_interceptsr

    output_matr_df = pd.DataFrame(fit_output_dict)
    output_matr_df.to_csv("output_matr_list.txt", sep="\t", index=False)
    print("Output matrices written to output_matr_list.txt")

    int_matr_dict = collections.OrderedDict()
    int_matr_dict["coef"] = "fit" + os.sep + args.output_name + "_coefs.p"
    int_matr_dict["coefr"] = "fit" + os.sep + args.output_name + "_coefsr.p"
    int_matr_dict["intercept"] = "fit" + os.sep + args.output_name + "_intercepts.p"
    # int_matr_dict["interceptr"] = "fit" + os.sep + args.output_name + "_interceptsr.p"

    int_matr_df = pd.DataFrame(int_matr_dict, index=[0])
    int_matr_df.to_csv("int_matr_list.txt", sep="\t", index=False)
    print("integrated matrices written to int_matr_list.txt")


    fit_result_dict = collections.OrderedDict()
    fit_result_dict["fit_result"] = fit_results
    fit_result_dict["fit_resultr"] = fit_resultsr

    output_df_df = pd.DataFrame(fit_result_dict)
    output_df_df.to_csv("output_df_list.txt", sep="\t", index=False)
    print("output dfs written to output_df_list.txt")


    int_df_dict = collections.OrderedDict()
    int_df_dict["fit_result"] = "fit" + os.sep + args.output_name + "_fit_result_df.txt"
    int_df_dict["fit_resultr"] = "fit" + os.sep + args.output_name + "_fit_result_dfr.txt"

    int_df_df = pd.DataFrame(int_df_dict, index=[0])
    int_df_df.to_csv("int_df_list.txt", sep="\t", index=False)
    print("Integrated dfs written to int_df_list.txt")


    with open("finish-none.sh", 'w') as ifile:
        ifile.write("#!/bin/bash\n")
        ifile.write("START=$(date)\n")
        ifile.write("set -e\n")
        ifile.write("time python integrate_outputs_rand_row.py -i output_matr_list.txt -o int_matr_list.txt " + (" -t m -a 1 " if args.only_array else " -t a "))
        ifile.write(" && " + \
                    "time python integrate_outputs_rand_row.py -i output_df_list.txt -o int_df_list.txt -t d " + "\n")
        ifile.write("time python get_result_coef.py -df " + data_file + " -rdf " + rand_data_file +\
                    " -lr " + str(args.load_reps) + \
                    " -bh " + "hyper" + os.sep + "best_hyper.p" + \
                    " -o " + \
                    args.output_name + " -cf " +  int_matr_dict["coef"] + " -if " + int_matr_dict["intercept"] + \
                    " -cfr " + int_matr_dict["coefr"] + " -fr " + \
                    int_df_dict["fit_result"] + " -frr " + int_df_dict["fit_resultr"] + " -l " + str(args.lag) + \
                    " -sb " + "n" + " -tn " + args.test_name + "\n")
        ifile.write("END=$(date)\n")
        ifile.write("echo " + "finish-none.sh" + ",$START,$END,$SECONDS >> " + resulttimefile + "\n")
        print("Finish script, stratby None, written to finish-none.sh")
        os.chmod("finish-none.sh", 0o777)

    with open("finish-effect.sh", 'w') as ifile:
        ifile.write("#!/bin/bash\n")
        ifile.write("START=$(date)\n")
        ifile.write("set -e\n")
        ifile.write("time python integrate_outputs_rand_row.py -i output_matr_list.txt -o int_matr_list.txt " + (" -t m -a 1 " if args.only_array else " -t a "))
        ifile.write(" && " + \
                    "time python integrate_outputs_rand_row.py -i output_df_list.txt -o int_df_list.txt -t d " + "\n")
        ifile.write("time python get_result_coef.py -df " + data_file + " -rdf " + rand_data_file +\
                    " -lr " + str(args.load_reps) + \
                    " -bh " + "hyper" + os.sep + "best_hyper.p" + \
                    " -o " + \
                    args.output_name + " -cf " +  int_matr_dict["coef"] + " -if " + int_matr_dict["intercept"] + \
                    " -cfr " + int_matr_dict["coefr"] + " -fr " + \
                    int_df_dict["fit_result"] + " -frr " + int_df_dict["fit_resultr"] + " -l " + str(args.lag) + \
                    " -sb " + "e" + " -tn " + args.test_name + "\n")
        ifile.write("END=$(date)\n")
        ifile.write("echo " + "finish-effect.sh" + ",$START,$END,$SECONDS >> " + resulttimefile + "\n")

        print("Finish script, stratby effect, written to finish-effect.sh")
        os.chmod("finish-effect.sh", 0o777)






    with open("plot_coef.sh", 'w') as ifile:
        ifile.write("#!/bin/bash\n")
        ifile.write("START=$(date)\n")
        ifile.write("time python get_result_coef.py -df " + data_file + " -rdf " + rand_data_file +\
                    " -lr " + str(args.load_reps) + \
                    " -bh " + "hyper" + os.sep + "best_hyper.p" + \
                    " -o " + \
                    args.output_name + " -cf " +  int_matr_dict["coef"] + " -if " + int_matr_dict["intercept"] + \
                    " -cfr " + int_matr_dict["coefr"]  + " -fr " + \
                    int_df_dict["fit_result"] + " -frr " + int_df_dict["fit_resultr"] + " -l " + str(args.lag) + \
                    " -sb " + "n" + " -tn " + args.test_name +  " -pcf 1 " + "\n")
        ifile.write("END=$(date)\n")
        ifile.write("echo " + "plot_coef.sh" + ",$START,$END,$SECONDS >> " + resulttimefile + "\n")

        print("Plot coef script written to plot_coef.sh")
        os.chmod("plot_coef.sh", 0o777)



    with open("cleanup_list.txt", 'w') as outfile:
        cleanup_list = row_filenames
        if args.cv:
            cleanup_list += cv_scripts + list(itertools.chain.from_iterable(list(hyper_output_dict.values())))

        cleanup_list += fit_scripts + fit_coefs + fit_intercepts + fit_results + fit_coefsr  + fit_resultsr
        for script in cleanup_list:
            outfile.write(script + "\n")
        print("Cleanup scripts written to cleanup_list.txt")


    with open("timing/timing_list.txt", 'w') as outfile:
        outfile.write(cvtimefile + "\n")
        outfile.write(fittimefile + "\n")
        outfile.write(resulttimefile + "\n")
    print("Timing files written to timing_list.txt")

    with open("summarize_time.sh", 'w') as outfile:
        outfile.write("python summarize_time.py -i timing/timing_list.txt -o timing/summary_time.csv -oo timing/overall_time.csv\n")
    os.chmod("summarize_time.sh", 0o777)
    print("Summarize timing script written to summarize_time.sh")






def main():
    run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()
