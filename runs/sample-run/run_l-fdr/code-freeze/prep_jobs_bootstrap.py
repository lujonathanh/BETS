__author__ = 'jlu96'

import prep_jobs as pj
import sys
import os
import pandas as pd
import math
import collections
import geneTSmunging as gtm
import csv

def get_parser():
    # Parse arguments
    import argparse

    description = 'Prepare cluster jobs by partitioning tests by rows and hyper-parameters.'
    parser = argparse.ArgumentParser(description=description)


    parser.add_argument('-d', '--data_file', required=True)

    parser.add_argument('-d2', '--rand_data_file', required=True, help="The effect genes")

    parser.add_argument('-lr', '--load_reps', required=True, type=int)

    parser.add_argument('-o', '--output_name', required=True)

    parser.add_argument('-t', '--test', required=True)

    parser.add_argument('-tn', '--test_name', required=True)

    parser.add_argument('-sn', '--script_num', type=int, default=3,
                        help="Number of scripts, total, that need to be run")

    parser.add_argument('-p', '--parallel_num', type=int, default=0)

    parser.add_argument('-l', '--lag', type=int, required=True)

    parser.add_argument('-n', '--null', type=str, required=True)

    parser.add_argument('-bn', '--bootstrap_num', type=int, required=True)

    parser.add_argument('-oa', '--only_array', type=int, default=1, help="Whehter to only save output coefs as arrays (1) or as whole matrices that are integrated by adding (0)")

    parser.add_argument('-waf', '--write_all_bootstrap_scripts_first',
                        type=int, default=1)

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




    # Make row files
    # Split up the rows according to number of input scripts
    partition_rows = pj.partition_inputs(list(range(n)), args.script_num)

    row_filenames = []


    print("*************")
    print("ROWS")
    print("*************")

    for partition_row, i in zip(partition_rows, list(range(len(partition_rows)))):

        row_filename = os.path.join("rows", args.output_name + "-row-" + str(i) + ".p")
        row_filenames.append(row_filename)

    print("Reading rows from format: ", row_filename)

    print("*************")
    print("BOOTSTRAP")
    print("*************")


    # Run the actual fit
    # Need an integration
    if not os.path.exists("bootstrap"):
        os.makedirs("bootstrap")

    # For the bootstrap individual fit scripts
    if not os.path.exists("bootstrap-fit-scripts"):
        os.makedirs("bootstrap-fit-scripts")


    # For the bootstrap finish scripts
    if not os.path.exists("bootstrap-finish-scripts"):
        os.makedirs("bootstrap-finish-scripts")

    # Finish, aggregating all the coefficients (stratification = none)
    if not os.path.exists(os.path.join("bootstrap-finish-scripts", "none")):
        os.makedirs(os.path.join("bootstrap-finish-scripts", "none"))

    # Finish, stratifying each coefficient by the effect gene (stratification = effect)
    if not os.path.exists(os.path.join("bootstrap-finish-scripts", "effect")):
        os.makedirs(os.path.join("bootstrap-finish-scripts", "effect"))








    # if args.write_all_bootstrap_scripts_first:

    print("WRITING ALL THE SCRIPTS INITIALLY!!!!!! NOTE the list will be written before all the files are written!!!")

    for b in range(args.bootstrap_num):
        if not os.path.exists(os.path.join("bootstrap-fit-scripts", str(b))):
            os.makedirs(os.path.join("bootstrap-fit-scripts", str(b)))

    all_bootstrap_scripts = [os.path.join("bootstrap-fit-scripts", str(b), args.output_name + "-bootstrap-" + str(b) + "-row-" + str(i) + ".sh")
                             for b in range(args.bootstrap_num) for i in range(len(row_filenames))]


    print("SCRIPTS")

    with open("bootstrap_script_list.txt", 'w') as outfile:
        for bootstrap_script in all_bootstrap_scripts:
            outfile.write("./" + bootstrap_script + "\n")
        print("bootstrap scripts written to bootstrap_script_list.txt")

        if args.parallel_num > 0:
            print("Parallel Number (# processes per job): " + str(args.parallel_num))

            script_groups = pj.partition_inputs(all_bootstrap_scripts, number=int(math.ceil(len(all_bootstrap_scripts) * 1.0/args.parallel_num)))

            print("Number of script groups ", len(script_groups))

            parallel_scripts = []
            for i, script_group in zip(list(range(len(script_groups))), script_groups):
                appended_script_filenames = ["./" + script_filename for script_filename in script_group]
                parallel_script = " & ".join(appended_script_filenames)
                parallel_scripts.append(parallel_script)

            with open("bootstrap_parallel_script_list.txt", 'w') as scriptfile:
                for parallel_script in parallel_scripts:
                    scriptfile.write(parallel_script + "\n")
                print("Parallel script list written to bootstrap_parallel_script_list.txt")









    # make one script for each...

    # all_bootstrap_scripts = set([])

    all_int_coefs = []
    all_int_intercepts = []

    finish_none_scripts = []
    finish_effect_scripts = []

    # record where the thresholded coefficients are written
    # For integrating these, later.
    fdrs = [0.01, 0.05, 0.1, 0.2]
    all_fdr_none_coefs_dict = dict([(x, []) for x in fdrs])
    all_fdr_effect_coefs_dict = dict([(x, []) for x in fdrs])

    all_fdr_none_intercepts_dict = dict([(x, []) for x in fdrs])
    all_fdr_effect_intercepts_dict = dict([(x, []) for x in fdrs])



    try:
        fittimefile = os.path.join("timing", "bootstrap_fit_time.csv")
        if not os.path.exists(fittimefile):
            with open(fittimefile, 'w') as csvfile:
                f = csv.writer(csvfile)
                f.writerow(["Name", "Start", "End", "Elapsed"])


        finishtimefile = os.path.join("timing", "bootstrap_finish_time.csv")
        if not os.path.exists(finishtimefile):
            with open(finishtimefile, 'w') as csvfile:
                f = csv.writer(csvfile)
                f.writerow(["Name", "Start", "End", "Elapsed"])

        resulttimefile = os.path.join("timing", "bootstrap_result_time.csv")
        if not os.path.exists(resulttimefile):
            with open(resulttimefile, 'w') as csvfile:
                f = csv.writer(csvfile)
                f.writerow(["Name", "Start", "End", "Elapsed"])

        with open(os.path.join("timing/timing_list.txt"), 'a') as f:
            f.write(fittimefile + "\n")
            f.write(finishtimefile + "\n")
            f.write(resulttimefile + "\n")


    except IOError:
        raise IOError("the timing folder does not exist. Please run ./prep_jobs_rand_cv.sh first.")


    for b in range(args.bootstrap_num):
        if b % 50 == 0:
            print("SEED/BOOTSTRAP NUM: ", b)

        bootstrap_outmost_name = args.output_name + "-bootstrap-" + str(b)

        bootstrap_folder = os.path.join("bootstrap", str(b))
        if not os.path.exists(bootstrap_folder):
            os.makedirs(bootstrap_folder)
        # print "Created folder: ", bootstrap_folder

        bootstrap_outmost_prefix = os.path.join(bootstrap_folder, bootstrap_outmost_name)



        if not os.path.exists(os.path.join("bootstrap-fit-scripts", str(b))):
            os.makedirs(os.path.join("bootstrap-fit-scripts", str(b)))


        # create scripts for bootstrap
        bootstrap_scripts = [os.path.join("bootstrap-fit-scripts", str(b), bootstrap_outmost_name + "-row-" + str(i) + ".sh")
                             for i in range(len(partition_rows))]
        bootstrap_row_prefixes = [bootstrap_outmost_prefix + "-row-" + str(i) for i in range(len(partition_rows))]

        command_template = "time python3 fit_bootstrap.py -d " + data_file + " -rd " + rand_data_file + " -lr " + str(args.load_reps) + \
                             " -o " + "bootstrap_row_prefixes[i]" + " -bh " + \
                            "hyper" + os.sep + "best_hyper.p" + " -t " + args.test + " -l " + str(args.lag) + " -rl " + \
                             "row_filename" + " -n " + args.null + " -s " + str(b) + " -oa " + str(args.only_array)

        for i, row_filename in zip(list(range(len(partition_rows))), row_filenames):

            # writing results to the bootstrap prefix

            command_string = command_template.replace("bootstrap_row_prefixes[i]", bootstrap_row_prefixes[i]).replace("row_filename", row_filename)

            with open(bootstrap_scripts[i], 'w') as outputfile:
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
                    outputfile.write("echo " + bootstrap_scripts[i] + ",$START,$END,$SECONDS >> " + fittimefile + "\n")
            os.chmod(bootstrap_scripts[i], 0o777)


        # print "Scripts made"

        # all_bootstrap_scripts = all_bootstrap_scripts.union(set(bootstrap_scripts))

        # Note the output files

        bootstrap_coefs = [bootstrap_row_prefix + "_coefs.p" for bootstrap_row_prefix in bootstrap_row_prefixes]
        bootstrap_intercepts = [bootstrap_row_prefix + "_intercepts.p" for bootstrap_row_prefix in bootstrap_row_prefixes]
        bootstrap_results = [bootstrap_row_prefix + "_fit_result_df.txt" for bootstrap_row_prefix in bootstrap_row_prefixes]
        bootstrap_coefsr = [bootstrap_row_prefix + "_coefsr.p" for bootstrap_row_prefix in bootstrap_row_prefixes]
        bootstrap_resultsr = [bootstrap_row_prefix + "_fit_result_dfr.txt" for bootstrap_row_prefix in bootstrap_row_prefixes]

        bootstrap_output_dict = collections.OrderedDict()
        bootstrap_output_dict["coef"] = bootstrap_coefs
        bootstrap_output_dict["coefr"] = bootstrap_coefsr
        bootstrap_output_dict["intercept"] = bootstrap_intercepts
        # bootstrap_output_dict["interceptr"] = bootstrap_interceptsr
        # rand intercepts aren't put above because if it's a local null fit, then too many possible intercepts for each effect gene

        output_matr_df = pd.DataFrame(bootstrap_output_dict)
        output_matr_file = os.path.join(bootstrap_folder, bootstrap_outmost_name + "_output_matr_list.txt")
        output_matr_df.to_csv(output_matr_file, sep="\t", index=False)
        # print "Raw parallelilized output matrices, before integration, written to", output_matr_file




        int_matr_dict = collections.OrderedDict()
        int_matr_dict["coef"] = bootstrap_outmost_prefix + "_coefs.p"
        int_matr_dict["coefr"] = bootstrap_outmost_prefix +  "_coefsr.p"
        int_matr_dict["intercept"] = bootstrap_outmost_prefix + "_intercepts.p"
        # int_matr_dict["interceptr"] = "bootstrap" + os.sep + bootstrap_outmost_name + "_interceptsr.p"

        # append these to the list of final bootstrapped coefficients
        all_int_coefs.append(int_matr_dict["coef"])
        all_int_intercepts.append(int_matr_dict["intercept"])

        int_matr_file = bootstrap_outmost_prefix +  "_int_matr_list.txt"
        int_matr_df = pd.DataFrame(int_matr_dict, index=[0])
        int_matr_df.to_csv(int_matr_file, sep="\t", index=False)
        # print "integrated matrices written to " + int_matr_file


        bootstrap_result_dict = collections.OrderedDict()
        bootstrap_result_dict["fit_result"] = bootstrap_results
        bootstrap_result_dict["fit_resultr"] = bootstrap_resultsr



        output_df_file = bootstrap_outmost_prefix + "_output_df_list.txt"
        output_df_df = pd.DataFrame(bootstrap_result_dict)
        output_df_df.to_csv(output_df_file, sep="\t", index=False)
        # print "output dfs file written to ", output_df_file

        int_df_dict = collections.OrderedDict()
        int_df_dict["fit_result"] = bootstrap_outmost_prefix + "_fit_result_df.txt"
        int_df_dict["fit_resultr"] = bootstrap_outmost_prefix + "_fit_result_dfr.txt"

        int_df_file = bootstrap_outmost_prefix + "_int_df_list.txt"
        int_df_df = pd.DataFrame(int_df_dict, index=[0])
        int_df_df.to_csv(int_df_file, sep="\t", index=False)
        # print "Integrated dfs file written to ", int_df_file



        # just need to put all of this into the outmost name


        finish_none_script = os.path.join("bootstrap-finish-scripts", "none", "finish-none-bootstrap-" + str(b) + ".sh")
        with open(finish_none_script, 'w') as ifile:
            ifile.write("set -e\n")
            ifile.write("START=$(date)\n")
            ifile.write("time python3 integrate_outputs_rand_row.py -i " + output_matr_file + " -o " + int_matr_file +  (" -t m -a 1 " if args.only_array else " -t a "))
            ifile.write(" && " + \
                        "time python3 integrate_outputs_rand_row.py -i " + output_df_file + " -t d -o " + int_df_file + "\n"
                        )
            ifile.write("time python3 get_result_coef.py -df " + data_file + " -rdf " + rand_data_file +\
                        " -lr " + str(args.load_reps) + \
                        " -bh " + "hyper" + os.sep + "best_hyper.p" + \
                        " -o " + \
                         bootstrap_outmost_name + " -cf " +  int_matr_dict["coef"] + " -if " + int_matr_dict["intercept"] + \
                        " -cfr " + int_matr_dict["coefr"] + " -fr " + \
                        int_df_dict["fit_result"] + " -frr " + int_df_dict["fit_resultr"] + " -l " + str(args.lag) + \
                        " -sb " + "n" + " -tn " + args.test_name + " -of " + bootstrap_folder + "\n")
            ifile.write("END=$(date)\n")
            ifile.write("echo " + finish_none_script + ",$START,$END,$SECONDS >> " + finishtimefile + "\n")
            # print "Finish script, stratby None, written to", finish_none_script
            os.chmod(finish_none_script, 0o777)

        finish_none_scripts.append(finish_none_script)


        finish_effect_script = os.path.join("bootstrap-finish-scripts", "effect", "finish-effect-bootstrap-" + str(b) + ".sh")
        with open(finish_effect_script, 'w') as ifile:
            ifile.write("set -e\n")
            ifile.write("START=$(date)\n")
            ifile.write("time python3 integrate_outputs_rand_row.py -i " + output_matr_file + " -o " + int_matr_file +  (" -t m -a 1 " if args.only_array else " -t a "))
            ifile.write(" && " + \
                        "time python3 integrate_outputs_rand_row.py -i " + output_df_file + " -t d -o " + int_df_file + "\n"
                        )
            ifile.write("time python3 get_result_coef.py -df " + data_file + " -rdf " + rand_data_file +\
                        " -lr " + str(args.load_reps) + \
                        " -bh " + "hyper" + os.sep + "best_hyper.p" + \
                        " -o " + \
                        bootstrap_outmost_name + " -cf " +  int_matr_dict["coef"] + " -if " + int_matr_dict["intercept"] + \
                        " -cfr " + int_matr_dict["coefr"] + " -fr " + \
                        int_df_dict["fit_result"] + " -frr " + int_df_dict["fit_resultr"] + " -l " + str(args.lag) + \
                        " -sb " + "e" + " -tn " + args.test_name  + " -of " + bootstrap_folder + "\n")
            ifile.write("END=$(date)\n")
            ifile.write("echo " + finish_effect_script + ",$START,$END,$SECONDS >> " + finishtimefile + "\n")

            # print "Finish script, stratby effect, written to", finish_effect_script
            os.chmod(finish_effect_script, 0o777)

        finish_effect_scripts.append(finish_effect_script)


        # get all the fdr files immediately

        for fdr in fdrs:
            all_fdr_none_coefs_dict[fdr].append(os.path.join(bootstrap_folder, "fdr-" + str(fdr) + "-" + "none",
                               bootstrap_outmost_name + "-union-fdr-" + str(fdr) + "-" + "none" +  "-coefs.p"))
            all_fdr_effect_coefs_dict[fdr].append(os.path.join(bootstrap_folder, "fdr-" + str(fdr) + "-" + "effect",
                                bootstrap_outmost_name + "-union-fdr-" + str(fdr) + "-" + "effect" +  "-coefs.p"))

            all_fdr_none_intercepts_dict[fdr].append(os.path.join(bootstrap_folder, "fdr-" + str(fdr) + "-" + "none",
                               bootstrap_outmost_name + "-union-fdr-" + str(fdr) + "-" + "none" +  "-intercepts.p"))
            all_fdr_effect_intercepts_dict[fdr].append(os.path.join(bootstrap_folder, "fdr-" + str(fdr) + "-" + "effect",
                                bootstrap_outmost_name + "-union-fdr-" + str(fdr) + "-" + "effect" +  "-intercepts.p"))




        # print "-----------"


    int_coef_file = "all_bootstrap_coefs.txt"
    with open(int_coef_file, 'w') as f:
        for b_coef in all_int_coefs:
            f.write(b_coef + "\n")
    print("All integrated bootstrapped coef files written to ", int_coef_file)

    int_intercept_file = "all_bootstrap_intercepts.txt"
    with open(int_intercept_file, 'w') as f:
        for b_intercept in all_int_intercepts:
            f.write(b_intercept + "\n")
    print("All integrated bootstrapped intercept files written to ", int_intercept_file)



    all_finish_effect_script = "finish-effect-bootstrap-all.sh"
    with open(all_finish_effect_script, 'w') as f:
        f.write("set -e\n")
        for s in finish_effect_scripts:
            f.write("./" + s + "\n")
    os.chmod(all_finish_effect_script, 0o777)

    print("All bootstrap effects scripts written to ", all_finish_effect_script)


    if args.parallel_num > 0:
        print("Parallel Number (# processes per job): " + str(args.parallel_num))

        script_groups = pj.partition_inputs(finish_effect_scripts, number=int(math.ceil(len(finish_effect_scripts) * 1.0/args.parallel_num)))

        print("Number of script groups ", len(script_groups))

        parallel_scripts = []
        for i, script_group in zip(list(range(len(script_groups))), script_groups):
            appended_script_filenames = ["./" + script_filename for script_filename in script_group]
            parallel_script = " & ".join(appended_script_filenames)
            parallel_scripts.append(parallel_script)

        with open("finish-effect-bootstrap_parallel_script_list.txt", 'w') as scriptfile:
            for parallel_script in parallel_scripts:
                scriptfile.write(parallel_script + "\n")
            print("Parallel script list written to finish-effect-bootstrap_parallel_script_list.txt")



    all_finish_none_script = "finish-none-bootstrap-all.sh"
    with open(all_finish_none_script, 'w') as f:
        f.write("set -e\n")
        for s in finish_none_scripts:
            f.write("./" + s + "\n")
    os.chmod(all_finish_none_script, 0o777)

    print("All bootstrap nones scripts written to ", all_finish_none_script)


    if args.parallel_num > 0:
        print("Parallel Number (# processes per job): " + str(args.parallel_num))

        script_groups = pj.partition_inputs(finish_none_scripts, number=int(math.ceil(len(finish_none_scripts) * 1.0/args.parallel_num)))

        print("Number of script groups ", len(script_groups))

        parallel_scripts = []
        for i, script_group in zip(list(range(len(script_groups))), script_groups):
            appended_script_filenames = ["./" + script_filename for script_filename in script_group]
            parallel_script = " & ".join(appended_script_filenames)
            parallel_scripts.append(parallel_script)

        with open("finish-none-bootstrap_parallel_script_list.txt", 'w') as scriptfile:
            for parallel_script in parallel_scripts:
                scriptfile.write(parallel_script + "\n")
            print("Parallel script list written to finish-none-bootstrap_parallel_script_list.txt")





    # integrate all the bootrastrapped FDR

    bootstrap_result_folder = os.path.join("bootstrap", "bootstrap-results")
    if not os.path.exists(bootstrap_result_folder):
        os.makedirs(bootstrap_result_folder)


    bootstrap_summary_file = "get_result_bootstrap.sh"
    with open(bootstrap_summary_file, 'w') as f:
        f.write("START=$(date)\n")
        f.write("time python3 get_result_bootstrap.py -df " + data_file + " -lr " + str(args.load_reps) + \
                             " -osf " + "bootstrap" + " -rsf " + bootstrap_result_folder + " -o " + args.output_name + " -l " + str(args.lag) + " -tn " + args.test + \
                " -b " + int_coef_file + " -da 1"+ " -tbf " + "bootstrap-transpose" + " -uabrd 0\n")
        f.write("time python3 get_intercept_bootstrap.py -b " + int_intercept_file + " -rsf " + bootstrap_result_folder + " -o " + args.output_name + "\n")
        f.write("END=$(date)\n")
        f.write("echo " + bootstrap_summary_file + ",$START,$END,$SECONDS >> " + resulttimefile + "\n")
    os.chmod(bootstrap_summary_file, 0o777)
    print("Script to analyze integrated bootstrapped coefs in", bootstrap_summary_file)


    # integrate in a lite version

    bootstrap_summary_file = "get_result_bootstrap_lite.sh"
    with open(bootstrap_summary_file, 'w') as f:
        f.write("START=$(date)\n")
        f.write("time python3 get_result_bootstrap.py -df " + data_file + " -lr " + str(args.load_reps) + \
                             " -osf " + "bootstrap" + " -rsf " + bootstrap_result_folder + " -o " + args.output_name + "_lite" + " -l " + str(args.lag) + " -tn " + args.test + \
                " -b " + int_coef_file + " -da 1"+ " -dl 1 -uabrd 0\n")
        f.write("time python3 get_intercept_bootstrap.py -b " + int_intercept_file + " -rsf " + bootstrap_result_folder + " -o " + args.output_name + "\n")
        f.write("END=$(date)\n")
        f.write("echo " + bootstrap_summary_file + ",$START,$END,$SECONDS >> " + resulttimefile + "\n")
    os.chmod(bootstrap_summary_file, 0o777)
    print("Script to analyze integrated bootstrapped coefs in", bootstrap_summary_file)



    for fdr in fdrs:
        print("*************************")
        print("Integrating bootstrap files for FDR ", fdr)

        print("****EFFECT***")

        bootstrap_result_folder = os.path.join("bootstrap", "bootstrap-results-fdr-" + str(fdr) + "-effect")
        if not os.path.exists(bootstrap_result_folder):
            os.makedirs(bootstrap_result_folder)


        # write the fdr file out
        bootstrap_fdr_effect_list_file = "all_bootstrap_coefs_fdr-" + str(fdr) + "-effect.txt"
        with open(bootstrap_fdr_effect_list_file, 'w') as f:
            for b_coef in all_fdr_effect_coefs_dict[fdr]:
                f.write(b_coef + "\n")

            print("All fdr effect written to ", bootstrap_fdr_effect_list_file)


        bootstrap_fdr_effect_intercept_list_file = "all_bootstrap_intercepts_fdr-" + str(fdr) + "-effect.txt"
        with open(bootstrap_fdr_effect_intercept_list_file, 'w') as f:
            for b_intercept in all_fdr_effect_intercepts_dict[fdr]:
                f.write(b_intercept + "\n")

            print("All fdr effect written to ", bootstrap_fdr_effect_intercept_list_file)


        bootstrap_fdr_effect_summary_script = "get_result_bootstrap-fdr-" + str(fdr) + "-effect.sh"

        with open(bootstrap_fdr_effect_summary_script, 'w') as f:
            f.write("START=$(date)\n")
            f.write("set -e\n")
            f.write("time python3 get_result_bootstrap.py -df " + data_file + " -lr " + str(args.load_reps) + \
                    " -osf " + "bootstrap" + " -rsf " + bootstrap_result_folder + " -o " + args.output_name + "-fdr-" + str(fdr) + "-effect" + " -l " + str(args.lag) + " -tn " + args.test + \
                    " -b " + bootstrap_fdr_effect_list_file +  " -da 0" + " -tbf " + "bootstrap-transpose" + "-fdr-" + str(fdr) + "-effect  -uabrd 1\n")
            # f.write("time python3 get_intercept_bootstrap.py -b " + int_intercept_file + " -rsf " + bootstrap_result_folder + " -o " + args.output_name + "\n")
            f.write("END=$(date)\n")
            f.write("echo " + bootstrap_fdr_effect_summary_script + ",$START,$END,$SECONDS >> " + resulttimefile + "\n")
        os.chmod(bootstrap_fdr_effect_summary_script, 0o777)
        print("Script to analyze integrated bootstrapped coefs in", bootstrap_fdr_effect_summary_script)


        bootstrap_fdr_effect_summary_script = "get_result_bootstrap-fdr-" + str(fdr) + "-effect_lite.sh"

        with open(bootstrap_fdr_effect_summary_script, 'w') as f:
            f.write("START=$(date)\n")
            f.write("set -e\n")
            f.write("time python3 get_result_bootstrap.py -df " + data_file + " -lr " + str(args.load_reps) + \
                    " -osf " + "bootstrap" + " -rsf " + bootstrap_result_folder + " -o " + args.output_name + "_lite" + "-fdr-" + str(fdr) + "-effect"  + " -l " + str(args.lag) + " -tn " + args.test + \
                    " -b " + bootstrap_fdr_effect_list_file +  " -da 0" + " -dl 1 -uabrd 1\n")
            f.write("END=$(date)\n")
            f.write("echo " + bootstrap_fdr_effect_summary_script + ",$START,$END,$SECONDS >> " + resulttimefile + "\n")
        os.chmod(bootstrap_fdr_effect_summary_script, 0o777)
        print("Script to analyze integrated bootstrapped coefs in", bootstrap_fdr_effect_summary_script)




        print("-----------------------")


        print("****NONE***")

        bootstrap_result_folder = os.path.join("bootstrap", "bootstrap-results-fdr-" + str(fdr) + "-none")
        if not os.path.exists(bootstrap_result_folder):
            os.makedirs(bootstrap_result_folder)


        bootstrap_fdr_none_list_file = "all_bootstrap_coefs_fdr-" + str(fdr) + "-none.txt"
        with open(bootstrap_fdr_none_list_file, 'w') as f:
            for b_coef in all_fdr_none_coefs_dict[fdr]:
                f.write(b_coef + "\n")

            print("All fdr none written to ", bootstrap_fdr_none_list_file)


        bootstrap_fdr_none_intercept_list_file = "all_bootstrap_intercepts_fdr-" + str(fdr) + "-none.txt"
        with open(bootstrap_fdr_none_intercept_list_file, 'w') as f:
            for b_intercept in all_fdr_none_intercepts_dict[fdr]:
                f.write(b_intercept + "\n")

            print("All fdr none written to ", bootstrap_fdr_none_intercept_list_file)


        bootstrap_fdr_none_summary_script = "get_result_bootstrap-fdr-" + str(fdr) + "-none.sh"

        with open(bootstrap_fdr_none_summary_script, 'w') as f:
            f.write("START=$(date)\n")
            f.write("set -e\n")
            f.write("time python3 get_result_bootstrap.py -df " + data_file + " -lr " + str(args.load_reps) + \
                    " -osf " + "bootstrap" + " -rsf " + bootstrap_result_folder + " -o " + args.output_name + "-fdr-" + str(fdr) + "-none" + " -l " + str(args.lag) + " -tn " + args.test + \
                    " -b " + bootstrap_fdr_none_list_file + " -da 0" + " -tbf " + "bootstrap-transpose" + "-fdr-" + str(fdr) + "-none -uabrd 1\n")
            f.write("END=$(date)\n")
            f.write("echo " + bootstrap_fdr_none_summary_script + ",$START,$END,$SECONDS >> " + resulttimefile + "\n")
        os.chmod(bootstrap_fdr_none_summary_script, 0o777)
        print("Script to analyze integrated bootstrapped coefs in", bootstrap_fdr_none_summary_script)



        bootstrap_fdr_none_summary_script = "get_result_bootstrap-fdr-" + str(fdr) + "-none_lite.sh"

        with open(bootstrap_fdr_none_summary_script, 'w') as f:
            f.write("START=$(date)\n")
            f.write("set -e\n")
            f.write("time python3 get_result_bootstrap.py -df " + data_file + " -lr " + str(args.load_reps) + \
                    " -osf " + "bootstrap" + " -rsf " + bootstrap_result_folder + " -o " + args.output_name + "_lite" + "-fdr-" + str(fdr) + "-none" + " -l " + str(args.lag) + " -tn " + args.test + \
                    " -b " + bootstrap_fdr_none_list_file + " -da 0" + " -dl 1 -uabrd 1\n")
            f.write("END=$(date)\n")
            f.write("echo " + bootstrap_fdr_none_summary_script + ",$START,$END,$SECONDS >> " + resulttimefile + "\n")
        os.chmod(bootstrap_fdr_none_summary_script, 0o777)
        print("Script to analyze integrated bootstrapped coefs in", bootstrap_fdr_none_summary_script)



        print()
    print("FDR DONE ")
    print(" *************************************")


    print("SCRIPTS")

    with open("bootstrap_script_list.txt", 'w') as outfile:
        # lEFT OFF HERE
        for bootstrap_script in sorted(all_bootstrap_scripts):
            outfile.write("./" + bootstrap_script + "\n")
        print("bootstrap scripts written to bootstrap_script_list.txt")

        if args.parallel_num > 0:
            print("Parallel Number (# processes per job): " + str(args.parallel_num))

            script_groups = pj.partition_inputs(all_bootstrap_scripts, number=int(math.ceil(len(all_bootstrap_scripts) * 1.0/args.parallel_num)))

            print("Number of script groups ", len(script_groups))

            parallel_scripts = []
            for i, script_group in zip(list(range(len(script_groups))), script_groups):
                appended_script_filenames = ["./" + script_filename for script_filename in script_group]
                parallel_script = " & ".join(appended_script_filenames)
                parallel_scripts.append(parallel_script)

            with open("bootstrap_parallel_script_list.txt", 'w') as scriptfile:
                for parallel_script in parallel_scripts:
                    scriptfile.write(parallel_script + "\n")
                print("Parallel script list written to bootstrap_parallel_script_list.txt")


    print("TIMING")




def main():
    run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()
