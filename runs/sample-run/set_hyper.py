__author__ = 'jlu96'


import causal_pipeline as cp
import sys
import pickle
import pandas as pd
import geneTSmunging as gtm
import os
import numpy as np

def get_parser():
    # Parse arguments
    import argparse

    description = 'Given the baseline, per gene hyperparameter fit results, choose the best hyperparameter'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-tn', '--test_name', default="")

    parser.add_argument('-ind', '--int_name_dfname', required=True)

    parser.add_argument('-r', '--result_dfname', required=True, help="Save result of hyperparameters")

    parser.add_argument('-o', '--output_name', required=True, help="Save the best hyperparameter")

    parser.add_argument('-hl', '--hyper_file', required=True, help="Pickle of the original list. Make sure this has same order of the row names in integration")

    parser.add_argument('-s', '--sort_by', default='mse_avg')

    return parser


def load_and_run(args):

    if args.test_name == "":
        name = ""
    else:
        name = args.test_name.capitalize() + " "


    hyperlist = pickle.load(open(args.hyper_file, 'rb'))

    int_name_df = pd.read_csv(args.int_name_dfname, sep="\t")

    print("Loading integrated")
    print(int_name_df.head())

    hyper_fit_dfs = [pd.read_csv(int_name_df[x].values[0], sep="\t")
                     if os.path.exists(int_name_df[x].values[0])  else None
                     for x in int_name_df]

    # Remove the Nones for which there is no information.
    remove_list = []
    for i in range(len(hyper_fit_dfs[:])):
        try:
            # Check if its empty
            if hyper_fit_dfs[i].empty:
                remove_list.append(i)
            # If it's equal to None will have an AttributeError here
        except AttributeError:
            remove_list.append(i)

    hyper_fit_dfs = [h for i, h in enumerate(hyper_fit_dfs) if i not in remove_list]
    hyperlist = [h for i, h in enumerate(hyperlist) if i not in remove_list]



    # Get the best hyper
    hyper_df = cp.summarize_hyper_fit_dfs(hyper_fit_dfs, hyperlist)

    best_hyper, best, hyper_df = cp.get_best_hyper(hyper_df, sort_by=args.sort_by)



    # Write the hypers out
    pickle.dump(best_hyper, open(args.output_name, 'wb'))
    hyper_df.to_csv(args.result_dfname, sep="\t", index=0)

    print("Test is ", name)
    print("Best hyper is ", best_hyper)
    print("Best hyper result is ", best)

    print("Best hyper written to ", args.output_name)
    print("Hyper result written to ",  args.result_dfname)


    if not os.path.exists("hyper"):
        os.makedirs("hyper")

    # Get correlations
    mse_vec = np.array([np.array(hyper_fit_df["mse"].values) for hyper_fit_df in hyper_fit_dfs])

    print(mse_vec.shape)

    mse_corr = np.corrcoef(mse_vec)
    gtm.save_gene_matrix("hyper" + os.sep + "mse_corr.txt", mse_corr, hyperlist)
    print("MSE Correlation:")
    print(mse_corr)
    print("MSE corr. matrix saved to ", "hyper" + os.sep + "mse_corr.txt")

    r2_vec = np.array([hyper_fit_df["r2"].values for hyper_fit_df in hyper_fit_dfs])
    r2_corr = np.corrcoef(r2_vec)
    gtm.save_gene_matrix("hyper" + os.sep + "r2_corr.txt", r2_corr, hyperlist)
    print("R2 Correlation")
    print(r2_corr)
    print("R^2 corr. matrix saved to ", "hyper" + os.sep + "r2_corr.txt")


    # Plot the hyperparameters
    if not os.path.exists("plots"):
        os.makedirs("plots")
    if not os.path.exists("plots" + os.sep + "hyper"):
        os.makedirs("plots" + os.sep + "hyper")

    cp.plot_corr_matrix(mse_corr, cp.hyperlist_to_labellist(hyperlist), title="MSE correlation among " + name + "hyperparams", filename="plots" + os.sep + "hyper" + os.sep + "mse_corr")
    cp.plot_corr_matrix(r2_corr, cp.hyperlist_to_labellist(hyperlist), title="$r^2$ correlation among " + name + "hyperparams", filename="plots" + os.sep + "hyper" + os.sep + "r2_corr")




    cp.plot_hyper_boxplot(cp.hyperlist_to_labellist(hyperlist), hyper_fit_dfs, "r2", xlabel=name + "Hyperparameter", ylabel="$r^2$", title=name + "Hyperparameter VS $r^2$", filename="plots" + os.sep + "hyper"+ os.sep + "hyperVSr2",
                         hyper_color_labels=[(cp.hyper_to_label(best_hyper), "k", "Best: " + cp.hyper_to_label(best_hyper) + ", $r^2$ = " + str(np.round(best["r2_avg"].values[0], 1)) )],
                         horizontal_line_color_labels=[(best["r2_avg"].values[0], 'k', None)])
    cp.plot_hyper_boxplot(cp.hyperlist_to_labellist(hyperlist), hyper_fit_dfs, "mse", xlabel=name + "Hyperparameter", ylabel="Mean-Squared Error", title=name + "Hyperparameter VS MSE", filename="plots" + os.sep + "hyper" + os.sep + "hyperVSmse",
                         hyper_color_labels=[(cp.hyper_to_label(best_hyper), "k", "Best: " + cp.hyper_to_label(best_hyper) + ", MSE = " + str(np.round(best["mse_avg"].values[0], 1)))],
                         horizontal_line_color_labels=[(best["mse_avg"].values[0], 'k', None)])
    cp.plot_hyper_boxplot(cp.hyperlist_to_labellist(hyperlist), hyper_fit_dfs, "avg_df", xlabel=name + "Hyperparameter", ylabel="Degrees of Freedom", title=name + "Hyperparameter VS df", filename="plots" + os.sep + "hyper" + os.sep + "hyperVSdof",
                         hyper_color_labels=[(cp.hyper_to_label(best_hyper), "k", "Best: " + cp.hyper_to_label(best_hyper) + ", df = " + str(int(np.round(best["df_avg"].values[0]))))],
                         horizontal_line_color_labels=[(best["df_avg"].values[0], 'k', None)])


    print("Correlation between hyperparameter results", "plots" + os.sep + "hyper")
    print("Hyper box plots of r^2, mse, avg d.o.f. written to  ", "plots" + os.sep + "hyper")



def main():
    load_and_run(get_parser().parse_args(sys.argv[1:]))




if __name__ == '__main__':
    main()