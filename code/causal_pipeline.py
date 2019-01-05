__author__ = 'jlu96'


import numpy as np
from sklearn.model_selection import LeaveOneOut
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
reload(plt)
import pandas as pd
import pickle
import collections
import fit_methods as fm
import geneTSmunging as gtm
import copy


global test2fit_method
test2fit_method = {'e': fm.fit_enet, "l": fm.fit_lasso, "r": fm.fit_ridge}

global args2stratify_by
args2stratify_by = {"e": "effect", "n": "none"}

def cross_validate(X_matr, lag, fit_method, hyperlist, rows=None, has_reps=False, **kwargs):
    """
    X_matr: n x T matrix of genes
    lag
    fit_method
    hyperlist: list of settings of "hyper"
    rows: optional restriction to a certain list of response variables

    For each hyperparam setting, use that hyperparam in fitting the whole data and evaluating prediction error
    Take the hyperparam with least avg mse.
    Return results from each hyperparam setting and the best hyperparam

    Return: best_hyper, best, hyper_df

    """

    n = X_matr.shape[0]
    T = X_matr.shape[1]

    if rows == None:
        rows = range(n)

    hyper_fit_dfs = []

    for hyper in hyperlist:

        fit_results = []

        loo = LeaveOneOut().split(X_matr)

        # Iterate over all the genes
        for train_index, test_index in loo:

            # Only test the responses in the appropriate row
            if test_index in rows:
                _, X_test = X_matr[train_index], X_matr[test_index]

                # there is a 3rd axis
                if has_reps:
                    Y_test = np.reshape(X_test, (1, T, X_test.shape[2]))
                else:
                    Y_test = np.reshape(X_test, (1, T))


                fit_result =  collections.OrderedDict()
                fit_result["hyper"] = hyper
                fit_result["row"] = test_index[0]

                # Since Y_test is not in X_train, replace_rows is None
                fit_result.update(fm.perform_loto_cv(X_matr=X_matr, Y_matr=Y_test, lag=lag, fit_method=fit_method,
                                                                             hyper=hyper, replace_row=test_index,
                                                                            has_reps=has_reps, **kwargs))

                fit_results.append(fit_result)

        fit_result_df = pd.DataFrame(fit_results)

        hyper_fit_dfs.append(fit_result_df)

        print "Hyper: ", hyper
        print fit_result_df.head(n=20)

    return hyper_fit_dfs


def summarize_hyper_fit_dfs(hyper_fit_dfs, hyperlist):
    """
    :param hyper_fit_dfs:  Result of fit per ttimepoint
    :param hyperlist: Hyperparameters
    :return: hyper_df: Summary of overall performance for each hyper_df
    """


    total_sselist = []
    avg_dflist = []
    std_dflist =[]
    avg_r2list = []
    std_r2list = []
    avg_nlist = []
    avg_mselist = []
    std_mselist = []

    for fit_result_df in hyper_fit_dfs:
        stats_df = fit_result_df[["sse", "avg_df", "r2", "n", "mse"]]

        total_df = stats_df.sum()

        avg_df = total_df / stats_df.shape[0]

        std_df = stats_df.std()

        total_sselist.append(total_df["sse"])
        avg_dflist.append(avg_df["avg_df"])
        std_dflist.append(std_df["avg_df"])
        avg_r2list.append(avg_df["r2"])
        std_r2list.append(std_df["r2"])
        avg_nlist.append(avg_df["n"])
        avg_mselist.append(avg_df["mse"])
        std_mselist.append(std_df["mse"])

    summary_dict = collections.OrderedDict()
    summary_dict["hyper"] = hyperlist
    summary_dict["n_avg"] = avg_nlist
    summary_dict["mse_avg"] = avg_mselist
    summary_dict["mse_std"] = std_mselist
    summary_dict["df_avg"] = avg_dflist
    summary_dict["df_std"] = std_dflist
    summary_dict["r2_avg"] = avg_r2list
    summary_dict["r2_std"] = std_r2list
    summary_dict["sse_total"] = total_sselist
    hyper_df = pd.DataFrame(summary_dict)

    return hyper_df


def get_best_hyper(hyper_df, sort_by="mse_avg", ascending=True):
    """
    :param hyper_df:
    :param sort_by:
    :param ascending:
    :return: the best hyper params
    """

    hyper_df.sort_values(sort_by, inplace=True, ascending=ascending)

    best = hyper_df.head(n=1)

    best_hyper = best["hyper"].values[0]

    return best_hyper, best, hyper_df




def run_cross_validate(geneTS, fit_method=fm.fit_lasso,
                        hyperlist=10**(-1 * np.arange(0, 4, 1.0)),
                        lag=2,
                       rows=None,
                       sort_by="mse_avg",
                        save_prefix=None,
                       has_reps=False):
    """
    Writes out the results for a given hyper-parameter list
    """


    print "Hyper-parameters for cross-validation"
    print hyperlist

    # Cross-validate

    hyper_fit_dfs = cross_validate(geneTS, lag, fit_method, hyperlist,
                                                               rows=rows,
                                   has_reps=has_reps)

    hyper_df = summarize_hyper_fit_dfs(hyper_fit_dfs, hyperlist)

    best_hyper, best, hyper_df = get_best_hyper(hyper_df, sort_by=sort_by)

    print "Hypers results:"
    print hyper_df


    return best_hyper, best, hyper_df, hyper_fit_dfs



def float_to_label(value):
    assert isinstance(value, float) or isinstance(value, int)
    return "%.0E" % value

def tuple_to_label(value):
    assert isinstance(value, tuple)
    return "(" + ", ".join([float_to_label(x) for x in value]) + ")"


def hyper_to_label(hyper):
    """
    :param hyper: hyperparameter
    :return: the corresponding label
    """
    if isinstance(hyper, float) or isinstance(hyper, int):
        return float_to_label(hyper)
    elif isinstance(hyper, tuple):
        return tuple_to_label(hyper)


def hyperlist_to_labellist(hyperlist):
    """
    :param hyperlist:
    :return: labellist, labels to use for plotting
    """
    return [hyper_to_label(hyper) for hyper in hyperlist]


def run_to_label(run):
    """
    :param run: string containing parameters
    :return: string designed for labels
    """

    label = run[:]

    replacedict = collections.OrderedDict()
    replacedict["original"] =  "orig"
    replacedict["integration"] = "int"
    replacedict["0mean-unnormalized"] = "ZU"
    replacedict["0mean-1var"] = "ZS"

    for key in replacedict:
        label = label.replace(key, replacedict[key])

    return label


#
# def hyperlist_to_namelist(hyperlist):
#     """
#     :param hyperlist:
#     :return: labellist, labels to use for plotting
#     """
#     hyper_value = hyperlist[0]
#
#     if isinstance(hyper_value, float) or isinstance(hyper_value, int):
#         return hyperlist_to_labellist(hyperlist)
#
#     elif isinstance(hyper_value, tuple):
#         return [hyper.replace(" ", "") for hyper in hyperlist_to_labellist(hyperlist)]




def plot_hyper_boxplot(hyperlist, fit_result_dfs, fit_result_key, xlabel="Hyperparameters", ylabel="Output parameter",
                     title="Hyperparameters VS Output parameters",filename=None, horizontal_line_color_labels=None, hyper_color_labels=None,
                       hyper_color_label_margin=0.4):
    """
    hyperlist: the list of hyperparameters
    fit_result_dfs: each a dataframe, assume ordered as hyperlist
    fit_result_key: the column name of fit_result_df to plot

    Plots the boxplot with 1.5 IQR
    """
    assert len(hyperlist) == len(fit_result_dfs)

    hyper_value = hyperlist[0]

    try:
        # Get length of the hyperparam
        label_length = len(hyper_value)

        plot_height = 5 + label_length * 0.2
    except TypeError:
        plot_height = 5

    if len(hyperlist) > 10:
        figsize = (len(hyperlist), plot_height)
    else:
        figsize = (8, plot_height)


    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=figsize)

    data = [fit_result_df[fit_result_key].values for fit_result_df in fit_result_dfs]
    pos = range(1, len(hyperlist) + 1)

    plt.boxplot(data, positions=pos, showmeans=True)

    axes.yaxis.grid(True)


#     avg = np.average(data, axis=1)
#     std = np.std(data, axis=1)
#     df = np.array([len(fit_result_df) - 1  for fit_result_df in fit_result_dfs ])
#     min_val = np.min(data)

#     plt.errorbar(pos, avg, yerr=stats.t.ppf(0.95, df) * std, color='r')



    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel(ylabel,fontsize=20)
    plt.title(title, fontsize=20)
    plt.xticks(pos, hyperlist)
    plt.yticks(fontsize=15)

    labels = axes.get_xticklabels()
    plt.setp(labels, rotation=90, fontsize=15)

    if hyper_color_labels != None:
        for hyper, color, label in hyper_color_labels:
            assert hyper in hyperlist
            loc = pos[hyperlist.index(hyper)]
            left_loc = loc - hyper_color_label_margin
            right_loc = loc + hyper_color_label_margin


            plt.axvline(left_loc, color=color, label=label)
            plt.axvline(right_loc, color=color)

            # plt.axvline(left_loc, color=color, label=label, linestyle='dashed')
            # plt.axvline(right_loc, color=color, linestyle='dashed')



    if horizontal_line_color_labels !=None:
        for line, color, label in horizontal_line_color_labels:
            plt.axhline(line, color=color, label=label)

            # plt.axhline(line, color=color,label=label, linestyle='dashed')




    plt.legend(loc='best')
    plt.tight_layout()

    if filename:
        fig.savefig(filename)
        print "Plot saved to ", filename

    plt.show()
    plt.close()




def plot_corr_matrix(corr_matr, labels, ylabels=None, title='Correlation matrix', cmap=copy.deepcopy(plt.cm.Blues),
                    xlabel="Hyperparameters", ylabel="Hyperparameters", filename=None,
                     min_thresh=0, max_thresh=1, figsize=None,
                     change_bad_col=True,
                     bad_col=(1, 0.7, 0.4, 1),
                     change_over_col=False,
                     over_col='k',
                     xticksize=12,
                     yticksize=12,
                     titlesize=20,
                     xlabelsize=20,
                     ylabelsize=20,
                     rotate_x=True,
                     rotate_x_by=90,
                     colorbar_label="",
                     colorbar_ticksize=12,
                     colorbar_labelsize=12,
                     change_x_pad=False,
                     x_pad=10,
                     change_y_pad=False,
                     y_pad=10,
                     do_special_marker=False,
                     special_marker_rows=[],
                     special_marker_columns=[],
                     special_marker_color='k',
                     special_marker_type='*',
                     special_marker_size=80,
                     ):
    """
    :param corr_matr:
    :param labels: For the x-axis. Can also be for y-axis
    :param ylabels: For y-axis
    :param title:
    :param cmap:
    :param xlabel:
    :param ylabel:
    :param filename:
    :param min_thresh:
    :param max_thresh:
    :return:
    """
    if len(labels) != corr_matr.shape[1]:
        raise ValueError("X Labels must have same shape as # corr matrix cols")
    try:
        len(ylabels)
    except TypeError:
        if ylabels == None:
            ylabels = labels
    if len(ylabels) != corr_matr.shape[0]:
        raise ValueError("Y labels must have same shape as # corr matrix rows")


    if figsize == None:
        if corr_matr.shape[0] < 16:
            figsize = ((8,8))
        else:
            figsize = ((corr_matr.shape[1]/2), (corr_matr.shape[0]/2))
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=figsize)

    if change_bad_col:
        cmap.set_bad(bad_col)
    if change_over_col:
        cmap.set_over(over_col)


    if do_special_marker:
        print "Doing special marker. Be wary of wrong turns."
        assert len(special_marker_rows) == len(special_marker_columns)
        plt.scatter(special_marker_columns, special_marker_rows, color=special_marker_color,
                    s=special_marker_size, marker=special_marker_type)
        plt.xlim((-0.5, len(labels) - 0.5))
        plt.ylim((-0.5, len(ylabels) -0.5))
        plt.gca().invert_yaxis()

    plt.imshow(corr_matr, interpolation='nearest', cmap=cmap)
    cb = plt.colorbar(extend='both')
    plt.clim(min_thresh, max_thresh)
    cb.ax.tick_params(labelsize=colorbar_ticksize)
    cb.set_label(colorbar_label, rotation=90, fontsize=colorbar_labelsize)
    #
    # axes.set_yticklabels(ylabels, ha='right')
    # axes.set_xticklabels(labels, va='top')

    # if do_special_marker:
    #     print "Doing special marker. Be wary of wrong turns."
    #     assert len(special_marker_rows) == len(special_marker_columns)
    #     plt.scatter(special_marker_columns, special_marker_rows, color=special_marker_color,
    #                 s=special_marker_size, marker=special_marker_type)
    #     plt.xlim((-0.5, len(labels) - 0.5))
    #     plt.ylim((-0.5, len(ylabels)-0.5))
    #     plt.imshow(corr_matr, interpolation='nearest', cmap=cmap)

    xtick_marks = np.arange(len(labels))
    plt.xticks(xtick_marks, labels, fontsize=xticksize, va='top')
    if rotate_x:
        plabels = axes.get_xticklabels()
        plt.setp(plabels, rotation=rotate_x_by, fontsize=xticksize, va='top', ha='center')
    ytick_marks = np.arange(len(ylabels))
    plt.yticks(ytick_marks, ylabels, fontsize=yticksize, ha='right', va='center')
    plt.title(title, fontsize=titlesize)


    plt.xlabel(xlabel, fontsize=xlabelsize)
    plt.ylabel(ylabel, fontsize=ylabelsize)

#     position = fig.add_axes([1, 0.01, 0.02, 0.85])
#     # cax, kw = mpl.colorbar.make_axes(fig_ax)

#     mpl.colorbar.colorbar_factory(position, im)

    plt.tight_layout()
    if change_x_pad:
        axes.tick_params(axis='x', which='major', pad=x_pad)
    if change_y_pad:
        axes.tick_params(axis='y', which='major', pad=y_pad)
    if filename:
        fig.savefig(filename, bbox_inches="tight", pad_inches=0.5)
        print "Plot saved to ", filename
    plt.show()
    plt.close()


def plot_coef(df, cause_gene, effect_gene, lag, coef, savefile=True, file_prefix="", title=None,
         cause_label=None, effect_label=None, legend_fontsize=None, ylabel=None, is_ensg=False,
              file_type=".pdf", **kwargs):
    if title == None:
        title = "Cause-Effect pair: (" + cause_gene + ", " + effect_gene + ")"

    if savefile == True:
        filename = file_prefix + cause_gene.replace(".", ",") + "-" + effect_gene.replace(".", ",") + "-lag-" + str(int(lag)) + "-coef-" + float_to_label(coef) + file_type
        print "Filename: ", filename
    else:
        filename = None

    # if is_ensg:
    #     label_spacing = 20
    # else:
    #     label_spacing = 10

    pco = np.round(coef, 2)

    if cause_label == None:
        cause_label = "Cause:   " + cause_gene+ " (Coef = " + str(pco) + ", Lag = " + str(lag) + ")"
        # cause_label = "Cause: " + ("{0:>" + str(label_spacing) + "}").format(cause_gene) + " (Coef = " + str(pco) + ", Lag = " + str(lag) + ")"
    if effect_label == None:
        effect_label = "Effect:  " + effect_gene

        # effect_label = "Effect: " + ("{0:>" + str(label_spacing) + "}").format( effect_gene)

    gtm.plot_genes(df, [cause_gene, effect_gene], title=title, filename=filename, gene_labels=[cause_label, effect_label], plot_outside=False,
                  legend_fontsize=legend_fontsize, ylabel=ylabel, **kwargs)



def plot_all_coef(acoefs, df, genes, min_coef=0.01, file_prefix=None, savefile=True, verbose=False, **kwargs):
    """
        acoefs: lag x n x n matrix, causes are rows, effets are columns
        df: original gene df to plot
        genes: list of genes
        min_coef: min entry size to be plotted
    """

    if verbose:
        print "Min Coef to plot: ", min_coef

    # iterate down the coefs of each effect gene. This is one effect gene per column.
    for j in range(acoefs.shape[2]):
        out_gene = genes[j]


        # acoef is a  lag x n causes matrix
        acoef = acoefs[:, :, j]

        # preds is a lag x n index matrix of the good coefficients
        preds = np.where(np.absolute(acoef) > min_coef)

        pcoefs = acoef[preds]

        cause_num = len(preds[0])

        if verbose:
            if cause_num > 0:
                print "Out gene: ", out_gene
                print "Causes: ", cause_num

        for i in range(cause_num):

            # the lag
            lag = preds[0][i] + 1
            cause_gene = genes[preds[1][i]]
            effect_gene = out_gene
            coef = pcoefs[i]

            plot_coef(df, cause_gene, effect_gene, lag, coef, savefile=savefile, file_prefix=file_prefix, **kwargs)





        # # the lags are encoded in the first index of acoefs
        # lags = preds[0] + 1
        # cgenes = genes[preds[1]]
        # cos = np.round(acoef, 2)[preds]
        #
        # if verbose:
        #     print "Out gene: ", out_gene
        #
        # for i in range(len(preds)):
        #
        #     plag = lags[i]
        #     cause_gene = cgenes[i]
        #     effect_gene = out_gene
        #     coef = cos[i]
        #
        #     plot_coef(df, cause_gene, effect_gene, plag, coef, savefile=savefile, file_prefix=file_prefix, **kwargs)





# def plot_all_coefs(acoefs, df, genes, lag, min_coef=0.01, num_per_plot=3, file_prefix=None, title_prefix="Causal genes for ", verbose=False,
#                        **kwargs):
#         """
#         acoefs: aligned coefs, of form (lag x n x n), acoefs[i] is lag_i+1
#         df: Original TS
#         genes: list of all genes
#         """
#
#         # iterate down the coefs of each effect gene. This is one effect gene per column.
#         for j in range(acoefs.shape[2]):
#             out_gene = genes[j]
#
#             acoef = acoefs[:, :, j]
#             preds = np.where(np.absolute(acoef) > min_coef)
#
#             lags = preds[0] + 1
#             cgenes = genes[preds[1]]
#             cos = np.round(acoef, 2)[preds]
#
#             num_plots = int(np.ceil(len(cgenes) * 1.0 / num_per_plot))
#
#             if verbose:
#                 print "Out gene: ", out_gene
#
#             for i in range(num_plots):
#                 plags = lags[i * len(cgenes)/num_plots: (i+ 1) * len(cgenes)/num_plots]
#                 pgenes = cgenes[i * len(cgenes)/num_plots: (i+ 1) * len(cgenes)/num_plots]
#                 pcos = cos[i * len(cgenes)/num_plots: (i+ 1) * len(cgenes)/num_plots]
#                 labels = ["{0:>20}".format(out_gene + ":") + " Coef, Lag" ] + ["{0:>20}".format(pgene + ":") + " " + str(pco) + ", " + str(plag) for pgene, pco, plag in zip(pgenes, pcos, plags)]
#
#                 if verbose:
#                     print "Part: ", i + 1
#                     print "Lag points: ", plags
#                     print "Pred genes:", pgenes
#                     print "Pred coefs: ", pcos
#                     print "Labels are ", labels
#
#                 plot_genes(df, [out_gene] + list(pgenes), title=title_prefix + out_gene + " , Part " + str(i+1),
#                               filename=None if file_prefix == None else file_prefix + out_gene.replace(".", ",") + "_lag-" + str(lag) + \
#                                                                         "-" + "-".join([x.replace(".", ",") for x in list(pgenes)]),
#                               gene_labels=labels, plot_outside=False,
#                            **kwargs)







def fit_all(X_matr, Y_matr, rows, lag, fit_method, save_prefix=None, save_XY=True, verbose=False,
            has_reps=False, bootstrap=False, seed=None, only_array=False, **kwargs):
    """
    X_matr: m x T of X's
    Y_matr: n x T of Y's
    rows: rows of X to replace in Y. There should be as many rows to replace in X_matr as there are in Y_matr
    lag:
    fit_method:
    bootstrap:
    seed:
    only_array: If True, return coefs as a matrix instead

    Perform the fit for each Y_matr
    Save each X_t and Y_t?


    Return: coefs: (m*(lag), m) with the columns n (corresponding w/ rows in X_matr) filled in, intercepts: (1, m), fit_result_df: (n, num_results)
    """

    assert len(rows) == Y_matr.shape[0]
    assert X_matr.shape[1] == Y_matr.shape[1]

    if has_reps:
        assert X_matr.shape[2] == Y_matr.shape[2]


    # m: number of predictor genes
    # n: number of response gene
    m = X_matr.shape[0]
    n = Y_matr.shape[0]
    T = X_matr.shape[1]

    coefs = np.zeros((m*lag, m))
    intercepts = np.zeros((1, m))

    fit_results = []

    for i, row in zip(range(Y_matr.shape[0]), rows):
        if has_reps:
            Y = np.reshape(Y_matr[i,], (1, T, Y_matr.shape[2]))

        else:
            Y = np.reshape(Y_matr[i,], (1, T))


        fit_result_dict = collections.OrderedDict()
        fit_result_dict["row"] = row
        if "hyper" in kwargs:
            fit_result_dict["hyper"] = kwargs["hyper"]


        X_t, Y_t, Y_pred, coef, intercept, fit_result = fm.perform_test(X_matr=X_matr, Y_matr=Y,
                                                                    lag=lag, fit_method=fit_method,
                                                                    replace_row=row,
                                                                        has_reps=has_reps,
                                                                     bootstrap=bootstrap,
                                                                     seed=seed,
                                                                        **kwargs)
        fit_result_dict.update(fit_result)

        # These are cause by effect
        coefs[:, row] = coef.flatten()
        intercepts[:, row] = intercept
        fit_results.append(fit_result_dict)




        if verbose:
            print i, row
            print "X: ", X_matr
            print "Y: ", Y
            print "X_t: ", X_t
            print "Y_t: ", Y_t
            print "Y_pred: ", Y_pred
            print "Checking Y_pred: ", fm.compute_fit(X_t, Y_t, coef, intercept)
            print "coef: ", coef
            print "intercept: ", intercept
            print "fit result: ", fit_result

    fit_result_df = pd.DataFrame(fit_results)

    print "Finished fitting all"


    # remember coefs are cause-by-effect.

    if only_array:
        return coefs[:, rows], intercepts[:, rows], fit_result_df

    else:
        return coefs, intercepts, fit_result_df



def fit_all_random(X_matr, rand_X_matr, Y_matr, rows, lag, fit_method, save_prefix=None, save_XY=True, verbose=False,
            has_reps=False, bootstrap=False, seed=None, only_array=False, **kwargs):
    """
    X_matr: m x T (x r) of X's
    rand_X_matr: m x T ( x r) of X's. Same as X_matr except each row (within each replicate) is permuted independently and randomly
    Y_matr: n x T (x r) of Y's
    rows: rows of X to replace in Y. There should be as many rows to replace in X_matr as there are in Y_matr
    lag:
    fit_method:
    only_array: just return the relevant rows of the array

    Perform the fit for each Y_matr
    Save each X_t and Y_t?


    Return: coefs: (m*(lag), m) with the columns n (corresponding w/ rows in X_matr) filled in, intercepts: (1, m), fit_result_df: (n, num_results)
    """

    assert len(rows) == Y_matr.shape[0]
    assert X_matr.shape[1] == Y_matr.shape[1]
    assert rand_X_matr.shape[1] == X_matr.shape[1]

    if has_reps:
        assert X_matr.shape[2] == Y_matr.shape[2]
        assert rand_X_matr.shape[2] == X_matr.shape[2]


    m = X_matr.shape[0]
    n = Y_matr.shape[0]
    T = X_matr.shape[1]

    coefs = np.zeros((m*lag, m))

    # We keep no intercepts
    # intercepts = np.zeros((1, m, m))

    fit_results = []

    for i, row in zip(range(n), rows):
        if has_reps:
            Y = np.reshape(Y_matr[i,], (1, T, Y_matr.shape[2]))

        else:
            Y = np.reshape(Y_matr[i,], (1, T))


        fit_result_dict = collections.OrderedDict()
        fit_result_dict["row"] = row
        if "hyper" in kwargs:
            fit_result_dict["hyper"] = kwargs["hyper"]




        coef, fit_result, coef_temp, coef_temps, intercept_temps = fm.perform_test_random(X_matr=X_matr, rand_X_matr=rand_X_matr,
                                                                        Y_matr=Y,
                                                                    lag=lag, fit_method=fit_method,
                                                                    replace_row=row,
                                                                        has_reps=has_reps,
                                                                            bootstrap=bootstrap, seed=seed,
                                                                        **kwargs)
        fit_result_dict.update(fit_result)

        # These are cause by effect
        coefs[:, row] = coef.flatten()

        fit_results.append(fit_result_dict)



        if verbose:
            print i, row
            print "X: ", X_matr
            print "Y: ", Y
#             print "X_t: ", X_t
#             print "Y_t: ", Y_t
#             print "Y_pred: ", Y_pred
#             print "Checking Y_pred: ", fm.compute_fit(X_t, Y_t, coef, intercept)
            print "coef: ", coef
            print "fit result: ", fit_result

    fit_result_df = pd.DataFrame(fit_results)

    if only_array:
        return coefs[:, rows], fit_result_df

    else:
        return coefs, fit_result_df


def summarize_fit(coefs, intercepts, fit_result_df, filename=None, hyper=None, test_name=None,
                  lag=None):

    summary_string = ""

    if test_name != None:
        summary_string += "Test: " + test_name + "\n"

    if lag != None:
        summary_string += "Lag: " + str(lag) + "\n"

    if hyper != None:
        summary_string += "Hyper: " + str(hyper) + "\n"

    summary_string += "Responses: " + str(coefs.shape[1]) + "\n"
    summary_string += "Coefs per response:" + str(coefs.shape[0]) + "\n"
    summary_string += "Total nonzero:" + str(len(np.nonzero(coefs)[0])) + "\n"
    summary_string += "Nonzero Coefs/Response:" + str(len(np.nonzero(coefs)[0]) * 1.0/ coefs.shape[1]) + "\n"
    summary_string += "Intercept avg: " + str(np.average(intercepts)) + "std: " + str(np.std(intercepts)) + "\n"
    summary_string += "Average fit result: " + "\n"
    summary_string += str(fit_result_df.mean()) + "\n"




    print summary_string

    if filename != None:
        with open(filename, 'w') as ifile:
            ifile.write(summary_string)

def get_union_network(dfs, suffixes, on_cols= ["Cause-Effect", "Cause", "Effect", "Test", "Lag"],  how="outer", on_to_front=True):
    """
    :param dfs: Network dfs. Each row is an edge
    :param on_cols: columns to join on
    :param suffixes: suffixes for the dfs
    :param how: how to join
    :param on_to_front: move the on columns to front
    :return: union_df
    """
    assert len(dfs) >= 2
    assert len(suffixes) == len(dfs)

    # Do the first union separately

    union_df = pd.merge(dfs[0], dfs[1], on=on_cols, how=how, suffixes=[suffixes[0], suffixes[1]])
    if on_to_front:
        cols = union_df.columns.values
        back_cols = [c for c in cols if c not in on_cols]
        union_df = union_df.reindex(columns=on_cols + back_cols)

    if len(dfs) == 2:
        return union_df
    else:
        i = 2

        while i < len(dfs):
            union_df = pd.merge(union_df, dfs[i], on=on_cols, how=how, suffixes=["", suffixes[i]])

            if on_to_front:
                cols = union_df.columns.values
                back_cols = [c for c in cols if c not in on_cols]
                union_df = union_df.reindex(columns=on_cols + back_cols)

            i += 1

        return union_df


# def get_max_network_old(dfs, max_col, index_col, other_cols=None):
#     """
#     :param dfs: Network dfs. Each row is an edge
#     :param max_col: Column to take the maximum over
#     :param index_col: Index to use
#     :return: max_df: The df of unioned indices where each max_col is its max value
#     """
#         # check that all networks have same column names
#     for df in dfs:
#         assert (df.columns.values == dfs[0].columns.values).all()
#
#
#     if other_cols == None:
#
#         cols = dfs[0].columns.values
#
#         # get the other columns
#
#         other_cols = [col for col in cols if col != max_col]
#
#     # Set the same indices:
#     for df in dfs:
#         df.index = df[index_col].values
#
#     indices = np.unique(np.concatenate(tuple(df.index.values for df in dfs)))
#
#
#     # tuples of form (df #, df max col value)
#     df_maxes = [max(zip(range(len(dfs)), [df[df[index_col] == index][max_col].values[0] if index in df[index_col].values
#                                           else np.NINF for df in dfs]),
#                                     key=lambda entry: entry[1]) for index in indices]
#
#
#     max_df = pd.DataFrame({index_col: indices})
#     max_df[max_col] = zip(*df_maxes)[1]
#     for other_col in other_cols:
#         max_df[other_col] = [dfs[df_max[0]][dfs[df_max[0]][index_col] == index][other_col].values[0]
#                              for df_max, index in zip(df_maxes, indices)]
#
#
#     return max_df



def get_max_network(dfs, max_col, index_col):
    """
    :param dfs: Network dfs. Each row is an edge
    :param max_col: Column to take the maximum over
    :param index_col: Index to use
    :return: max_df: The df of unioned indices where each max_col is its max value and each row comes from the max value df
    """
        # check that all networks have same column names
    for df in dfs:
        if not (set(list(df.columns.values)) == set(list(dfs[0].columns.values))):
            raise ValueError("COLUMNS DONT MATCH: " + str(df.columns.values) + " " +
                             str(dfs[0].columns.values))


    whole_dfs = []

    for df in dfs:
        whole_df = df.set_index([index_col],  drop=False)

        whole_dfs.append(whole_df)

    indices = np.unique(np.concatenate(tuple(whole_df.index.values for whole_df in whole_dfs)))
    for i in range(len(whole_dfs)):
        whole_df = whole_dfs[i]

        missing_indices = np.setdiff1d(indices, whole_df.index)
        missing_df = pd.DataFrame(data=dict(zip(whole_df.columns.values, [np.NINF for x in whole_df.columns.values])),
                                            index=missing_indices)
        whole_dfs[i] = pd.concat((whole_df, missing_df)).sort_index()
        # at this point all dfs should have same ordered indices

    df_maxes = np.argmax(np.array([whole_df[max_col].values for whole_df in whole_dfs]), axis=0)
    #print df_maxes
    df_indices = [df_maxes == i for i in range(len(whole_dfs))]


    max_df = pd.concat(tuple(whole_df[df_index] for df_index, whole_df in zip(df_indices, whole_dfs)))

    return max_df





def run(geneTS, geneTSr, hyper, fit_method=fm.fit_lasso,
        lag=2,
        rows=None,
        save_prefix=None,
        has_reps=False,
        null="l",
        only_array=False
    ):
    """
    geneTS: n x T original gene matrix
    geneTSr: n x T randomized gene matrix where each independently shuffled across time
    hyper: hyper-parameter
    fit_method: fit_method to use
    lag
    rows: geneTS will be the response variable here. This determines which rows of geneTS to select for response variables
        Defaulting None will just use all response variables as y
    save_prefix: how to save
    has_reps: if data has replicates
    null: "l" if to use the local, "g" if to use the global
    """

    if rows == None:
        rows = range(geneTS.shape[0])

    if null not in {"l", "g"}:
        raise ValueError("Null must either be l (local) or g (global)")

    geneTSy = geneTS[rows,:]

    print "Y shape: ", geneTSy.shape

    coefs, intercepts, fit_result_df = fit_all(geneTS, geneTSy, rows, lag, fit_method, hyper=hyper, has_reps=has_reps,
                                               only_array=only_array)


    # use local null
    if null == "l":
        coefsr, fit_result_dfr = fit_all_random(geneTS, geneTSr, geneTSy, rows, lag, fit_method, hyper=hyper, has_reps=has_reps,
                                                only_array=only_array)

    # use global null
    elif null == "g":
        coefsr, interceptsr, fit_result_dfr = fit_all(geneTSr, geneTSy, rows, lag, fit_method, hyper=hyper, has_reps=has_reps,
                                                      only_array=only_array)
    else:
        raise ValueError("This shouldn't be reached")


    if save_prefix != None:
        pickle.dump(coefs, open(save_prefix + "_coefs.p", 'wB'))
        pickle.dump(intercepts, open(save_prefix + "_intercepts.p", 'wB'))
        pickle.dump(coefsr, open(save_prefix + "_coefsr.p", 'wB'))
        # pickle.dump(interceptsr, open(save_prefix + "_interceptsr.p", 'wB'))
        fit_result_df.to_csv(save_prefix + "_fit_result_df.txt", sep="\t", index=0)
        fit_result_dfr.to_csv(save_prefix + "_fit_result_dfr.txt", sep="\t", index=0)


        print "Coefs saved to ", save_prefix + "_coefs.p"
        print "Intercepts saved to ", save_prefix + "_intercepts.p"
        print "Fit result saved to ", save_prefix + "_fit_result_df.txt"
        print "Coefs-rand. saved to ", save_prefix + "_coefsr.p"
        # print "Intercepts-rand. saved to ", save_prefix + "_interceptsr.p"
        print "Fit result-rand. saved to ", save_prefix + "_fit_result_dfr.txt"




    return coefs, intercepts, fit_result_df, coefsr, fit_result_dfr


def run_bootstrap(geneTS, geneTSr, hyper, fit_method=fm.fit_lasso,
        lag=2,
        rows=None,
        save_prefix=None,
        has_reps=False,
        null="l", seed=None,
        only_array=True
    ):
    """
    Run one fit solely with the bootstrapped version. No null here. Since it's the same seed, each column will be subsampled in the same way.

    geneTS: n x T original gene matrix
    hyper: hyper-parameter
    fit_method: fit_method to use
    lag
    rows: geneTS will be the response variable here. This determines which rows of geneTS to select for response variables
        Defaulting None will just use all response variables as y
    save_prefix: how to save
    has_reps: if data has replicates
    null: "l" if to use the local, "g" if to use the global
    seed: Random seed

    Basically the same as run, except with bootstrapped columns.
    """

    if rows == None:
        rows = range(geneTS.shape[0])

    if null not in {"l", "g"}:
        raise ValueError("Null must either be l (local) or g (global)")

    geneTSy = geneTS[rows,:]

    print "Y shape: ", geneTSy.shape

    coefs, intercepts, fit_result_df = fit_all(geneTS, geneTSy, rows, lag, fit_method, hyper=hyper, has_reps=has_reps,
                                              bootstrap=True, seed=seed,
                                               only_array=only_array)


    # use local null
    if null == "l":
        coefsr, fit_result_dfr = fit_all_random(geneTS, geneTSr, geneTSy, rows, lag, fit_method, hyper=hyper, has_reps=has_reps,
                                               bootstrap=True, seed=seed,
                                                only_array=only_array)

    # use global null
    elif null == "g":
        coefsr, interceptsr, fit_result_dfr = fit_all(geneTSr, geneTSy, rows, lag, fit_method, hyper=hyper, has_reps=has_reps,
                                                     bootstrap=True, seed=seed,
                                                      only_array=only_array)
    else:
        raise ValueError("This shouldn't be reached. Null better be l or g")


    if save_prefix != None:
        pickle.dump(coefs, open(save_prefix + "_coefs.p", 'wB'))
        pickle.dump(intercepts, open(save_prefix + "_intercepts.p", 'wB'))
        pickle.dump(coefsr, open(save_prefix + "_coefsr.p", 'wB'))
        # pickle.dump(interceptsr, open(save_prefix + "_interceptsr.p", 'wB'))
        fit_result_df.to_csv(save_prefix + "_fit_result_df.txt", sep="\t", index=0)
        fit_result_dfr.to_csv(save_prefix + "_fit_result_dfr.txt", sep="\t", index=0)


        print "Coefs saved to ", save_prefix + "_coefs.p"
        print "Intercepts saved to ", save_prefix + "_intercepts.p"
        print "Fit result saved to ", save_prefix + "_fit_result_df.txt"
        print "Coefs-rand. saved to ", save_prefix + "_coefsr.p"
        # print "Intercepts-rand. saved to ", save_prefix + "_interceptsr.p"
        print "Fit result-rand. saved to ", save_prefix + "_fit_result_dfr.txt"




    return coefs, intercepts, fit_result_df, coefsr, fit_result_dfr



def get_bootstrap_results(bootstrap_matrs):
    """
    bootstrap_matrs: An iterable of coefficient matrices from bootstrapped data. All must have same dimension: n x n


    What about lags?

    Output:
    bootstrap_mean: cause x effect (n x n)
    bootstrap_std:
    bootstrap_freq: relative frequency of these
    """
    all_bootstrap_matr = np.stack(tuple(bootstrap_matrs), axis=0)

    bootstrap_mean = all_bootstrap_matr.mean(axis=0)
    bootstrap_std = all_bootstrap_matr.std(axis=0, ddof=1)
    bootstrap_freq = (all_bootstrap_matr != 0).mean(axis=0)





#     bootstrap_2std_ci = np.stack((bootstrap_mean - 2 * bootstrap_std,
#                                  bootstrap_mean + 2 * bootstrap_std),
#                                  axis=-1)

    # Turn it into numpy array

    return bootstrap_mean, bootstrap_std, bootstrap_freq


def bootstrap_matrices_iter_free(filenames):
    """
    note you don't have to align them first, can bootstrap
    :param filenames: filenames of numpy matrices

    recursive formula of getting mean, variance

    our parameter of interest, though is: 2.5%, 5%, 95% and 97.5%

    :return: first_matr
    """

    first_matr = pickle.load(open(filenames[0], 'rU'))

    mean_matr = first_matr.copy()
    count_matr = (first_matr !=0).astype(int)

    print "Bootstrap calculation by incremental adding + freeing"

    n = 2
    for i in range(1, len(filenames)):
        print i
        next_matr = pickle.load(open(filenames[i], 'rU'))

        if next_matr.shape == mean_matr.shape:

            if n == 2:
                var_matr = 1.0/n * (next_matr - mean_matr) ** 2
            else:
                var_matr = (n-2.0)/(n-1.0) * var_matr + 1.0 /n * ( next_matr - mean_matr)**2

            mean_matr = (1 - 1.0/n) * mean_matr + 1.0/n * next_matr

            count_matr += (next_matr != 0).copy()

            n += 1
        else:

            raise ValueError("Wrong shape: " +  str(next_matr.shape) + " From " + filenames[i])



        del next_matr

    std_matr = np.sqrt(var_matr)

    freq_matr = count_matr / (1.0 * (n-1))


    stats_matr_dict = collections.OrderedDict()
    stats_matr_dict["mean"] = mean_matr
    stats_matr_dict["std"] = std_matr
    stats_matr_dict["freq"] = freq_matr

    return stats_matr_dict

# Commented out below code since it will soon be obsolete. See the README_cpipeline for
# how to run the actual pipeline

# def get_parser():
#     # Parse arguments
#     import argparse
#
#     description = 'Run causal pipeline on input datasets'
#     parser = argparse.ArgumentParser(description=description)
#
#     parser.add_argument('-df', '--data_file', required=True)
#
#     parser.add_argument('-rdf', '--rand_data_file', required=True)
#
#     parser.add_argument('-o', '--out_prefix', required=True)
#
#     parser.add_argument('-hl', '--hyper_file', required=True)
#
#     parser.add_argument('-pa', '--plot_all', type=int, default=0)
#
#     #parser.add_argument('-nf', '--num_folds', type=int, required=True)
#
#     parser.add_argument('-t', '--test', required=True)
#
#     parser.add_argument('-l', '--lag',  type=int, required=True)
#
#     parser.add_argument('-sb', '--stratify_by', type=str, required=True)
#
#     return parser
#
# def load_and_run(args):
#
#     data_file = args.data_file
#     rand_data_file = args.rand_data_file
#     save_prefix = args.out_prefix
#
#
#     assert args.test in {'e', 'l', 'r'}
#     fit_method = test2fit_method[args.test]
#
#     lag = args.lag
#     #nfolds = args.num_folds
#     hyperlist = pickle.load(open(args.hyper_file, 'rB'))
#
#     assert args.stratify_by in {"e", "n"}
#
#     stratify_by = args2stratify_by[args.stratify_by]
#
#
#
#     # Load data file
#     df = pd.read_csv(data_file, sep="\t")
#     genes, geneTS = gtm.get_gene_TS(df)
#
#     dfr = pd.read_csv(rand_data_file, sep="\t")
#     genesr, geneTSr = gtm.get_gene_TS(dfr)
#
#     assert (geneTS.shape == geneTSr.shape)
#     assert (genes == genesr).all()
#
#
#
#     # run cross validate
#     best_hyper, best, hyper_df, hyper_fit_dfs = run_cross_validate(geneTS, fit_method=fit_method, hyperlist=hyperlist, lag=lag)
#
#
#     print "Best hyper is : ", best_hyper
#     print "Best result : ", best
#
#     print "Hyper df: "
#     print hyper_df
#
#
#     if not os.path.exists("hyper"):
#         os.makedirs("hyper")
#
#     pickle.dump(best_hyper, open("hyper" + os.sep + save_prefix + "_best_hyper.p", 'wB'))
#     best.to_csv("hyper" + os.sep + save_prefix + "_best_df.txt", sep="\t", index=0)
#     hyper_df.to_csv("hyper" + os.sep + save_prefix + "_hyper_df.txt", sep="\t", index=0)
#     for hyper, hyper_fit_df in zip(hyperlist, hyper_fit_dfs):
#         if isinstance(hyper, int) or isinstance(hyper, float):
#             suffix = "%.0E" % hyper
#         elif isinstance(hyper, tuple):
#             suffix = str(("%.0E" % h for h in hyper))
#
#         else:
#             raise ValueError("Hyper not float or tuple, how to plot?")
#
#         hyper_fit_df.to_csv("hyper" + os.sep + save_prefix + "_hyper_fit_df_" + suffix + ".txt", sep="\t", index=0)
#
#     print "Best hyper written to ", "hyper" + os.sep + save_prefix + "_best_hyper.p"
#     print "Best result written to ", "hyper" + os.sep + save_prefix + "_best_df.txt"
#     print "Hyper result written to ", "hyper" + os.sep + save_prefix + "_hyper_df.txt"
#     print "Hyper fit results written under ", "hyper" + os.sep + save_prefix + "_hyper_fit_df_"
#
#
#
#     if not os.path.exists("plots"):
#         os.makedirs("plots")
#     if not os.path.exists("plots" + os.sep + "hyper"):
#         os.makedirs("plots" + os.sep + "hyper")
#
#     # Plot the hyperparameters
#     plot_hyper_boxplot(hyperlist_to_labellist(hyperlist), hyper_fit_dfs, "r2", xlabel="Hyperparameters", ylabel="$r^2$", title="Hyperparameters VS $r^2$", filename="plots" + os.sep + "hyper"+ os.sep + "hyperVSr2")
#     plot_hyper_boxplot(hyperlist_to_labellist(hyperlist), hyper_fit_dfs, "mse", xlabel="Hyperparameters", ylabel="Mean-Squared Error", title="Hyperparameters VS MSE", filename="plots" + os.sep + "hyper" + os.sep + "hyperVSmse")
#     plot_hyper_boxplot(hyperlist_to_labellist(hyperlist), hyper_fit_dfs, "avg_df", xlabel="Hyperparameters", ylabel="Degrees of Freedom", title="Hyperparameters VS d.o.f.", filename="plots" + os.sep + "hyper" + os.sep + "hyperVSdof")
#
#
#     # Plot correlations
#     mse_vec = np.array([hyper_fit_df["mse"].values for hyper_fit_df in hyper_fit_dfs])
#     mse_corr = np.corrcoef(mse_vec)
#     gtm.save_gene_matrix("hyper" + os.sep + "mse_corr.txt", mse_corr, hyperlist)
#     print "MSE Correlation:"
#     print mse_corr
#     print "MSE corr. matrix saved to ", "hyper" + os.sep + "mse_corr.txt"
#
#     r2_vec = np.array([hyper_fit_df["r2"].values for hyper_fit_df in hyper_fit_dfs])
#     r2_corr = np.corrcoef(r2_vec)
#     gtm.save_gene_matrix("hyper" + os.sep + "r2_corr.txt", r2_corr, hyperlist)
#     print "R2 Correlation"
#     print r2_corr
#     print "R^2 corr. matrix saved to ", "hyper" + os.sep + "r2_corr.txt"
#
#
#     plot_corr_matrix(mse_corr, hyperlist_to_labellist(hyperlist), title="MSE correlation among hyperparams", filename="plots" + os.sep + "hyper" + os.sep + "mse_corr")
#     plot_corr_matrix(r2_corr, hyperlist_to_labellist(hyperlist), title="$r^2$ correlation among hyperparams", filename="plots" + os.sep + "hyper" + os.sep + "r2_corr")
#
#
#
#
#     # Run all
#
#     coefs, intercepts, fit_result_df, coefsr, interceptsr, fit_result_dfr = run(geneTS,
#                                                                            geneTSr,
#                                                                             hyper=best_hyper,
#                                                                             fit_method=fit_method,
#                                                                             lag=lag,
#                                                                             save_prefix=save_prefix
#                                                                            )
#
#
#
#     print "RESULTS"
#     print "*************************"
#     print "NORMAL: "
#     summarize_fit(coefs, intercepts, fit_result_df)
#
#
#     print "*************************"
#     print "RANDOM:"
#     summarize_fit(coefsr, interceptsr, fit_result_dfr)
#
#
#
#     # Align the coefs
#
#     acoefs = lc.align_coefs(coefs, lag)
#     acoefsr = lc.align_coefs(coefsr, lag)
#
#
#     print "Nonzero for:"
#     print "Original nonzero: ", len(np.nonzero(acoefs)[0])
#     print "Random nonzero: ", len(np.nonzero(acoefsr)[0])
#
#
#     # Save the gene matrices
#     for i in range(acoefs.shape[0]):
#         gtm.save_gene_matrix(filename=save_prefix + "_acoefs-lag-" + str(i+1) + ".txt", matrix=acoefs[i], genes=genes)
#         gtm.save_gene_matrix(filename=save_prefix + "_acoefsr-lag-" + str(i+1) + ".txt", matrix=acoefsr[i], genes=genes)
#
#
#         print "Lag ", i+1
#         print "Matrices written to:"
#         print save_prefix + "_acoefs-lag-" + str(i+1) + ".txt", ": Original"
#         print save_prefix + "_acoefsr-lag-" + str(i+1) + ".txt", ": Randomized"
#
#
#
#     if not os.path.exists("plots" + os.sep + "betas"):
#         os.makedirs("plots" + os.sep + "betas")
#
#     # Plot the betas
#     for i in range(acoefs.shape[0]):
#
#         fc.plot_betas(acoefs[i][np.nonzero(acoefs[i])].flatten(), acoefsr[i][np.nonzero(acoefsr[i])].flatten(), filename="plots" + os.sep + "betas" + os.sep + "beta_nonzero_lag-" + str(i+1), title="Causal coefs, Lag " + str(i+1))
#         fc.plot_betas(acoefs[i][np.nonzero(acoefs[i])].flatten(), acoefsr[i][np.nonzero(acoefsr[i])].flatten(), filename="plots" + os.sep + "betas" + os.sep + "beta_nonzero_lag-" + str(i+1) + "_zoom-in-95", zoom_in_percentile=95, title="Causal coefs, Lag " + str(i+1))
#
#         fc.plot_betas(np.absolute(acoefs[i][np.nonzero(acoefs[i])].flatten()), np.absolute(acoefsr[i][np.nonzero(acoefsr[i])].flatten()), filename="plots" + os.sep + "betas" + os.sep + "beta_abs_lag-" + str(i+1), title="Absolute causal coefs, Lag " + str(i+1))
#         fc.plot_betas(np.absolute(acoefs[i][np.nonzero(acoefs[i])].flatten()), np.absolute(acoefsr[i][np.nonzero(acoefsr[i])].flatten()), filename="plots" + os.sep + "betas" + os.sep + "beta_abs_lag-" + str(i+1) + "_zoom-in-95", zoom_in_percentile=95, title="Absolute causal coefs, Lag " + str(i+1))
#
#
#         print "Lag ", i+1
#         print "Plots written to: plots" + os.sep + "betas"
#
#
#
#
#
#     # get FDRS
#     fdrs = [0.01, 0.05, 0.1, 0.2]
#
#
#     acoefs_fdrs = []
#
#     for fdr in fdrs:
#
#         if not os.path.exists("fdr-" + str(fdr)):
#             os.makedirs("fdr-" + str(fdr))
#
#         fdr_prefix = "fdr-" + str(fdr) + os.sep + save_prefix
#
#         acoefs_fdr = np.zeros(acoefs.shape)
#
#         print "*************"
#         print "FDR = ", fdr
#         print "-----"
#         for i in range(acoefs.shape[0]):
#             print "Lag ", i + 1
#             acoefs_fdr[i], threshes = fc.get_abs_thresh(acoefs[i], acoefsr[i], fdr, stratify_by=stratify_by)
#             print "Threshes", threshes
#
#             matr_df = gtm.save_gene_matrix(filename=fdr_prefix + "_acoefs-lag-" + str(i+1) + "-fdr-" + str(fdr) + ".txt", matrix=acoefs_fdr[i], genes=genes)
#             pickle.dump(threshes, open(fdr_prefix + "_acoefs-lag-" + str(i+1) + "-fdr-" + str(fdr) + "-threshes.p", 'wB'))
#
#             net_df = nh.matr_to_net(matr_df, fdr_prefix + "-lag-" + str(i+1) + "-fdr-" + str(fdr) + "-sb-" + stratify_by, make_pair=False)
#             net_df.to_csv(fdr_prefix + "_acoefs-lag-" + str(i+1) + "-fdr-" + str(fdr) + "-network.txt", sep="\t", index=False)
#
#
#             print "Matrix written to:", fdr_prefix + "_acoefs-lag-" + str(i+1) + "-fdr-" + str(fdr) + ".txt", ": FDR at ", fdr
#             print "Threshes written to:", fdr_prefix + "_acoefs-lag-" + str(i+1) + "-fdr-" + str(fdr) + "-threshes.p"
#             print "Network written to:", fdr_prefix + "_acoefs-lag-" + str(i+1) + "-fdr-" + str(fdr) + "-network.txt"
#
#
#         print "FDR = ", fdr, " nonzero: ", len(np.nonzero(acoefs_fdr)[0])
#
#         acoefs_fdrs.append(acoefs_fdr.copy())
#
#
#
#     with open("matrices_done.txt", 'w') as donefile:
#         donefile.write("done\n")
#
#
#
#
#
#
#     for i, fdr in zip(range(len(fdrs)), fdrs):
#         acoefs_fdr = acoefs_fdrs[i]
#
#         if not os.path.exists("plots" + os.sep + "fdr-" + str(fdr)):
#             os.makedirs("plots" + os.sep + "fdr-" + str(fdr))
#
#         gtm.plot_coefs(acoefs_fdr, df, genes, lag, file_prefix="plots" + os.sep + "fdr-" + str(fdr) + os.sep + save_prefix+ "-",
#                   title_prefix="Lag " + str(lag) + ", fdr " + str(fdr) + " causal genes for ")
#
#         print "FDR plots written to: ", "plots" + os.sep + "fdr-" + str(fdr)
#
#
#
#     # Plot all the coefs
#     # NOTE: this will take a long time!
#     if args.plot_all:
#         if not os.path.exists("plots" + os.sep + "original"):
#             os.makedirs("plots" + os.sep + "original")
#         gtm.plot_coefs(acoefs, df, genes, lag, file_prefix="plots" + os.sep + "original" + os.sep + save_prefix + "-",
#                   title_prefix="Lag " + str(lag) + " causal genes for ")
#         print "Original plots written to: ", "plots" + os.sep + "original"
#
#         if not os.path.exists("plots" + os.sep + "randomized"):
#             os.makedirs("plots" + os.sep + "randomized")
#         gtm.plot_coefs(acoefsr, dfr, genes, lag, file_prefix="plots" + os.sep + "randomized" + os.sep + save_prefix+ "-",
#                   title_prefix="Lag " + str(lag) + " rand. causal genes for ")
#
#         print "Randomized plots written to: ", "plots" + os.sep + "randomized"
#
#
#
#
# def main():
#     load_and_run(get_parser().parse_args(sys.argv[1:]))
#
#
#
#
# if __name__ == '__main__':
#     main()





#
# def old_run_cross_validate(geneTS, fit_method=fm.fit_lasso,
#                         hyperlist=10**(-1 * np.arange(0, 4, 1.0)),
#                         lag=2,
#                         nfolds=5,
#                         save_prefix=None):
#     """
#     Writes out the results for a given hyper-parameter list
#     """
#
#
#     print "Hyper-parameters for cross-validation"
#     print hyperlist
#
#     # Cross-validate
#
#     best_hyper, best, hyper_df = cross_validate(geneTS, lag, fit_method, hyperlist, nfolds=nfolds)
#
#     print "Hypers list:"
#     print hyper_df
#
#     if save_prefix != None:
#         save_name = save_prefix + "_hyper_df.txt"
#         hyper_df.to_csv(save_name, sep="\t", index=0)
#         print "Hyper list written to ", save_name
#
#     return best_hyper, best, hyper_df


#
# def old_cross_validate(X_matr, lag, fit_method, hyperlist, nfolds=5, **kwargs):
#     """
#     X_matr: n x T matrix of genes
#     lag
#     test
#     hyperlist: list of settings of "hyper"
#
#     For each hyperparam setting, use that hyperparam in fitting the training (out-of-fold) on the test
#     Compute avg r^2. Among those above a r^2 threshold, take the hyperparam that results in sparsest.
#     Return results from each hyperparam setting and the best hyperparam
#
#     """
#
#     kf = KFold(X_matr.shape[0], n_folds=nfolds, shuffle=True, random_state=123)
#
#
#     total_sselist = []
#     avg_dflist = []
#     avg_r2list = []
#     avg_nlist = []
#
#     for hyper in hyperlist:
#
#         fit_results = []
#
#         for train, test in kf:
#             X_train = X_matr[train]
#             X_test = X_matr[test]
#
#             for i in range(X_test.shape[0]):
#                 Y_test = np.reshape(X_test[i,], (1, X_test.shape[1]))
#
#                 X_t, Y_t, Y_pred, coef, intercept, fit_result = fm.perform_test(X_matr=X_train, Y_matr=Y_test,
#                                                                              lag=lag, fit_method=fit_method,
#                                                                              hyper=hyper, **kwargs)
#
#                 fit_results.append(fit_result)
#
#         fit_result_df = pd.DataFrame(fit_results)
#
#         total_df = fit_result_df.sum()
#
#         avg_df = total_df / fit_result_df.shape[0]
#
#         total_sselist.append(total_df["sse"])
#         avg_dflist.append(avg_df["df"])
#         avg_r2list.append(avg_df["r2"])
#         avg_nlist.append(avg_df["n"])
#
#
#
#     hyper_df = pd.DataFrame({"hyper": hyperlist, "total_sse": total_sselist,
#                             "avg_df": avg_dflist, "avg_r2": avg_r2list,
#                             "lag": lag, "nfolds": nfolds,
#                             "avg_n": avg_nlist})
#
#     hyper_df.sort_values("total_sse", inplace=True)
#
#     best_hyper, best = get_best_hyper(hyper_df)
#
#     return best_hyper, best, hyper_df

# def get_best_hyper(hyper_df, min_r2=0.99):
#     hyper_df.sort_values("total_sse", inplace=True)
#
#     best_df = hyper_df[hyper_df["avg_r2"] > min_r2] if (hyper_df["avg_r2"] > min_r2).any() else hyper_df
#
#     best_df.sort_values("avg_df", inplace=True)
#
#     best = best_df.head(n=1)
#
#     best_hyper = best["hyper"].values[0]
#     return best_hyper, best