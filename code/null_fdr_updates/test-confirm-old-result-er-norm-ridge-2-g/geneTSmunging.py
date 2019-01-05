__author__ = 'jlu96'


try:
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    HAS_PLT = True
except ImportError, RuntimeError:
    HAS_PLT = False

import scipy.stats as stats
import pandas as pd
import csv
import collections
import numpy as np
import os

global timekeys
timekeys = ["t00", "t05", "t1_", "t2_", "t3_", "t4_", "t5_", "t6_", "t7_", "t8_", "t10", "t12"]

global repkeys
repkeys = ['t00_rep1', 't00_rep2plusextra', 't00_rep3plusextra', 't00_rep4', 't05_rep1', 't05_rep2', 't05_rep3', 't05_rep4', 't1_rep1', 't1_rep2', 't1_rep3', 't1_rep4', 't2_rep1', 't2_rep2', 't2_rep3', 't2_rep4', 't3_rep1', 't3_rep2', 't3_rep3', 't3_rep4', 't4_rep1', 't4_rep2', 't4_rep3', 't4_rep4', 't5_rep2', 't5_rep3', 't5_rep4', 't6_rep2', 't6_rep3', 't6_rep4', 't7_rep1', 't7_rep2', 't7_rep3', 't7_rep4', 't8_rep1', 't8_rep2', 't8_rep3', 't8_rep4', 't10_rep1', 't10_rep2', 't10_rep3', 't10_rep4', 't12_rep1', 't12_rep2', 't12_rep3', 't12_rep4']

rep2index = {}
for repkey in repkeys:
    rep2index[repkey] = timekeys.index(repkey[:3])



# the default scale of the plotting
global ylim_common_dist
ylim_common_dist = 7


def load_rep_file_list(rep_file):
    """
    :param rep_file: Lists which files in its directory to load
    :return: dfs, genes, geneTS

    Load a list of files at the same time and return them as replicates
    """

    dfs = []

    dir = os.path.dirname(rep_file)

    with open(rep_file, 'rU') as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")

        prev_df = None
        for row in reader:

            df = pd.read_csv(os.path.join(dir,row[0]), sep="\t")

            if prev_df != None:
                assert df["gene"] == prev_df["gene"]

            dfs.append(df)

    # subtract one to ignore the gene column
    geneTS = np.zeros((df.shape[0], df.shape[1] - 1, len(dfs)))

    for i, df in zip(range(len(dfs)), dfs):

        genes, geneTS[:, :, i] = get_gene_TS(df)


    timekeys, num_per_keys = get_shared_timekeys(dfs, index="gene")

    # join all the replicates by gene, then let index be normal
    df = dfs[0].set_index('gene').join([x.set_index('gene') for x in dfs[1:]])
    df.insert(0, "gene", df.index)
    df.reset_index(drop=True, inplace=True)
    df = avg_and_std(df, keys=timekeys)

    return dfs, genes, geneTS, df, timekeys, num_per_keys


def get_shared_timekeys(dfs, index=None):
    """
    Assume all the dfs are just independent but equally sized samples
    Get the time keys that are shared among the columns
    :param dfs:
    :param index: column to index by. None just uses normal index.
    :return: timekeys, num_per_keys
    """

    # [[col1A, col1B], [col2A, col2B], [col3A, col3B]]
    if index != None:
        collists = [list(x) for x in np.array([df.drop(index, 1).columns.values for df in dfs]).T]
    else:
        collists = [list(x) for x in np.array([df.columns.values for df in dfs]).T]

    timekeys = [os.path.commonprefix(collist) for collist in collists]

    join_df = dfs[0].set_index('gene').join([x.set_index('gene') for x in dfs[1:]])

    num_per_keys = [len([x for x in join_df.columns.values if x.startswith(t)]) for t in timekeys]

    return timekeys, num_per_keys



def avg_and_std(df, keys=timekeys):
    """
    Get the avg and std of data over time (axis 1)
    :param df: loaded file of where each col is a timepoint replciate, eg. t00_rep1
    :param keys: List of column prefixes to avg/std over
    :return: Updated version
    """

    for key in keys:
        cols = [col for col in list(df.columns.values) if col[:len(key)] == key]
        df[key] = pd.Series(sum([df[col] for col in cols]) * 1.0 / len(cols), index=df.index)
        std = np.std([df[col] for col in cols], axis=0)
        df[key + 'std'] = std
    return df








# Try not to use the below anymore
def load_file_and_avg(filename, keys = timekeys, all_avg=True, fold_change=False,
                      diff=False, normal_diff=False, verbose=False):
    raise ValueError("Rewrite thie function before calling.")
    data = pd.read_csv(filename, sep="\t")
    num_per_key = []

    data = avg_and_std(data, keys)

    if all_avg:
        data["avg"] = np.average(data[keys], axis=1)

    # normalized
    # look at differences

    if fold_change:
        print "Fold change loaded"
        fold_keys = []
        for i in range(len(keys) - 1):
            key1 = keys[i]
            key2 = keys[i + 1]
            fold_key = key1 + "-" + key2 + " fold"
            data[fold_key] = data[key2] * 1.0/ data[key1]
            fold_keys.append(fold_key)
    if diff:
        print "Diff loaded"
        diff_keys = []
        for i in range(len(keys) - 1):
            key1 = keys[i]
            key2 = keys[i + 1]
            diff_key = key1 + "-" + key2 + " diff"
            data[diff_key] = data[key2] - data[key1]
            diff_keys.append(diff_key)



    if normal_diff:
        print "normalized diff loaded"
        normal_diff_keys = []
        data["Mean_diff"] = data[diff_keys].mean(axis=1)
        data["Std_diff"] = data[diff_keys].std(axis=1)

        for i, diff_key in zip(range(len(keys) - 1), diff_keys):
            key1 = keys[i]
            key2 = keys[i + 1]
            normal_diff_key = key1 + "-" + key2 + " normal_diff"
            data[normal_diff_key] = (data[diff_key] - data["Mean_diff"])/(data["Std_diff"])
            normal_diff_keys.append(normal_diff_key)



        
    if verbose:
        for key, key_num in zip(keys, num_per_key):
            print key, "has", key_num, "data points"

    return data

def get_gene_TS(data, genes=None, gene_col="gene"):
    """
    :param data: pandas dataframe with all keys as timepoint except the gene col
    :param genes: genes to extract
    :param gene_col:
    :return: geneTS, the matrix of timepoints
    """
    if genes == None:
        genes = data["gene"]
    cols = list(data.columns.values)
    cols.remove(gene_col)

    index = data[gene_col].isin(genes)

    found_data = data[index]

    found_genes = found_data[gene_col].values

    return found_genes, found_data[cols].as_matrix()


# Old get_gene_TS code
# def get_gene_TS(data, genes, keys=timekeys):
#     """Return genes found
#     and geneTS an n x t matrix of expression levels where n is number of genes and t is number of timepoints, i.e. keys
#     """
#
#     gene_set = set(genes)
#
#     key_indices = np.where([value in keys for value in data.columns.values])[0]
#
#
#     geneTS = data[data.gene.isin(gene_set)][key_indices].values
#     found_genes = data[data.gene.isin(gene_set)]['gene'].values
#
#     return found_genes, geneTS

if HAS_PLT:
    def plot_genes(data, genes, gene_col="gene", keys = timekeys, title="Gene Expression From a few species",
                   num_per_keys  = None, line_color_labels=None, ylabel="Expression Level",
                   horizontal_line_color_labels=None, xlim=None, ylim=None, filename=None, plot_bar=True, plot_reps=False, repkeys=repkeys,
                   gene_labels=None,
                   plot_outside=False,
                   verbose=False,
                   legend_fontsize=None,
                   linewidth=3,
                   capsize=8,
                   capwidth=2,
                   bar_linestyle="--",
                   show_plot=True,
                   xlabel="Time (hr)",
                    legend=True,
                   no_axis=False,
                   figsize=(12,8)):
        """
        Plots genes optionally with replicates

        Error bar code pulled from http://stackoverflow.com/questions/20033396/how-to-visualize-95-confidence-interval-in-matplotlib


        :param data: dataframe, with one column as the gene
        :param genes: list of genes to plot
        :param keys: timepoint keys, in order, to plot. Their entries should have the timepoint avg.
        They themselves should be prefixes of the replicates, if there are any.
        :param title:
        :param num_per_keys:
        :param line_color_labels:
        :param ylabel:
        :param horizontal_line_color_labels:
        :param xlim:
        :param ylim:
        :param filename:
        :param plot_reps:
        :param repkeys:
        :param gene_labels:
        :param plot_outside:
        :param verbose:
        :param legend_fontsize:
        :return:
        """

        if plot_bar:

            if num_per_keys == None:
                raise ValueError("Need to include the num per key for plot_bar")
            else:
                if len(num_per_keys) != len(timekeys):
                    raise ValueError("Need the # reps for each timekey")

            # print "Rep # per timepoint unset. Inferring by counting 'rep'"
            #
            # num_per_keys = np.array([len([x for x in data.columns.values if "rep" in x and x.startswith(timekey)])
            #                 for timekey in timekeys])
            # print list(num_per_keys)


        fig = plt.figure(figsize=figsize)
        ax = plt.subplot(111)
        for gene, i in zip(genes, range(len(genes))):
            # time points, average for the index, std for index
            index = data[gene_col] == gene

            if gene_labels != None:
                label = gene_labels[i]
            else:
                label = gene

            gene_avg = data[index][keys].values.flatten()

            if verbose:
                print "Gene: ", gene
                print gene_avg

            # try to plot the bar here
            if plot_bar:
                try:
                    gene_std = data[index][[key + 'std' for key in keys]].values.flatten()


                    (_, caps, eb) = plt.errorbar(range(len(keys)), gene_avg,
                             yerr=stats.t.ppf(0.95, np.array(num_per_keys) - 1) * gene_std, label=label,
                                 linewidth=linewidth, capsize=capsize)
                    for cap in caps:
                        cap.set_markeredgewidth(capwidth)

                    eb[0].set_linestyle(bar_linestyle)

                except KeyError:
                    plt.errorbar(range(len(keys)), gene_avg, label=label,
                                 linewidth=linewidth)
            else:
                plt.errorbar(range(len(keys)), gene_avg, label=label,
                                 linewidth=linewidth)



            if plot_reps:
                plt.scatter([rep2index[r] for r in repkeys], data[index][repkeys])

        # Shrink current axis by 20%




        if line_color_labels != None:
            for line, color, label in line_color_labels:
                plt.axvline(line, color=color,label=label)
        if horizontal_line_color_labels !=None:
            for line, color, label in horizontal_line_color_labels:
                plt.axhline(line, color=color,label=label, linestyle='dashed')

        if xlim != None:
            plt.xlim(xlim[0], xlim[1])
        else:
            plt.xlim(-1, 12)
        if ylim != None:
            plt.ylim(ylim[0], ylim[1])




        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.xlabel(xlabel, fontsize=20)
        plt.ylabel(ylabel, fontsize=20)
        plt.title(title, fontsize=20)

        #plt.tight_layout(w_pad=5)

        # x0, x1, y0, y1 = plt.axis()
        # plt.axis((x0 - x_margin,
        #           x1 + x_margin,
        #           y0 - y_margin,
        #           y1 + y_margin))

        if no_axis:
            plt.axis('off')

        if legend:
            if plot_outside:
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
                        # Put a legend to the right of the current axis
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=legend_fontsize)
            else:
                plt.legend(loc='best', fontsize=legend_fontsize)

        if filename:
            fig.savefig(filename)
            print "Plot saved to ", filename

        if show_plot:
            plt.show()
        plt.close()





    def plot_gene_pairs(data, gene_pairs, keys = timekeys,
                        title="Paired Gene Expression From a few species",
                   num_per_keys  = np.array([4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4]), no_error = False,
                        ylabel="Gene Expression", filename=None):
        """Assumes dataframe has gene and time point averages` in first row."""
        colors = ['red', 'blue', 'green', 'cyan', 'magenta', 'yellow']
        fig = plt.figure(figsize=(12,8))
        ax = plt.subplot(111)
        for gene_pair, i in zip( gene_pairs, range(len(gene_pairs))):
            # time points, average for the index, std for index

            gene1, gene2 = gene_pair

            print gene1, gene2
            gene1_avg = data[data['gene'] == gene1][keys].values.flatten()
            gene2_avg = data[data['gene'] == gene2][keys].values.flatten()
            if no_error:
                plt.errorbar(range(len(keys)), gene1_avg, label=gene1, color=colors[i])
            else:
                try:
                    gene1_std = data[data['gene'] == gene1][[key + 'std' for key in keys]].values.flatten()

                    plt.errorbar(range(len(keys)), gene1_avg,
                             yerr=stats.t.ppf(0.95, num_per_keys - 1) * gene1_std, label=gene1, color=colors[i])
                except KeyError:
                    plt.errorbar(range(len(keys)), gene1_avg, label=gene1, color=colors[i])

            if no_error:
                plt.errorbar(range(len(keys)), gene2_avg,  color=colors[i], linestyle='dashed')
            else:
                try:
                    gene2_std = data[data['gene'] == gene2][[key + 'std' for key in keys]].values.flatten()

                    plt.errorbar(range(len(keys)), gene2_avg,
                             yerr=stats.t.ppf(0.95, num_per_keys - 1) * gene2_std,  color=colors[i], linestyle='dashed')
                except KeyError:
                    plt.errorbar(range(len(keys)), gene2_avg,  color=colors[i], linestyle='dashed')


        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])





        # Put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        plt.xlabel("Time points", fontsize=20)
        plt.ylabel(ylabel, fontsize=20)
        plt.title(title, fontsize=20)
        plt.xlim(-1, 12)
        if filename:
            fig.savefig(filename)
            print "Plot saved to ", filename

        plt.show()
        plt.close()

    def plot_multi_pairs(data, pairs, title_prefix, ylabel="Log-TPM Expression", min_plot=None, max_plot=None, num_plot=3, step_size=None, file_prefix=None):
        """Plot several pair-plots
        """
        if not min_plot:
            min_plot = 0
        if not max_plot:
            max_plot = len(pairs)
        if not step_size:
            step_size = num_plot


        plot_indices = range(min_plot, max_plot, step_size)


        for i, num in zip(plot_indices, range(1, len(plot_indices)+1)):
            title = title_prefix + "," + str(num)
            if file_prefix:
                filename = file_prefix + "-" + str(num)

            plot_pairs = pairs[i: min(i+num_plot, max_plot)]
            plot_pairs = [tuple(x) for x in plot_pairs]
            if file_prefix:
                plot_gene_pairs(data, plot_pairs, ylabel=ylabel, title=title, filename=filename)
            else:
                plot_gene_pairs(data, plot_pairs, ylabel=ylabel, title=title)



    def plot_genes_grid(data, genes, keys = timekeys, title="Gene Expression From a few species",
               num_per_keys  = np.array([4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4]), line_color_labels=None, ylabel="Expression Level",
               horizontal_line_color_labels=None, xlim=None, ylim=None, file_prefix=None, xlabel="Timepoints",
                   colors = ['r', 'g', 'b'], styles = ["solid", "solid", "solid"], figsize=(8,6),
                        sort_names=True, plot_reps=True):
        """Assumes dataframe has gene and time point averages` in first row."""
        if len(genes) > 4:
            raise ValueError("# genes must <= 4")

        # assign colors

        geneToColor = dict(zip(genes, colors))
        print "Gene to color is: ", geneToColor


        geneToStyle = dict(zip(genes, styles))
        print "Gene to style is ", geneToStyle

        # get the min and max of the data
        ymin_data = np.min(data[data['gene'].isin(genes)][keys].values)
        ymax_data = np.max(data[data['gene'].isin(genes)][keys].values)
        yavg_data = np.average(data[data['gene'].isin(genes)][keys].values)
        ystd_data = np.std(data[data['gene'].isin(genes)][keys].values)
        print "Y min is : ", ymin_data
        print "Y max is : ", ymax_data
        print "Y avg is : ", yavg_data
        print "Y std is : ", ystd_data



        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex="col", figsize=figsize)


        # Main plot
        for gene in genes:
            # time points, average for the index, std for index
            gene_avg = data[data['gene'] == gene][keys].values.flatten()
            try:
                gene_std = data[data['gene'] == gene][[key + 'std' for key in keys]].values.flatten()

                ax1.errorbar(range(len(keys)), gene_avg,
                         yerr=stats.t.ppf(0.95, num_per_keys - 1) * gene_std, color=geneToColor[gene], linestyle=geneToStyle[gene])
            except KeyError:
                ax1.errorbar(range(len(keys)), gene_avg, label=gene, color=geneToColor[gene], linestyle=geneToStyle[gene])



        if line_color_labels != None:
            for line, color, label in line_color_labels:
                ax1.axvline(line, color=color,label=label)
        if horizontal_line_color_labels !=None:
            for line, color, label in horizontal_line_color_labels:
                ax1.axhline(line, color=color,label=label, linestyle="dashed")

        if xlim != None:
            ax1.set_xlim(xlim[0], xlim[1])
        else:
            ax1.set_xlim(-1, 12)
        if ylim != None:
            ax1.set_ylim(ylim[0], ylim[1])
        else:
            ycenter = (ymin_data + ymax_data)/2.0
            ax1.set_ylim(ycenter - ylim_common_dist * 0.5, ycenter + ylim_common_dist * 0.5)

        ax1.legend(loc="best", fontsize=10)
        ax1.set_title("All genes")

        # Subplots
        for gene, ax in zip(genes, [ax2, ax3, ax4]):
            gene_avg = data[data['gene'] == gene][keys].values.flatten()
            try:
                gene_std = data[data['gene'] == gene][[key + 'std' for key in keys]].values.flatten()

                ax.errorbar(range(len(keys)), gene_avg,
                         yerr=stats.t.ppf(0.95, num_per_keys - 1) * gene_std, color=geneToColor[gene], linestyle=geneToStyle[gene])
            except KeyError:
                ax.errorbar(range(len(keys)), gene_avg, label=gene, color=geneToColor[gene], linestyle=geneToStyle[gene])
            if plot_reps:
                index = data['gene'] == gene
                ax.scatter([rep2index[r] for r in repkeys], data[index][repkeys], color=geneToColor[gene], s=5)




            ax.set_title(gene)
            if xlim != None:
                ax.set_xlim(xlim[0], xlim[1])
            else:
                ax.set_xlim(-1, 12)


        fig.text(0.54, 0.03, xlabel, ha='center', fontsize=20)
        fig.text(0.03, 0.5, ylabel, va='center', rotation='vertical', fontsize=20)

        # Put a legend to the right of the current axis
        fig.suptitle(title,  fontsize=20, y=0.98, x=0.54)

        plt.tight_layout()
        fig.subplots_adjust(bottom=0.11, left=0.11, top=0.89)

        if file_prefix:
            if sort_names:
                gene_string = "-".join(sorted(genes))
            else:
                gene_string = "-".join(genes)

            filename = file_prefix + "-" + str(int(yavg_data)) + "-" + str(int(ystd_data)) + "-" + gene_string


            print "Figure saved to ", filename
            fig.savefig(filename)

        plt.show()
        plt.close()


    def plot_multi_genes_grid(data, genes, title, ylabel="Log-TPM Expression", min_plot=None, max_plot=None, num_plot=3, step_size=None, file_prefix=None,
                              horizontal_line_color_labels=None):
        """Plot several pair-plots
        """
        if not min_plot:
            min_plot = 0
        if not max_plot:
            max_plot = len(genes)
        if not step_size:
            step_size = num_plot


        plot_indices = range(min_plot, max_plot, step_size)


        for i, num in zip(plot_indices, range(1, len(plot_indices)+1)):

            plot_genes = genes[i: min(i+num_plot, max_plot)]


            if file_prefix:
                plot_genes_grid(data, plot_genes, ylabel=ylabel, title=title, horizontal_line_color_labels=horizontal_line_color_labels, file_prefix=file_prefix)
            else:
                plot_genes_grid(data, plot_genes, ylabel=ylabel, title=title, horizontal_line_color_labels=horizontal_line_color_labels)



    def plot_pairs_grid(data, pairs, keys = timekeys, title="Gene Expression From a few species",
               num_per_keys  = np.array([4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4]), line_color_labels=None, ylabel="Expression Level",
               horizontal_line_color_labels=None, xlim=None, ylim=None, file_prefix=None, xlabel="Timepoints",
                   colors = ['r', 'g', 'b'], figsize=(8,6), no_error=False):
        """Assumes dataframe has gene and time point averages` in first row."""
        if len(pairs) > 4:
            raise ValueError("# pairs must <= 4")

        # assign colors

        genes = set()
        for pair in pairs:
            genes.add(pair[0])
            genes.add(pair[1])


        pairToColor = dict(zip(pairs, colors))
        print " pair to color is: ",  pairToColor

        print pairs
        print genes

        # get the min and max of the data
        ymin_data = np.min(data[data['gene'].isin(genes)][keys].values)
        ymax_data = np.max(data[data['gene'].isin(genes)][keys].values)
        yavg_data = np.average(data[data['gene'].isin(genes)][keys].values)
        ystd_data = np.std(data[data['gene'].isin(genes)][keys].values)
        print "Y min is : ", ymin_data
        print "Y max is : ", ymax_data
        print "Y avg is : ", yavg_data
        print "Y std is : ", ystd_data



        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex="col", figsize=figsize)



        for pair, i in zip(pairs, range(len(pairs))):
            # time points, average for the index, std for index

            gene1, gene2 = pair

            print gene1, gene2
            gene1_avg = data[data['gene'] == gene1][keys].values.flatten()
            gene2_avg = data[data['gene'] == gene2][keys].values.flatten()
            if no_error:
                ax1.errorbar(range(len(keys)), gene1_avg,  color=pairToColor[pair])
            else:
                try:
                    gene1_std = data[data['gene'] == gene1][[key + 'std' for key in keys]].values.flatten()

                    ax1.errorbar(range(len(keys)), gene1_avg,
                             yerr=stats.t.ppf(0.95, num_per_keys - 1) * gene1_std,  color=pairToColor[pair])
                except KeyError:
                    ax1.errorbar(range(len(keys)), gene1_avg,  color=pairToColor[pair])

            if no_error:
                ax1.errorbar(range(len(keys)), gene2_avg,  color=pairToColor[pair], linestyle='dashed')
            else:
                try:
                    gene2_std = data[data['gene'] == gene2][[key + 'std' for key in keys]].values.flatten()

                    ax1.errorbar(range(len(keys)), gene2_avg,
                             yerr=stats.t.ppf(0.95, num_per_keys - 1) * gene2_std,  color=pairToColor[pair], linestyle='dashed')
                except KeyError:
                    ax1.errorbar(range(len(keys)), gene2_avg,  color=pairToColor[pair], linestyle='dashed')
        if line_color_labels != None:
            for line, color, label in line_color_labels:
                ax1.axvline(line, color=color,label=label)
        if horizontal_line_color_labels !=None:
            for line, color, label in horizontal_line_color_labels:
                ax1.axhline(line, color=color,label=label, linestyle="dashed")

        if xlim != None:
            ax1.set_xlim(xlim[0], xlim[1])
        else:
            ax1.set_xlim(-1, 12)
        if ylim != None:
            ax1.set_ylim(ylim[0], ylim[1])
        else:
            ycenter = (ymin_data + ymax_data)/2.0
            ax1.set_ylim(ycenter - ylim_common_dist * 0.5, ycenter + ylim_common_dist * 0.5)

        ax1.legend(loc="best", fontsize=10)
        ax1.set_title("All genes")

        # Subplots
        for pair, ax in zip(pairs, [ax2, ax3, ax4]):
            gene1, gene2 = pair

            print gene1, gene2
            gene1_avg = data[data['gene'] == gene1][keys].values.flatten()
            gene2_avg = data[data['gene'] == gene2][keys].values.flatten()
            if no_error:
                ax.errorbar(range(len(keys)), gene1_avg,  label=gene1, color=pairToColor[pair])
            else:
                try:
                    gene1_std = data[data['gene'] == gene1][[key + 'std' for key in keys]].values.flatten()

                    ax.errorbar(range(len(keys)), gene1_avg,
                             yerr=stats.t.ppf(0.95, num_per_keys - 1) * gene1_std,  label=gene1,color=pairToColor[pair])
                except KeyError:
                    ax.errorbar(range(len(keys)), gene1_avg, label=gene1, color=pairToColor[pair])

            if no_error:
                ax.errorbar(range(len(keys)), gene2_avg, label=gene2, color=pairToColor[pair], linestyle='dashed')
            else:
                try:
                    gene2_std = data[data['gene'] == gene2][[key + 'std' for key in keys]].values.flatten()

                    ax.errorbar(range(len(keys)), gene2_avg,
                             yerr=stats.t.ppf(0.95, num_per_keys - 1) * gene2_std, label=gene2, color=pairToColor[pair], linestyle='dashed')
                except KeyError:
                    ax.errorbar(range(len(keys)), gene2_avg, label=gene2, color=pairToColor[pair], linestyle='dashed')
            if xlim != None:
                ax.set_xlim(xlim[0], xlim[1])
            else:
                ax.set_xlim(-1, 12)

            ax.set_title(gene1 + ", " + gene2)
            ax.legend(loc="best", fontsize=10)


        fig.text(0.54, 0.03, xlabel, ha='center', fontsize=20)
        fig.text(0.03, 0.5, ylabel, va='center', rotation='vertical', fontsize=20)

        # Put a legend to the right of the current axis
        fig.suptitle(title,  fontsize=20, y=0.98, x=0.54)

        plt.tight_layout()
        fig.subplots_adjust(bottom=0.11, left=0.11, top=0.89)

        if file_prefix:

            gene_string = "-".join(sorted(genes))

            filename = file_prefix + "-" + str(int(yavg_data)) + "-" + str(int(ystd_data)) + "-" + gene_string


            print "Figure saved to ", filename
            fig.savefig(filename)

        plt.show()
        plt.close()

    def plot_multi_pairs_grid(data, pairs, title, ylabel="Log-TPM Expression", min_plot=None, max_plot=None, num_plot=3, step_size=None, file_prefix=None,
                                  horizontal_line_color_labels=None):
            """Plot several pair-plots
            """
            if not min_plot:
                min_plot = 0
            if not max_plot:
                max_plot = len(pairs)
            if not step_size:
                step_size = num_plot


            plot_indices = range(min_plot, max_plot, step_size)


            for i, num in zip(plot_indices, range(1, len(plot_indices)+1)):

                plot_pairs = pairs[i: min(i+num_plot, max_plot)]


                if file_prefix:
                    plot_pairs_grid(data, plot_pairs, ylabel=ylabel, title=title, horizontal_line_color_labels=horizontal_line_color_labels, file_prefix=file_prefix)
                else:
                    plot_pairs_grid(data, plot_pairs, ylabel=ylabel, title=title, horizontal_line_color_labels=horizontal_line_color_labels)


    def plot_timepoint_histogram(matr, x_label="Value", title_prefix="Histogram of timepoint ", bins=30,
                                time_names=timekeys,
                                same_axes=True, percentile_zoom = None, xlim=None, save_prefix=None, line_color_labels=None,
                                 horizontal_line_color_labels=None, ylim=None):
        """
        Given matr, an N x T matrix, plot the histogram of values at each timepoint t.
        percentile_zoom: removes the tails
        """

        matr = np.array(matr)
        T = matr.shape[1]

        if time_names == None:
            time_names = [str(t) for t in range(T)]

        trans_matr = matr.T


        if percentile_zoom != None:
            trans_matr = []
            for t in range(T):
                top_percentile = stats.scoreatpercentile(matr[:, t], 100 - percentile_zoom)
                bottom_percentile = stats.scoreatpercentile(matr[:, t], percentile_zoom)

                new_array = matr[:,t]
                new_indices = np.where((new_array < top_percentile) & (new_array > bottom_percentile))[0]
                trans_matr.append(new_array[new_indices])

            trans_matr = np.array(trans_matr)

        if same_axes:
            if not xlim:
                min_value = np.min(np.concatenate(tuple(trans_matr)))
                max_value = np.max(np.concatenate(tuple(trans_matr)))
            else:
                min_value = xlim[0]
                max_value = xlim[1]
            bins = np.linspace(min_value, max_value, bins)
            print bins

        figs = []
        for t, time_name in zip(range(T), time_names):
            title = title_prefix + time_name
            fig = plt.figure(figsize=(8,8))
            plt.hist(trans_matr[t], bins=bins)

            plt.title(title, fontsize=20)
            plt.xlabel(x_label, fontsize=20)
            plt.ylabel("Frequency", fontsize=20)
            if line_color_labels != None:
                for line, color, label in line_color_labels:
                    plt.axvline(line, color=color,label=label)
                plt.legend(loc='best')

            if horizontal_line_color_labels !=None:
                for line, color, label in horizontal_line_color_labels:
                    plt.axhline(line, color=color,label=label)
                plt.legend(loc='best')


            if save_prefix:
                filename = save_prefix + time_name
                fig.savefig(filename)
                print "Plot saved to ", filename
            figs.append(fig)
            plt.show()
            plt.close()

        if save_prefix:
            print "All images saved with prefix:", save_prefix

        return figs




    # MOVED TO CPIPELINE
    # def plot_all_coefs(acoefs, df, genes, lag, min_coef=0.01, num_per_plot=3, file_prefix=None, title_prefix="Causal genes for ", verbose=False,
    #                    **kwargs):
    #     """
    #     acoefs: aligned coefs, of form (lag x n x n), acoefs[i] is lag_i+1
    #     df: Original TS
    #     genes: list of all genes
    #     """
    #
    #     # iterate down the coefs of each effect gene. This is one effect gene per column.
    #     for j in range(acoefs.shape[2]):
    #         out_gene = genes[j]
    #
    #         acoef = acoefs[:, :, j]
    #         preds = np.where(np.absolute(acoef) > min_coef)
    #
    #         lags = preds[0] + 1
    #         cgenes = genes[preds[1]]
    #         cos = np.round(acoef, 2)[preds]
    #
    #         num_plots = int(np.ceil(len(cgenes) * 1.0 / num_per_plot))
    #
    #         if verbose:
    #             print "Out gene: ", out_gene
    #
    #         for i in range(num_plots):
    #             plags = lags[i * len(cgenes)/num_plots: (i+ 1) * len(cgenes)/num_plots]
    #             pgenes = cgenes[i * len(cgenes)/num_plots: (i+ 1) * len(cgenes)/num_plots]
    #             pcos = cos[i * len(cgenes)/num_plots: (i+ 1) * len(cgenes)/num_plots]
    #             labels = ["{0:>20}".format(out_gene + ":") + " Coef, Lag" ] + ["{0:>20}".format(pgene + ":") + " " + str(pco) + ", " + str(plag) for pgene, pco, plag in zip(pgenes, pcos, plags)]
    #
    #             if verbose:
    #                 print "Part: ", i + 1
    #                 print "Lag points: ", plags
    #                 print "Pred genes:", pgenes
    #                 print "Pred coefs: ", pcos
    #                 print "Labels are ", labels
    #
    #             plot_genes(df, [out_gene] + list(pgenes), title=title_prefix + out_gene + " , Part " + str(i+1),
    #                           filename=None if file_prefix == None else file_prefix + out_gene.replace(".", ",") + "_lag-" + str(lag) + \
    #                                                                     "-" + "-".join([x.replace(".", ",") for x in list(pgenes)]),
    #                           gene_labels=labels, plot_outside=False,
    #                        **kwargs)




def get_sig_gene_pairs(sig_matr, genes):
    """
    :param sig_matr: Matrix of indices where significant
    :param genes: List of genes
    :return: List of pairs of gene by indices
    """

    rows, cols = np.where(sig_matr)

    grows, gcols = [genes[row] for row in rows], [genes[col] for col in cols]

    return zip(grows, gcols)



def compare_sig_matr(sig_matr_list):
    """
    Returns:
    true-values across all matrices
    tuples of # in, same, sorted
    total # sig
    """

    num = len(sig_matr_list)

    sig_matr_sum = np.sum(np.array(sig_matr_list), axis=0)

    all_sig_matr = sig_matr_sum >= num

    all_sig_num = len(np.where(all_sig_matr.flatten())[0])

    not_sig_num = len(np.where(sig_matr_sum.flatten())[0]) - all_sig_num

    return all_sig_matr, all_sig_num, not_sig_num

def randomize_geneTS(geneTS, seed=123, replace=False, num_genes=None, by_time=True):
    """
    Given an n x T (T is timepoints) matrix of gene expression levels over time, return a new matrix with randomization.

    Replace=True means don't sample with replacement the new values
    """

    np.random.seed(seed)

    if num_genes==None:
        num_genes = geneTS.shape[0]
        T = geneTS.shape[1]

    if by_time:

        newgeneTS = np.array([np.random.choice(geneTS[i], size=T, replace=replace) for i in range(num_genes)])

        return newgeneTS

    else:
        TbyG = geneTS.T


        newTbyG = np.array([np.random.choice(TbyG[i], size=num_genes, replace=replace) for i in range(TbyG.shape[0])])

        newgeneTS = newTbyG.T

        return newgeneTS

def make_and_save_randomized_data(data, seed=123, filename=None):
    """
    :param data: A pandas dataframe with a column named "gene" and remaining columns are timepoints
    :param filename:
    :return:
    """
    genes, geneTS = get_gene_TS(data, data["gene"])
    rand_geneTS = randomize_geneTS(geneTS, seed=seed)

    keys = list(data.columns.values)
    keys.remove("gene")
    keys = np.array(keys)

    rand_data_dict = {}
    rand_data_dict['gene'] = data['gene'].values
    for i, key in zip(range(len(keys)), keys):
        rand_data_dict[key] = rand_geneTS[:, i]
    rand_data = pd.DataFrame(data=rand_data_dict, columns=["gene"] + list(keys))

    if filename == None:
        print "Randomized not written"
    else:
        rand_data.to_csv(filename, sep='\t', index=False)
        print "Randomized written to ", filename

    return rand_data


def make_gene_matrix(matrix, genes):
    data_dict = collections.OrderedDict()
    for i, gene in zip(range(len(genes)), genes):
        data_dict[gene] = matrix[:, i]
    df = pd.DataFrame(data_dict)
    df["gene"] = genes

    df = df.set_index("gene")

    return df

def save_gene_matrix(filename, matrix, genes):

    df = make_gene_matrix(matrix, genes)

    df.to_csv(filename,sep="\t")

    return df


