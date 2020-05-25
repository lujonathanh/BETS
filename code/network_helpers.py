__author__ = 'jlu96'
import csv
import numpy as np
import json
import pandas as pd
import collections

global data_dir
data_dir = "/Users/jlu96/v-causal-snps/data/GeneExpressionData/GGR_Network/"

global genecols
genecols = ["Official Symbol Interactor A", "Official Symbol Interactor B"]


def load_hg_ensg_old():
    hg_ensg_file = data_dir +  "HGNC-to-Ensembl.txt"

    hg2ensg = {}
    with open(hg_ensg_file, 'rU') as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")

        for row in reader:
            if row["Ensembl Gene ID"]:
                hg2ensg[row["Approved Symbol"]] = row["Ensembl Gene ID"]
    ensg2hg = dict([(item[1], item[0]) for item in list(hg2ensg.items()) if item[0] != ""])

    return hg2ensg, ensg2hg


def load_hg_ensg(ensg_hg_file = data_dir +  "gencode.v22.gene_id_to_gene_name.json"
):
    with open(ensg_hg_file , 'rU') as jsonfile:
        ensg2hg = json.load(jsonfile)

    hg2ensg = dict([(item[1], item[0]) for item in list(ensg2hg.items()) if item[0] != ""])

    return hg2ensg, ensg2hg




def hg2ensg(gene, verbose=False):
    try:
        conv = hg2ensg.conv
    except AttributeError:
        hg2ensg.conv, _ = load_hg_ensg()
        conv = hg2ensg.conv

    try:
        return conv[gene]
    except KeyError:
        if verbose:
            print("Gene ", gene, " missing from hg2ensg")
        return ""

def ensg2hg(gene, verbose=False):
    try:
        conv = ensg2hg.conv
    except AttributeError:
        _, ensg2hg.conv = load_hg_ensg()
        conv = ensg2hg.conv

    try:
        return conv[gene]
    except KeyError:
        if verbose:
            print("Gene ", gene, "missing from ensg2hg")
        return ""



def load_hg_prot():
    prot_hg_file = "raw_files/ProteinToGene.txt"

    prot_hg_df = pd.read_csv(prot_hg_file, sep="\t")
    genes = prot_hg_df["Gene_Name"]
    prots = prot_hg_df["Protein"]

    prot2hg = dict(list(zip(prots, genes)))

    hg2prot = dict(list(zip(genes, prots)))

    return hg2prot, prot2hg



def hg2prot(gene, verbose=False):
    try:
        conv = hg2prot.conv
    except AttributeError:
        hg2prot.conv, _ = load_hg_prot()
        conv = hg2prot.conv

    try:
        return conv[gene]
    except KeyError:
        if verbose:
            print("Gene ", gene, " missing from hg2prot")
        return ""

def prot2hg(gene, verbose=False):
    try:
        conv = prot2hg.conv
    except AttributeError:
        _, prot2hg.conv = load_hg_prot()
        conv = prot2hg.conv

    try:
        return conv[gene]
    except KeyError:
        if verbose:
            print("Gene ", gene, "missing from prot2hg")
        return ""

def load_syn_hg():
    synonym_file = "HGNC-to-Ensembl.txt"

    syn_df = pd.read_csv(synonym_file, sep="\t")


    syn2hg = {}

    prev_syns = []
    for i in np.where(pd.notnull(syn_df["Previous Symbols"]))[0]:
        syns = syn_df["Previous Symbols"][i].split(",")
        hg = syn_df["Approved Symbol"][i]
        for syn in syns:
            if syn not in syn2hg:
                syn2hg[syn] = {hg}
            else:
                syn2hg[syn].add(hg)
        prev_syns.extend(syns)

    just_syns = []
    for i in np.where(pd.notnull(syn_df["Synonyms"]))[0]:
        syns = syn_df["Synonyms"][i].split(",")
        hg = syn_df["Approved Symbol"][i]
        for syn in syns:
            if syn not in syn2hg:
                syn2hg[syn] = {hg}
            else:
                syn2hg[syn].add(hg)
        just_syns.extend(syns)

    extra_syns = [syn for syn in syn2hg if len(syn2hg[syn]) > 1]

    print("Num syns:", len(list(syn2hg.keys())))
    print("Num symbols: ", len(prev_syns), "set: ", len(set(prev_syns)))
    print("Num Synonyms: ", len(just_syns), "set: ", len(set(just_syns)))
    print("Num more than one: ", len(extra_syns))

    unique_syn2hg = syn2hg.copy()
    for syn in syn2hg:
        if len(syn2hg[syn]) == 1:
            unique_syn2hg[syn] = list(syn2hg[syn])[0]
        else:
            del unique_syn2hg[syn]

    return syn2hg, unique_syn2hg

def load_hg_gr():
    gr_file = "raw_files/Glucocorticoid_receptor_regulatory_network-gene-conversions.txt"

    df = pd.read_csv(gr_file, sep="\t")

    df = df[df["Gene Symbol"] != ""]

    gr = df["GR Gene"].values
    genes = df["Gene Symbol"].values

    gene2gr = dict(list(zip(genes, gr)))
    gr2gene = dict(list(zip(gr, genes)))


    return gene2gr, gr2gene

def gr2hg(gene, verbose=False):
    try:
        conv = gr2hg.conv
    except AttributeError:
        _, gr2hg.conv = load_hg_gr()
        conv = gr2hg.conv

    try:
        return conv[gene]
    except KeyError:
        if verbose:
            print("Gene ", gene, "missing from gr2hg")
        return ""


def load_genes(filename, header=False, verbose=True, as_set=True):
    genes = []
    with open(filename, 'rU') as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        if header:
            next(reader)
        for row in reader:
            genes.append(row[0])

    if as_set:
        genes = set(genes)

    if verbose:
        print("# Genes in ", filename, ": ", len(genes))
    return genes

def write_genes(filename, genes):
    with open(filename, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        for gene in genes:
            writer.writerow([gene])
    print(len(genes), " genes writen to: ", filename)



def annotate_cols(df, cols, genes, name):
    """Return the sum across a row's cols of genes in the gene"""

    df[name] = np.sum([df[col].isin(genes).astype(int) for col in cols], axis=0)
    return df

def make_pair_col(df, cols, col_name = "Pair", reorder=False):
    assert len(cols) == 2
    pair_array = [str(tuple(sorted(x))) for x in df[cols].values]
    df[col_name] = pair_array
    if reorder:
        df[cols[0]] = [eval(x)[0] for x in pair_array]
        df[cols[1]] = [eval(x)[1] for x  in pair_array]
    return df

def make_cause_col(df, cause_col, effect_col, col_name = "Cause-Effect"):

    df[col_name] = [tuple(x) for x in zip(df[cause_col],df[effect_col])]

    return df

def make_type_col(df, col_name):
    """
    :param df:
    :param name:
    :param col_name:
    :return: Add a new column with a column of 1s indicating the type of the column.
    """

    df[col_name] = 1

    return df



def annotate_cols_summary(df, cols, col_name = "Type", do_reverse=False):
    """Assumes there is a non-zero, non-null in each of the cols. Summarizes them one-by-one"""
    assert len(set(cols).intersection(set(df.columns.values))) == len(cols)

    sum_df = df.copy()

    type_array = np.empty(len(df), dtype=object)

    type_array[:] = ""

    rev_cols = cols[:]
    if do_reverse:
        rev_cols.reverse()

    for col in rev_cols:
        indices = np.where(np.logical_and(pd.notnull(sum_df[col]), sum_df[col] != 0))[0]
        print("FOr col", col, ":", len(indices))
        type_array[indices] += col + ","

    sum_df[col_name] = type_array

    #print sum_df[sum_df[col_name] != ""][col_name]
    return sum_df


def get_annotated(df, names=[], cols=genecols,
                  suffixes = ["_A", "_B"]):
    """  Given annotation type
    :param df:
    :param name:
    :return:  df with annotations, and with a col of all rows containing at least one
    """

    file_dir                                    =   "/Users/jlu96/v-causal-snps/data/GeneExpressionData/GGR_Network/raw_files/"
    name2filename = {}
    name2filename["Diabetes_highconf"]          =  file_dir + "GSEA-diabetes-high-conf-7_29_16-genes.txt"
    name2filename["Diabetes_lowconf"]           =  file_dir + "GSEA-diabetes-low-conf-7_29_16-genes.txt"
    name2filename["Metabolism_GSEA_highconf"]   =  file_dir + "GSEA-gluconeo-glyc-gluc-metab-8_2_16-genes.txt"
    name2filename["Metabolism_GSEA_lowconf"]    =  file_dir + "GSEA-gluconeo-am-glyc-lip-gluc-metab-7_29_16-genes.txt"
    name2filename["Immune"]                     =  file_dir + "immune-GSEA-GOC-genes.txt"
    name2filename["Immune_GOC-Priority"]        =  file_dir + "Immune-GOC-Priority.txt"
    name2filename["Immune_GSEA"]                =  file_dir + "GSEA-immune-inflammatory-7_29_16-genes.txt"
    name2filename["GR"]                         =  file_dir + "GR_GO_direct_candidates_union.txt"
    name2filename["GR_direct-candidate"]        =  file_dir + "GR_direct_candidates-HGNC.txt"
    name2filename["GR_GO"]                      =  file_dir + "GR_Pathway-GO-HGNC.txt"
    name2filename["DEG_edgeR-0.05"]             =  file_dir + "sig_genes_reg_fdr-0.05-all.txt"

    default_loads = ["GR", "DEG_edgeR-0.05", "Diabetes_highconf", "Diabetes_lowconf", "Immune", "Metabolism_GSEA_highconf",
                     "Metabolism_GSEA_lowconf"]

    if names == []:
        names = default_loads
    else:
        for name in names:
            if name not in name2filename:
                raise ValueError("Name " + name + " not in list of annotations")


    df_genes = get_genes(df, genecols=cols)

    for name in names:
        filename = name2filename[name]

        annot_genes = load_genes(filename)

        # Get the df's genes. See # annotated

        both_genes = df_genes.intersection(annot_genes)

        print("# ", name, " genes in df: ", len(both_genes))


        newcols = []
        # Annotate each column if extra genes
        for col, suffix in zip(cols, suffixes):
            newcol = name + suffix
            newcols.append(newcol)
            df = annotate_cols(df, [col], annot_genes, newcol)


        # Annotate the total number in that row
        df[name] = np.sum([df[newcol] for newcol in newcols], axis=0)

        print("# ", name, " edges in df: ", len(np.where(df[name])[0]))

    return df






def load_network_df(filename, cols=genecols, make_pairs=True):
    df = pd.read_csv(filename, sep="\t")

    df = make_pair_col(df, cols)

    print("Initial pairs: ", len(df))
    df.drop_duplicates(subset="Pair", keep="first", inplace=True)
    print("After dup drops: ", len(df))

    # Drop self-interactions

    df.drop(df[cols[0]] == df[cols[1]])

    print("After self drops: ", len(df))

    df.index = list(range(len(df)))

    return df

def load_causal_df(filename, cause_col, effect_col):
    cols= [cause_col, effect_col]
    df = pd.read_csv(filename, sep="\t")

    df = make_pair_col(df, cols, reorder=False)

    df = make_cause_col(df, cause_col, effect_col)

    #drop self-interactions
    df.drop(df[cols[0]] == df[cols[1]])

    df.index = list(range(len(df)))

    return df

def count_rows(df, name):
    return len(np.where(np.logical_and(pd.notnull(df[name]), df[name] != 0))[0])

def get_genes(df, genecols=genecols):
    genes = set()
    for genecol in genecols:
        genes.update(set(df[genecol].values))
    return genes

def get_genes_in(df, name, genecols=["Gene"], verbose=False):
    """
    :param df:
    :param name:
    :return: Genes where the colname is not empty
    """

    in_df = df[np.logical_and(pd.notnull(df[name]), df[name] != 0)]

    genes = get_genes(in_df, genecols=genecols)

    if verbose:
        print("Genes in ", name, ":", len(genes))

    return genes

def get_cause_to_effect(pairs):
    cause2effects = {}
    for pair in pairs:
        cause, effect = pair
        if cause not in cause2effects:
            cause2effects[cause] = set()
        cause2effects[cause].add(effect)

    return cause2effects

def filter_pairs(pairs, genes):
    """
    :param pairs:
    :param genes:
    :return: Limit only to pairsi n the genes
    """
    return [p for p in pairs if p[0] in genes and p[1] in genes]


# def get_cause_plot_triples(cause2effects, sort_dict=None):
#     """
#     :param cause2effects: Dictionary of the causes and effects
#     :param sort_dict: Dictionary returning a key to sort the effects by
#     :return: a list of plot_triples, cause at beginning
#     """
#     plot_triples_list = []
#     for cause in cause2effects:
#         effects = sorted(cause2effects[cause], key = lambda entry: sort_dict[entry], reverse=True)
#         effect_list = pj.partition_inputs(list(effects), int(round(len(effects)/2.0)))
#
#
#         plot_triples_list.extend([[cause] + e for e in effect_list])
#
#     print "Plot triples: "
#     print plot_triples_list[0:20]
#
#     return plot_triples_list

def limit_to_genes_all(df, genes, cols=genecols):
    """
    :param df:
    :param genes:
    :param cols:
    :return: df where all values in cols are in genes
    """

    num_cols = len(cols)
    indices = np.sum([df[col].isin(genes).values for col in cols], axis=0) >= num_cols

    new_df = df[indices].copy()

    new_df.index = list(range(len(new_df)))

    return new_df


# def df_to_graph_causal(df, key_col, cause_col, source_col=None, target_col=None, type_cols=[]):
#     """
#     :param df: Input dataframe
#     :param key_col: Column containing source-target pairs
#     :param cause_col: Column saying the type of causal relation. If None, assume just PPI
#     :param source_col: The source column of the causal relation if it exists
#     :param target_col: Target column
#     :param type_cols: Other attributes to annotate the edge with
#     :return: A Digraph where each edge has attributes: source, target (None if Causal Type is None)
#     Causal Type, and other type_col annotations
#     """
#     G = nx.Graph()
#
#     for i in range(len(df)):
#         if pd.notnull(cause_col):
#
#
#         type_dict = {}
#         for type_col in type_cols:
#             type_dict[type_col] = df[type_col][i]
#
#         G.add_edge(source, target, attr_dict=type_dict)
#
#     return G


def matr_to_net(matr_df, edge_name=None, abs_name=None, cause_effect_col = "Cause-Effect", colnames=None, make_pair=False,
                make_type=True, name=None, sort_by=None, extra_dict=None,
                no_abs=False, do_sort=True):
    """
    Convert a coefficient matrix to a network.

    :param matr_df: rows and columns are the genes
    :param edge_name: The name to give the column of edge values (from matrix)
    :param extra_dict: Dictionary to update the rest with.


    :param cause_effect_col:
    :param colnames: Customize cause effect colnames
    :param make_pair: Make a pair column?
    :param extra_dict: an extra dictionary of attributes you want to specify


    :return: net_df, the network from all matrix nonzero entries
    """
    if colnames == None:
        colnames = ["Cause", "Effect"]

    if edge_name == None:
        edge_name = "Weight"

    if abs_name == None:
        abs_name = "AbsWeight"

    if sort_by == None:
        sort_by = abs_name

    matr = matr_df.values

    genes = matr_df.columns.values
    indices = np.where(matr != 0)
    betas = matr[indices]


    net_dict = collections.OrderedDict()

    net_dict[cause_effect_col] = ["-".join(x) for x in zip(genes[indices[0]],genes[indices[1]])]
    net_dict[colnames[0]] = genes[indices[0]]
    net_dict[colnames[1]] = genes[indices[1]]
    net_dict[edge_name] = betas

    if not no_abs:
        net_dict[abs_name] = np.absolute(betas)

    if extra_dict != None:
        net_dict.update(extra_dict)


    net_df = pd.DataFrame.from_dict(net_dict)

    if make_pair:
        net_df = make_pair_col(net_df, colnames)
    if make_type:
        net_df["Type"] = name

    if do_sort:
        net_df.sort_values(sort_by, ascending=False, inplace=True)
    print("New network (edges, attributes) = ", net_df.shape)

    return net_df


def matr_file_to_net_file(matr_file, name, net_file=None, conv_to_hg=True, add_pair=True):
    """Convert a matrix file to a network file"""

    if not net_file:
        net_file = matr_file[:-4] + "-network.txt"

    print(name)

    cause_name = name + " Cause"
    effect_name = name + " Effect"

    matr_df = pd.read_csv(matr_file, sep="\t", header=0, index_col=0)

    print(matr_df.head())

    net_df = matr_to_net(matr_df, name, colnames=[cause_name, effect_name])

    print(net_df.head())

    if conv_to_hg:
        net_df[cause_name] = [ensg2hg(gene) for gene in net_df[cause_name]]
        net_df[effect_name] = [ensg2hg(gene) for gene in net_df[effect_name]]

        print("Post conversion: ")
        print(net_df.head())

    if add_pair:
        net_df = make_pair_col(net_df, [cause_name, effect_name], "Pair")

    print()
    print("FINAL:")
    print(net_df.head())
    print("Writing to ", net_file)
    net_df.to_csv(net_file, sep="\t", index=False)

def overlay_dfs(old_df, over_df, key = "Pair", over_cols=[], fill_empty=False, fill_genecols=genecols, how='outer'):
    """
    Overlay dfs where the key is Pair and using the known genecols.
    :param old_df:
    :param over_df: Df to overlay
    :param key: Key to use to match up rows. Unwrap this into the genecols
    :param over_cols: The columns to merge over
    :param cols: Columns to overlay with
    :return:
    """
    if len(genecols) != 2:
        raise ValueError("There must be 2 genecols to unwrap the pair.")

    if over_cols == []:
        over_cols = over_df.columns.values


    # columns to merge over
    add_df = over_df[over_cols]

    print("Over cols: ", over_cols)

    #
    df = old_df.merge(add_df, left_on=key, right_on=key, how=how)

    if fill_empty:
        df = fill_empty_genecols(df, fill_genecols, key)

    return df






def fill_empty_genecols(df, genecols=genecols, key="Pair"):
        # for the places where genecols are empty, rewrite the annotation
    indices = np.where(np.logical_or(pd.isnull(df[genecols[0]]), pd.isnull(df[genecols[1]])))[0]
    print("Num missing genecols: ", len(indices))
    print(indices)

    pairs = [eval(p) for p in df[key][indices].values]

    zip_cols = list(zip(*pairs))
    df[genecols[0]][indices] = zip_cols[0]
    df[genecols[1]][indices] = zip_cols[1]

    return df

def graph_add_df_both_ways(G, df, col1, col2, edge_cols):

    edge_set = set(edge_cols)
    for i, r in df.iterrows():
        gene1 = r[col1]
        gene2 = r[col2]

        edge_dict = dict((x, r[x]) for x in edge_set if x in r)

        G.add_edge(gene1, gene2, attr_dict=edge_dict)
        G.add_edge(gene2, gene1, attr_dict=edge_dict)

    return G

def graph_add_df_one_way(G, df, col1, col2, edge_cols):

    edge_set = set(edge_cols)
    for i, r in df.iterrows():
        gene1 = r[col1]
        gene2 = r[col2]

        edge_dict = dict((x, r[x]) for x in edge_set if x in r)

        G.add_edge(gene1, gene2, attr_dict=edge_dict)

    return G

def get_feedforward(G, Ps=None, Xs=None, Ts=None):
    """
    :param G: A Digraph
    :return: a list of tuples (P, X, C) where P -> X, X -> T, P -> T
    """

    if Ps == None:
        Ps = set(G.nodes())
    if Xs == None:
        Xs = set(G.nodes())
    if Ts == None:
        Ts = set(G.nodes())

    feedforward_set = set()

    for X in Xs:
        for T in set(G.successors(X)).intersection(Ts):
            this_Ps = set(G.predecessors(X)).intersection(set(G.predecessors(T)).intersection(Ps))

            if len(this_Ps) > 0:
                for P in this_Ps:
                    feedforward_set.add((P, X, T))

    print("Num feedforward: ", len(feedforward_set))
    return feedforward_set