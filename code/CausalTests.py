__author__ = 'jlu96'

from sklearn.linear_model import LinearRegression, ElasticNetCV, ElasticNet
from sklearn.cross_validation import StratifiedKFold


import scipy.stats as stats
import numpy as np
import pandas as pd
import csv



########### PAIRWISE GRANGER CAUSALITY



def load_kwargs_file(argsfile):
    with open(argsfile, 'rU') as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        newdict = {}
        for row in reader:
            arg = row["Argument"]
            value = row["Value"]
            try:
                newdict[arg] = eval(value)
            except NameError:
                newdict[arg] = value

    return newdict



def get_lagged_vector(X_list, index, tau, E):
    """Generate lagged vector of X: (X[index], X[index - tau], ..., X[index - tau*E - 1]).

    :param X_list: list of X observations at equally spaced time intervals.
    :param index: index to start from
    :param tau: time lag
    :param E: dimension of lagged vector

    :return: (X[index], X[index - tau], ..., X[index - tau*E - 1])
    """
    assert index >= tau * (E - 1)
    return [X_list[index - e * tau] for e in range(E)]


def pairwise_granger_causality(i, j, model_order):
    """
    Computes p-value of fit to determine if i granger causes j.
    i and j must have length at least model_order (to construct lagged vectors)
    Returns: pvalue
    """
    if len(i) != len(j):
        raise ValueError("i and j must have same length")

    m = len(i) - model_order
    if m <= 2 * model_order:
        raise ValueError("i and j must have length > 2*model_order + model_order, so that target has legnth 2model_order" +
                "with 2model_order predictors")

    T = len(i)
    #lagged_i goes from p up to T, the length of i. It is a T-p matrix
    lagged_i = np.array([get_lagged_vector(i, index, 1, model_order) for index in range(model_order-1, T-1)])

    lagged_j = np.array([get_lagged_vector(j, index, 1, model_order) for index in range(model_order-1, T-1)])

    target = np.array(j[model_order:T])


    # Restricted classification
    clf_j = LinearRegression()
    clf_j.fit(lagged_j, target)
    pred_j = clf_j.predict(lagged_j)

    SSE_j = np.sum(np.power(target - pred_j, 2))

    # Unrestricted classification
    clf_ij = LinearRegression()
    lagged_ij = np.concatenate((lagged_i, lagged_j), axis=1)
    clf_ij.fit(lagged_ij, target)
    pred_ij = clf_ij.predict(lagged_ij)

    SSE_ij = np.sum(np.power(target - pred_ij, 2))


    beta_i = clf_ij.coef_[0:lagged_i.shape[1]]

    # Formula is taken from: http://pages.uoregon.edu/aarong/teaching/G4075_Outline/node4.html
    # The variance explained by i itself compared to variance of just j
    # IF i explains more variance then we reject
    #DOF of first is p = 2p - p = # in ij fit - # in j fit (# of params added)
    # DOF of bottom is is n - p - 1 (where p is the number of parameters in the bottom model)
    F = ((SSE_j - SSE_ij) * 1.0 / model_order)/(SSE_ij/(len(target) - 2 * model_order - 1))

    p_value = 1.0 - stats.f.cdf(F, model_order, len(target) - 2 * model_order - 1)

    return p_value, beta_i


def pairwise_granger_causality_row(matr1, matr2, rows, model_order):
    """
    :param matr1: Cause matrix, n x T (n genes, T timepoients)
    :param matr2: Effect matrix, n x T (n genes, T timepoints)
    :param rows: The rows of the effect matrix to test
    :param model_order: Model order of granger causality
    :return: betas, a n x len(rows) x model_order matrix, where betas_ij_ is causal coefficient of gene i onto gene j.
    """
    n, T = np.array(matr1).shape
    assert (np.array(matr1).shape == np.array(matr2).shape)


    betas = np.zeros(shape=(n, len(rows), model_order))

    # row is the effect gene
    for j, row in zip(range(len(rows)), rows):
        # i is the cause gene
        for i in range(n):
            _, betas[i][j] = pairwise_granger_causality(matr1[i], matr2[row], model_order)

    return betas








########### ELASTIC NET



def get_Xy_lagged(matr1, matr2, row, lag):
    """
    :param matr1: Input matrix, n genes x T timepoints
    :param matr2: Effect matrix, n genes x T timeponts
    :param row: row of matr 2 to use. Will replace corresponding in matr1
    :param lag: lag
    :return: X_t (T - lag) x (n * lag) matrix of predictors
            y_t (T - lag) matrix of outputs
    """
    assert (np.array(matr1).shape == np.array(matr2).shape)

    n, T = np.array(matr1).shape


    y = matr2[row].T

    X = np.concatenate((matr1[:row], [matr2[row]], matr1[row+1:])).T

    # X is T x n
    # y is T x 1
    # Generate lagged vectors, X_t, y_t
    # X_t is (T - lag) x (n*lag)
    # y_t is (T - lag) x 1


    # X is a matrix of T - lag rows where for each row j, i = j + lag-1
    # [A_i, B_i, .... Z_i, A_{i-1}, B_{i-1}, .... Z_{i-1}......  A_{i-lag+1}, B_{i-lag+1}, .... Z_{i-lag+1}]
    # Where A, B, Z are genes



    X_t = np.zeros((T - lag, n * lag))

    for j in range(T - lag):
        for k in range(lag):
            X_t[j, n*k:n*(k+1)] = X[j + lag - k - 1, ]

    # Y is a vector of T -lag rows where for row j, the value of j is i = j + lag
    y_t = y[lag:T]


    return X_t, y_t


def fit_Xy(matr1, matr2, beta, lag):
    """
    NOTE: ADJUST SO CAN BE CALLED BY JUST 1 MATR 2 OR MULTIPLE ROWS
    :param matr1: Cause matrix, n x T timepoints
    :param matr2: Effect matrix, n x T timepoint
    :param beta: Matrix of coefficients
    :param lag: Output matrix
    :return: Y_test, Y_pred, r^2, Durbin-watson fits for each output y
    """



def enet_granger_causality_row_cv(matr1, matr2, rows,  model_order=2, min_alpha=0.01, max_alpha=0.99, num_alpha=20, min_lambda=0.01, max_lambda=1, num_lambda=None, choose_by="test_err", choose_bottom=True, top_num=None, top_perc=10, n_folds=3,
                               max_iter=1000000):
    """
    :param matr1: Cause matrix
    :param matr2: Effect matrix
    :param rows: Which rows (i.e. genes) to test as effects
    :param min_alpha: value of alpha to use for regularization
    :param max_alpha:
    :param num_alpha:
    :param choose_by: Which column to use to select from test_df
    :param choose_bottom: Choose the smallest of choose_by. If false, choose largest
    :param top_num: The top several from cross-validation
    :param top_perc: The top % of params from cross-validation
    :param n_folds: # Cross-validation folds
    :param model_order: # model order of fit
    :return: beta_tuple  (ngenes x nrows x model_order) (where the t-1 coefficient is first, t-2 is econd, t-3 is 3rd...), all_res_df (all cvs), use_df (df of the one that was used)
    """

    assert (np.array(matr1).shape == np.array(matr2).shape)

    n, T = np.array(matr1).shape
    lag = model_order


    # Each row is the betas for that effect gene.
    betas = np.zeros(shape=(len(rows), n * model_order))


    # Initialize the error and other value storage later
    header = None

    alphas = np.linspace(min_alpha, max_alpha, num_alpha)

    if num_lambda != None:
        lambdas = np.linspace(min_lambda, max_lambda, num_lambda)
        print "Lambdas are ", lambdas
    else:
        lambdas = None

    for i, row in zip(range(len(rows)), rows):

        X_t, y_t = get_Xy_lagged(matr1, matr2, row, lag)

        # Generate Cross-validation

        cv = get_TS_cv(len(y_t), n_folds)

        print "CV is ", cv

        # CHECK CORRECT ARGS/RETURN VALUES
        top_df = enet_granger_causality_cv(X_t, y_t, cv, alphas, lambdas=lambdas,top_num=top_num, top_perc=top_perc, max_iter=max_iter )

        res_df, test_betas = enet_granger_causality_test(X_t, y_t, top_df, max_iter=max_iter)

        res_df["Row"] = row

        if choose_bottom:
            index = np.where(res_df[choose_by] <= min(res_df[choose_by]))[0][0]
        else:
            index = np.where(res_df[choose_by] >= max(res_df[choose_by]))[0][0]

        # Add the best betas
        betas[i] = test_betas[index]

        if header == None:
            header = res_df.columns.values

            all_res_index = range(len(res_df) * len(rows))
            all_res_df = pd.DataFrame(index=all_res_index, columns=header)

            use_index = range(len(rows))
            use_df = pd.DataFrame(index=use_index, columns=header)


        all_res_df.iloc[len(res_df) * i : len(res_df) * (i + 1)] = res_df.sort_values(choose_by, ascending=choose_bottom).as_matrix()
        use_df.iloc[i] = res_df.iloc[index]

    beta_tuple = np.swapaxes(np.dstack(tuple(np.split(betas, lag, axis=1))), 0, 1)

    assert beta_tuple.shape == (n, len(rows), model_order)
    print "Shape of beta tuple is ", beta_tuple.shape

    return beta_tuple, all_res_df, use_df


def enet_granger_causality_row_load(matr1, matr2, rows, param_df, **kwargs):
    """
    :param matr1: Cause matrix
    :param matr2: Effect matrix
    :param rows:
    :param param_df: parameters
    :param model_order:
    :param max_iter:
    :return:    # beta_tuple (n x nrows x model_order)
    """

    if "choose_by" in kwargs:
        choose_by = kwargs["choose_by"]
    else:
        choose_by = "test_err"

    if "choose_bottom" in kwargs:
        choose_bottom = kwargs["choose_bottom"]
    else:
        choose_bottom = True

    if "model_order" in kwargs:
        model_order = kwargs["model_order"]
    else:
        model_order = 2

    if "max_iter" in kwargs:
        max_iter = kwargs["max_iter"]
    else:
        max_iter = 1000000


    assert (np.array(matr1).shape == np.array(matr2).shape)

    n, T = np.array(matr1).shape
    lag = model_order






    # Each row is the betas for that gene.
    betas = np.zeros(shape=(len(rows), n * model_order))


    # Initialize the error and other value storage later
    header = None


    for i, row in zip(range(len(rows)), rows):

        X_t, y_t = get_Xy_lagged(matr1, matr2, row, lag)

        # Generate Cross-validation

        use_param_df = param_df[param_df["Row"] == row]

        res_df, test_betas = enet_granger_causality_test(X_t, y_t, use_param_df, max_iter=max_iter)


        if choose_bottom:
            index = np.where(res_df[choose_by] <= min(res_df[choose_by]))[0][0]
        else:
            index = np.where(res_df[choose_by] >= max(res_df[choose_by]))[0][0]

        # Add the best betas
        betas[i] = test_betas[index]

        if header == None:
            header = res_df.columns.values

            all_res_index = range(len(res_df) * len(rows))
            all_res_df = pd.DataFrame(index=all_res_index, columns=header)

            use_index = range(len(rows))
            use_df = pd.DataFrame(index=use_index, columns=header)


        all_res_df.iloc[len(res_df) * i : len(res_df) * (i + 1)] = res_df.sort_values(choose_by, ascending=choose_bottom).as_matrix()
        use_df.iloc[i] = res_df.iloc[index]


    beta_tuple = np.swapaxes(np.dstack(tuple(np.split(betas, lag, axis=1))), 0, 1)

    assert beta_tuple.shape == (n, len(rows), model_order)
    print "Shape of beta tuple is ", beta_tuple.shape

    return beta_tuple, all_res_df, use_df





def enet_granger_causality_cv(X_t, y_t, cv, alphas, top_num=None, top_perc=4,max_iter=100, lambdas=None):
    """
    :param X_t: Lagged matrix: (T - lag ) x (n * lag)
    :param y_t: Output vector: (T - lag) x 1
    :param cv: The sklearn cross-validation generator
    :param alphas:
    :param top_num: # of top parameter results to keep
    :param top_perc: % of top param results to keep, only use if top_num = None
    :param max_iter: max iter in Enet fit
    :param lambdas:
    :return:
    """


    # alpha is the l1_ratio in sklearn
    if lambdas != None:
        use_lambdas = np.tile(lambdas, len(alphas)).reshape(len(alphas), len(lambdas))
        enet = ElasticNetCV(l1_ratio=alphas, alphas=use_lambdas, cv=cv, max_iter=max_iter)
        fit = enet.fit(X_t, y_t)

        use_lambdas = fit.alphas_
        use_lambdas = np.tile(use_lambdas, len(alphas)).reshape(len(alphas), len(lambdas))
        print "Used lambdas"
        print use_lambdas

    else:
        enet = ElasticNetCV(l1_ratio=alphas,  cv=cv, max_iter=max_iter)
        fit  = enet.fit(X_t, y_t)
        use_lambdas = fit.alphas_


    # lambdas is a matrix

    cv_mses = enet.mse_path_.sum(axis=2).flatten()


    cv_alphas = np.repeat(alphas, use_lambdas.shape[1])
    cv_lambdas = use_lambdas.flatten()

    if top_num == None:
        print "Num cv alphas: ", len(cv_alphas)

        top_num = int(len(cv_alphas) * top_perc / 100.0)
        print "Top num ", top_num

    # this will keep the smallest
    top_indices, top_mses = get_min_k(cv_mses, top_num)

    top_lambdas = cv_lambdas[top_indices]
    top_alphas = cv_alphas[top_indices]

    top_df = pd.DataFrame(data={"lambda.min": top_lambdas, "alpha": top_alphas, "error.min": top_mses})

    return top_df

def enet_granger_causality_test(X_t, y_t, top_df, max_iter=10000000):
    """
    Return the cv-parameters tested across the whole data
    :param X_t:
    :param y_t:
    :param top_df:
    :return: res_df, test_betas
    """

    test_errs = np.zeros(len(top_df))
    scores = np.zeros(len(top_df))
    dfs = np.zeros(len(top_df))

    test_coefs = np.zeros((len(top_df), X_t.shape[1]))
    for i in range(len(top_df)):
        alpha = top_df.iloc[i]["alpha"]
        lambda_min = top_df.iloc[i]["lambda.min"]
        enet = ElasticNet(l1_ratio=alpha, alpha=lambda_min, max_iter=max_iter)
        enet.fit(X_t, y_t)
        y_pred = enet.predict(X_t)
        test_errs[i] = np.average((y_t - y_pred)**2)
        scores[i] = enet.score(X_t, y_t)
        test_coefs[i] = enet.coef_

        dfs[i] = len(np.where(enet.coef_)[0])

    top_df["test_err"] = test_errs
    top_df["score"] = scores
    top_df["df"] = dfs


    return top_df, test_coefs


def get_TS_cv(length, n_folds):
    """
    :param length: Length of TS
    :param n_folds: # of folds to produce
    :return: indices where 1 in every n_folds is a test index
    """

    labels = [i/n_folds for i in range(length)]
    cv = StratifiedKFold(labels, n_folds)

    return cv

def get_min_k(array, k):
    """
    :param array: A numpyarray of numerics
    :return: topk_indices, topk_array
    the sorted top (or bottom) k of indices
    """

    a = np.array(array)
    topk_indices = sorted(np.argpartition(a, k)[:k], key = lambda entry: a[entry])

    return topk_indices, a[topk_indices]

















# ***********************************************************************
# THE BELOW VERSIONS OF PAIRWISE GRANGER CAUSALITY OR CCM ARE NO LONGER USED!!!!
# ***********************************************************************


#
# def main():
#     tstart = time.time()
#
#
#     input_file = args.input_file
#     out_file_prefix = args.out_file_prefix
#
#
#     start_index = args.start_index
#     end_index = args.end_index
#
#
#     df = gtm.load_file_and_avg(input_file)
#
#     genes = df['gene'][start_index:end_index].values
#
#     found_genes, geneTS = gtm.get_gene_TS(df, genes)
#
#     cause_type = args.cause_type
#
#     if cause_type == 'g':
#         model_orders = range(args.model_order_min, args.model_order_max + 1)
#
#         threshold = args.p_threshold
#
#         p_matr_list = []
#         sig_matr_list = []
#
#         for model_order in model_orders:
#             t_gc = time.time()
#             p_matr = pairwise_granger_causality_all(geneTS, model_order=model_order, use_processes=args.use_processes, procnum=args.procnum)
#             print "Time for granger causality", time.time() - t_gc
#
#
#             sig_matr = p_matr < threshold
#
#             p_matr_list.append(p_matr)
#             sig_matr_list.append(sig_matr)
#
#
#
#         all_sig_matr, all_sig_num, not_sig_num = gtm.compare_sig_matr(sig_matr_list=sig_matr_list)
#
#         print "Total number of significant pairs ", all_sig_num + not_sig_num
#         print "Pairs significant across all matrices ", all_sig_num, all_sig_num * 1.0 / (all_sig_num + not_sig_num)
#
#
#         out_file_name = out_file_prefix + "_GC.p"
#         pickle.dump([model_orders, p_matr_list, sig_matr_list, (all_sig_matr, all_sig_num, not_sig_num)], open(out_file_name, "w"))
#
#         print "Results written  to", out_file_name
#
#
#
#     # compare the significant matrices
#
#     # save the output p matrices
#
#     print "Total time used ", time.time() - tstart
#
# if __name__ == '__main__':
#     main()





# def pairwise_granger_causality_all(matr, pairs, model_order=2, use_processes=False, procnum=32,
#                                    add_ones=False):
#     """
#     Assume matr is an n by T matrix, where n is number of genes, T is # timepoints.
#     """
#     n, T = np.array(matr).shape
#
#     p_matr = np.zeros(shape=(n,n))
#
#     beta_matr = np.zeros(shape=(n,n,model_order))
#
#     if pairs == None:
#        pairs = itertools.permutations(range(n), 2)
#
#     if use_processes:
#         pairs = list(pairs)
#
#         print pairs[0:10], pairs[-10:]
#
#         function = pairwise_granger_causality_process
#         args = [matr, model_order, pairs]
#         input_list = pairs
#         input_index = 2
#         partition_input_function = pac.partition_inputs
#         join_functions_dict = {0: join_pvalue_matrices, 1: join_pvalue_matrices}
#
#         p_matr, beta_matr = pac.parallel_compute_new(function, args, input_list, input_index, partition_input_function,
#                                            join_functions_dict, number=procnum, multi_return_values=True)
#
#
#
#
#
#     else:
#         for pair in pairs:
#             i, j = pair
#             p_matr[i][j], beta_matr[i][j] = pairwise_granger_causality(matr[i], matr[j], model_order)
#
#     if add_ones:
#         p_matr += np.diagflat(np.ones(p_matr.shape[0]))
#
#     return (p_matr, beta_matr)
#
# def pairwise_granger_causality_process(matr, model_order, pairs):
#     """Return a matrix of zeros with the causal p-values return for the pairs.
#     """
#
#     n, T = np.array(matr).shape
#
#     p_matr = np.zeros(shape=(n,n))
#     beta_matr = np.zeros(shape=(n,n,model_order))
#
#
#     for pair in pairs:
#         i, j = pair
#         p_matr[i][j], beta_matr[i][j] = pairwise_granger_causality(matr[i], matr[j], model_order)
#
#     return p_matr, beta_matr
#
# def join_pvalue_matrices(p_matr_a, p_matr_b):
#     return p_matr_a + p_matr_b
#
#
# def pairwise_granger_causality_all_list(matr, pairs, model_order=2, use_processes=False, procnum=32):
#     """
#     Assume matr is an n by T matrix, where n is number of genes, T is # timepoints.
#     """
#     n, T = np.array(matr).shape
#
#     ps = np.zeros(len(pairs))
#
#     betas = np.zeros(shape=(len(pairs), model_order))
#
#     for index, pair in zip(range(len(pairs)), pairs):
#         i, j = pair
#         ps[index], betas[index] = pairwise_granger_causality(matr[i], matr[j], model_order)
#
#     return (ps, betas)
#
#
# def pairwise_granger_causality_matr_pair(matr1, matr2, pairs=None, model_order=2, use_processes=False, procnum=32):
#     """
#     Assume both matrices are n by T matrices, where n is number of genes, T is # timepoints.
#     """
#     n, T = np.array(matr1).shape
#     assert (np.array(matr1).shape == np.array(matr2).shape)
#
#     if pairs == None:
#         ind = range(n)
#         ind1 = np.repeat(ind, len(ind))
#         ind2 = np.tile(ind, len(ind))
#         pairs = zip(ind1, ind2)
#
#     ps = np.zeros(len(pairs))
#
#     betas = np.zeros(shape=(len(pairs), model_order))
#
#     for index, pair in zip(range(len(pairs)), pairs):
#         i, j = pair
#         ps[index], betas[index] = pairwise_granger_causality(matr1[i], matr2[j], model_order)
#
#     return (ps, betas)
#
#
# def get_CCM(X_list, Y_list, L, tau=1, E=3, test_indices=None, num_test=100, use_same=True):
#     """
#     Compute the correlation coeffiicent of using X's cross-mapped estimates to
#     predict Y at the test_indices.
#
#     :param X_list: list of X observations
#     :param Y_list: list of Y observations
#     :param L: Library length, i.e. number of X/Y observations to use to construct manifold for estiamtion
#     :param  tau: Time lag of lagged vecotr
#     :param E: dimesnion of lagged vector
#     :param test_indices: indices of Y to estimate. All shou be > than L. Default is random.
#     :param num_test: number of indeces to test if random
#     :param use_same: If we find an X that has the exact value of the test X, then use only those Xs to generate estiamtes of Y
#
#     :return: rho, the correlation coefficient of the estimates and true alues
#
#     Workflow
#     1 ) make all lagged vectors
#     2) assign first L - first_lag to train, in kdtree
#     3) predict by taking test_indices, convert
#     4) find closets tin lag tree
#     5) calculate weights
#     6) get indices in old
#     7) make estimate
#     """
#
#     length = len(X_list)
#     first_lag = tau * (E - 1)
#
#     train_indices = np.arange(first_lag, L)
#     other_indices = np.arange(L, length)
#     all_indices = np.arange(first_lag, length)
#     if test_indices == None:
#         test_indices = np.random.choice(other_indices, num_test, replace=False)
#
#     # make all lagged vectors
#     x_lag_list = np.array([get_lagged_vector(X_list, i, tau, E) for i in all_indices])
#     y_lag_list = np.array([get_lagged_vector(X_list, i, tau, E) for i in all_indices])
#
#
#     # put training X and Y (used for estimation) into kdtree for nearest neighbors seaerch
#     x_lag_train = x_lag_list[np.ix_(train_indices - first_lag)]
#
#     x_lagtree = cKDTree(x_lag_train, leafsize=100)
#
#
#
#
#     Y_target = Y_list[test_indices]
#     Y_ests = []
#
#     # generate each estimate
#     for k in range(len(test_indices)):
#         test_index = test_indices[k]
#         # for each t, find contemporaneous x[t]
#         x_lag = x_lag_list[test_index - first_lag]
#
#
#
#         # Find e+1 nearest neighbors, Calculate distances
#         distances, indices = x_lagtree.query(x_lag, k=E+1)
#         min_dist = min(distances)
#
#
#         # Case 1: we find an X that has the exact same value as the test X
#         # In this case, use only those Xs that have same value as test X, and take
#         # the average of Ys
#         if (use_same and min_dist == 0) or all([dist == 0 for dist in distances]):
#             zipped = zip(*[(dist, i) for dist, i in zip(distances, indices) if dist == 0])
#             distances, indices = zipped[0], zipped[1]
#             weights = [1.0 / len(distances)] * len(distances)
#             Y_est = sum([weight * Y_list[i + first_lag] for weight, i in zip(weights, indices)])
#
#         # Case 2: all Xs are different
#         # Use exponential weighting as in paper
#         else:
#             zipped = zip(*[(dist, i) for dist, i in zip(distances, indices) if dist != 0])
#             distances, indices = zipped[0], zipped[1]
#             min_dist = min(distances)
#
#             # Calculate weights
#             weights = np.array([np.exp(-dist * 1.0/ min_dist) for dist in distances])
#
#             weights /= sum(weights)
#
#             # generate estimate y^
#             Y_est = sum([weight * Y_list[i + first_lag] for weight, i in zip(weights, indices)])
#
#         Y_ests.append(Y_est)
#
#     rho = np.corrcoef(Y_ests, Y_target)[0][1]
#
#
#     return rho
#
#
# def get_CCM_complete(X_list, Y_list, Ls, tau=1, E=3, test_indices=None, num_test=100, use_same=True):
#     """
#     :param X_list: list of X observations
#     :param Y_list: list of Y observations
#     :param Ls: Library lengths, i.e. number of X/Y observations to use to construct manifold for estiamtion
#     :param  tau: Time lag of lagged vecotr
#     :param E: dimesnion of lagged vector
#     :param test_indices: indices of Y to estimate. All shou be > than L. Default is random.
#     :param num_test: number of indeces to test if random
#     :param use_same: If we find an X that has the exact value of the test X, then use only those Xs to generate estiamtes of Y
#
#     :return: rhos, the correlation coefficient for each value of L
#     """
#
#
#     rhos = [get_CCM(X_list, Y_list, L, tau=tau, E=E, test_indices=test_indices, num_test=num_test, use_same=use_same)
#            for L in Ls]
#
#     for i in range(len(rhos)):
#         if np.isnan(rhos[i]):
#             rhos[i] = 0
#
#     return rhos, Ls


