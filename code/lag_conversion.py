__author__ = 'jlu96'


import numpy as np


def conv_to_X(matr, lag, has_reps=False, verbose=False):
    """
    matr: n x T (x r ) matrix of gene expression

    return: X_t: (T - lagged) x (lag * n)

    Sample by row
    X_t is a matrix of T - lag replace_rows where for each row j, i = j + lag-1
    [A_i, B_i, .... Z_i, A_{i-1}, B_{i-1}, .... Z_{i-1}......  A_{i-lag+1}, B_{i-lag+1}, .... Z_{i-lag+1}]

    """

    n = matr.shape[0]
    T = matr.shape[1]

    # if has reps, just get the X_t for each individually, and concatenate them

    if has_reps:
        r = matr.shape[2]

        X_t = np.zeros(((T - lag) * r, n * lag))

        for i in range(r):
            X_t[(T - lag) * i: (T - lag) * (i+1), :] = conv_to_X(matr[:, :, i], lag, has_reps=False)


        return X_t


    if verbose:
        print("Shape of input is ", matr.shape)

    X = matr.T # T x n

    # Generate lagged vectors, X_t
    # X_t is (T - lag) x (n*lag)
    X_t = np.zeros((T - lag, n * lag))

    if verbose:
        print(X.shape)
        print("Shape of out is ", X_t.shape)

    # X_t is a matrix of T - lag replace_rows where for each row j, i = j + lag-1
    # [A_i, B_i, .... Z_i, A_{i-1}, B_{i-1}, .... Z_{i-1}......  A_{i-lag+1}, B_{i-lag+1}, .... Z_{i-lag+1}]
    # Where A, B, Z are genes
    for j in range(T - lag):
        for k in range(lag):
            # k = 0 -> first columns, X[j + lag - 1, ]
            # Range of X indices: j + lag - 1 to j. Decreasing. So the first is the i-1, second is i-2, etc.
            X_t[j, n*k:n*(k+1)] = X[j + lag - k - 1, ]

    return X_t


def conv_to_Y(matr, lag, has_reps=False, verbose=False):
    """
    matr: n x T ( x r) matrix of gene expression
    returns: Y_t: (T - lag)* r x n matrix of expression
    """

    n = matr.shape[0]
    T = matr.shape[1]

    # if has reps, just get the X_t for each individually, and concatenate them

    if has_reps:
        r = matr.shape[2]

        Y_t = np.zeros(((T - lag) * r, n ))

        for i in range(r):
            Y_t[(T - lag) * i: (T - lag) * (i+1), :] = conv_to_Y(matr[:, :, i], lag, has_reps=False)

        return Y_t

    Y = matr.T # T x n

    Y_t = np.zeros((T - lag, n))

    for i in range(n):
        Y_t[:, i] = Y[lag:T, i]


    if verbose:
        print("input shape", matr.shape)
        print("input", matr)
        print("Transpose shape", Y.shape)
        print("output shape", Y_t.shape)
        print("output", Y_t)
    return Y_t




def get_XY_lagged(X_matr, Y_matr, lag, replace_row, has_reps=False, verbose=False,
                 bootstrap=False, seed=None):
    """
    X_matr: n x T ( x r) of predictor genes, r is # replicates
    Y_matr: 1 x T of (x r) reponse genes
    lag:
    replace_row: row of X to replace by Y_matr
    return: X_t (r*(T - lag), n * lag) , Y_t (r*(T - lag), 1)
    """

    if X_matr.shape[1] != Y_matr.shape[1]:
        raise ValueError("X and Y must both have T timepoints")

    if Y_matr.shape[0] != 1:
        raise ValueError("Y_matr must have exactly 1 gene")

    if has_reps and X_matr.shape[2] != Y_matr.shape[2]:
        raise ValueError("X and Y must have same third axis size")



    n = X_matr.shape[0]
    T = X_matr.shape[1]

    wanted_rows = list(range(len(X_matr)))
    wanted_rows.remove(replace_row)


    if has_reps:
        new_X_matr = np.zeros((n, T, X_matr.shape[2]))
    else:
        new_X_matr = np.zeros((n, T))


    new_X_matr[replace_row] = Y_matr
    new_X_matr[wanted_rows] = X_matr[wanted_rows]

    X_t = conv_to_X(new_X_matr, lag, has_reps=has_reps)

    Y_t = conv_to_Y(Y_matr, lag, has_reps=has_reps)

    if bootstrap:
        if seed == None:
            raise ValueError("Bootstrap: random seed not set. Must be integer.")

        # print "Performing bootstrap resampling of the samples in lagged matrix."
        # print "Seed = ", seed

        np.random.seed(seed)
        ind = np.random.choice(X_t.shape[0], X_t.shape[0], replace=True)
        X_t = X_t[ind]
        Y_t = Y_t[ind]


    if verbose:
        print("X matr:", X_matr.shape)
        print(X_matr[0:5,])
        print("Y matr:", Y_matr.shape)
        print(Y_matr)
        print("new_X_matr:", new_X_matr.shape)
        print(new_X_matr[0:5,])
        print("rows removed: ", replace_row)
        print("X_t: ", X_t[:,])
        print("Y_t:", Y_t)

    return X_t, Y_t



# def old_get_XY_lagged(X_matr, Y_matr, lag, replace_rows=None, has_reps=False, verbose=False):
#     """
#     X_matr: m x T ( x r) of predictor genes, r is # replicates
#     Y_matr: 1 x T of (x r) reponse genes
#     lag:
#     replace_rows: rows of X to ignore due to prescence in Y_matr
#     Will put the Y at the very beginning here.
#     return: X_t (r*(T - lag), (m - rows + 1) * lag) , Y_t (r*(T - lag), 1)
#     """
#
#     if X_matr.shape[1] != Y_matr.shape[1]:
#         raise ValueError("X and Y must both have T timepoints")
#
#     if Y_matr.shape[0] != 1:
#         raise ValueError("Y_matr must have exactly 1 gene")
#
#     if has_reps and X_matr.shape[2] != Y_matr.shape[2]:
#         raise ValueError("X and Y must have same third axis size")
#
#
#
#
#     T = X_matr.shape[1]
#
#     wanted_rows = range(len(X_matr))
#
#     if replace_rows != None:
#         # The number of replace rows should be same as Y_matr here
#         assert len(replace_rows) == Y_matr.shape[0]
#
#         for row in replace_rows:
#             wanted_rows.remove(row)
#
#
#     if has_reps:
#         new_X_matr = np.zeros((len(wanted_rows) + 1, T, X_matr.shape[2]))
#         # LEFT OFF HERE -JLu 12/15/16
#     else:
#         new_X_matr = np.zeros((len(wanted_rows) + 1, T))
#
#     new_X_matr[0, ] = Y_matr
#     new_X_matr[1:, ] = X_matr[wanted_rows]
#
#     X_t = conv_to_X(new_X_matr, lag, has_reps=has_reps)
#
#     Y_t = conv_to_Y(Y_matr, lag, has_reps=has_reps)
#
#
#     if verbose:
#         print "X matr:", X_matr.shape
#         print X_matr[0:5,]
#         print "Y matr:", Y_matr.shape
#         print Y_matr
#         print "new_X_matr:", new_X_matr.shape
#         print new_X_matr[0:5,]
#         print "rows removed: ", replace_rows
#         print "X_t: ", X_t[:,]
#         print "Y_t:", Y_t
#
#     return X_t, Y_t



# obsolete since no coef formed anymore
# def old_align_coefs(coefs, lag, verbose=False):
#     """
#     coefs: n* lag x n. the i * lagth row of col k is the auto-coefficient of k on k
#
#
#     # align the coefs:
#     # 1) split into same lag
#     # 2) For each lag: for out gene k, take the first coef and move to position k
#
#     return: lag x n x n. the i * lagth row is no longer the same as the response
#     """
#
#     n = coefs.shape[1]
#
#     coef_list = np.split(coefs, lag)
#
#     out_coef = np.zeros((lag, n, n))
#
#     for i in range(len(coef_list)):
#         coef = coef_list[i]
#
#
#         new_coefs = np.zeros(coef.shape)
#
#         for j in range(coef.shape[0]):
#             # For the jth row, the autocorrelation goes into the jth index
#             insert_coef = coef[0, j]
#             other_coef = coef[1:, j]
#
#
#
#             new_coef = np.insert(other_coef, j, insert_coef)
#
#             if verbose:
#                 print "Where to insert: ", j
#                 print "Insert coef: ", insert_coef
#                 print "Old coef: ", other_coef
#                 print "New coef: ", new_coef
#
#             new_coefs[:, j] = new_coef
#
#
#         out_coef[i, :, :] = new_coefs
#
#     return out_coef

def align_coefs(coefs, lag, verbose=False):
    """
    coefs: n* lag x n. the i * lagth row of col k is the auto-coefficient of k on k


    # align the coefs:
    # split into same lag
    return: lag x n x n. the
    """
    assert len(coefs.shape) == 2
    assert coefs.shape[0] / lag == coefs.shape[1]

    n = coefs.shape[1]

    coef_list = np.split(coefs, lag)

    out_coef = np.zeros((lag, n, n))

    for i in range(len(coef_list)):
        coef = coef_list[i]

        out_coef[i, :, :] = coef

    return out_coef

def unalign_coefs(coefs, lag):
    """
    The inverse of align_coefs

    coefs: lag x n x n
    # merge the lags to get the original operation

    return: n* lag x n. the i * lagth row of col k is the auto-coefficient of k on k

    """

    assert len(coefs.shape) == 3
    assert coefs.shape[0] == lag
    assert coefs.shape[1] == coefs.shape[2]

    n = coefs.shape[1]


    out_coef = np.zeros((n * lag, n))

    for i in range(lag):
        out_coef[i*n:(i+1)*n, :] = coefs[i, :, :]

    return out_coef


def remove_alphas(acoefs, lag):
    """
    :param acoefs: aligned coefs, lag x n x n, where for each lag k, ith row and jth column is coef of beta_k of
    cause gene i on effect gene j
    :param lag: the lag
    :return: acoefs_noalph: the same except diagonal entries are 0 (i.e. coef of gene j on gene j is 0 for all lags)
    """

    assert lag == acoefs.shape[0]
    assert acoefs.shape[1] == acoefs.shape[2]


    acoefs_noalph = acoefs.copy()

    for i in range(lag):
        np.fill_diagonal(acoefs_noalph[i, :, :], 0)

    return acoefs_noalph