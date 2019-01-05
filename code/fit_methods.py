__author__ = 'jlu96'


import numpy as np
from sklearn.linear_model import Lasso, Ridge, ElasticNet, LinearRegression
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.model_selection import LeaveOneOut
import collections
from lag_conversion import get_XY_lagged
import pandas as pd
import time


def fit_ols(X, Y, verbose=False,**kwargs):
    """
    X: n x p matrix
    Y: n x 1 matrix

    hyperparams: regularization used for lasso
    return:
    coef: p x 1
    """

    assert X.shape[0] == Y.shape[0]

    n = X.shape[0]
    p = X.shape[1]

    linreg = LinearRegression()

    linreg.fit(X, Y)

    coef = np.reshape(linreg.coef_, (p, 1))
    intercept = linreg.intercept_



    Y_pred, fit_result = compute_fit(X, Y, coef, intercept)

    if verbose:
        print "Diff in prediction"
        print linreg.predict(X).shape
        print Y_pred.shape
        print Y_pred - np.reshape(linreg.predict(X), (n,1))

    return Y_pred, coef, intercept, fit_result

def fit_lasso(X, Y, verbose=False, max_iter = 50000, **kwargs):
    """
    X: n x p matrix
    Y: n x 1 matrix

    hyper: regularization used for lasso
    return:
    coef: p x 1
    """

    assert X.shape[0] == Y.shape[0]

    n = X.shape[0]
    p = X.shape[1]

    alpha = kwargs["hyper"]

    lasso = Lasso(alpha=alpha, max_iter=max_iter, selection='random')

    lasso.fit(X, Y)

    coef = np.reshape(lasso.coef_, (p, 1))
    intercept = lasso.intercept_



    Y_pred, fit_result = compute_fit(X, Y, coef, intercept)

    fit_result["dual_gap"] = lasso.dual_gap_

    # Temp adding 2/6/17 -JLu
    if verbose:
        print "Tolerance: "
        print lasso.tol

        print "Duality gap: "
        print lasso.dual_gap_
        print coef

        # print "Diff in prediction"
        # print lasso.predict(X).shape
        # print Y_pred.shape
        # print Y_pred - np.reshape(lasso.predict(X), (n,1))

    return Y_pred, coef, intercept, fit_result

def fit_ridge(X, Y, verbose=False, max_iter=50000, **kwargs):
    """
    X: n x p matrix
    Y: n x 1 matrix

    hyper: regularization used for ridge
    return:
    coef: p x 1
    """

    assert X.shape[0] == Y.shape[0]

    n = X.shape[0]
    p = X.shape[1]

    alpha = kwargs["hyper"]

    ridge = Ridge(alpha=alpha, max_iter=max_iter)

    ridge.fit(X, Y)

    coef = np.reshape(ridge.coef_, (p, 1))
    intercept = ridge.intercept_



    Y_pred, fit_result = compute_fit(X, Y, coef, intercept)

    if verbose:
        print "Diff in prediction"
        print ridge.predict(X).shape
        print Y_pred.shape
        print Y_pred - np.reshape(ridge.predict(X), (n,1))

    return Y_pred, coef, intercept, fit_result

def fit_enet(X, Y, verbose=False, max_iter=30000, **kwargs):
    """
    X: n x p matrix
    Y: n x 1 matrix

    hyper: (alpha, l1_ratio)

    Minimizes the objective function:
    1 / (2 * n_samples) * ||y - Xw||^2_2
    + alpha * l1_ratio * ||w||_1
    + 0.5 * alpha * (1 - l1_ratio) * ||w||^2_2
    return:
    coef: p x 1
    """

    assert X.shape[0] == Y.shape[0]

    n = X.shape[0]
    p = X.shape[1]

    alpha, l1_ratio = kwargs["hyper"]

    enet = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, max_iter=max_iter,
                      selection='random')

    enet.fit(X, Y)

    coef = np.reshape(enet.coef_, (p, 1))
    intercept = enet.intercept_



    Y_pred, fit_result = compute_fit(X, Y, coef, intercept)

    fit_result["dual_gap"] = enet.dual_gap_

    if verbose:
        print "Diff in prediction"
        print enet.predict(X).shape
        print Y_pred.shape
        print Y_pred - np.reshape(enet.predict(X), (n,1))

    return Y_pred, coef, intercept, fit_result


def durbin_watson(resid, lag=1, verbose=False):
    """
    Calculates the Durbin-Watson statistic:

    DW = sum([resid_{i} - resid_{i-lag}]^2)/sum(resid_{i}^2)

    Around 2 is no correlation-- residuals well captured by model.
    Close to 1 is positive correlation.
    """
    resid = resid.flatten()

    # sum of squared auto-residuals
    ssar = np.sum(np.diff(resid,n=lag)**2)
    ssr = np.sum(resid**2)

    if verbose:
        print np.diff(resid,n=lag), resid
        print ssar, ssr

    return ssar * 1.0 / ssr


def predict(X, coef, fit_intercept):
    assert X.shape[1] == coef.shape[0]

    return np.dot(X, coef) + fit_intercept

def compute_fit(X, Y, coef, fit_intercept, dw_lags=[1]):
    """
    X: n x p matrix
    Y: n x o matrix
    coef: p x o
    """
    assert X.shape[1] == coef.shape[0]
    assert X.shape[0] == Y.shape[0]
    assert Y.shape[1] == coef.shape[1]

    Y_pred = predict(X, coef, fit_intercept)


    resid = Y_pred - Y



    fit_result = collections.OrderedDict()
    fit_result["r2"] = r2_score(Y, Y_pred)
    fit_result["mse"] = mean_squared_error(Y, Y_pred)
    fit_result["sse"] = np.sum((resid)**2)
    fit_result["n"] = Y.shape[0]
    fit_result["df"] = len(np.nonzero(coef)[0])
    for dw_lag in dw_lags:
        fit_result["DW:Lag" + str(dw_lag)] = durbin_watson(resid, lag=dw_lag)


    return Y_pred, fit_result


def perform_test(X_matr, Y_matr, lag, fit_method, replace_row,
                 has_reps=False, bootstrap=False, seed=None, **kwargs):
    """
    X_matr: n x T (x r) matrix of input genes
    Y_matr: 1 x T matrix of output gene
    lag: lag
    fit_method: one of the fit_ methods, e.g. fit_lasso
    hyper: hyperparams for calling test
    replace_rows: which row Y_matr is in in X
    bootstrap: do it?
    seed: random seed for bootstrap

    return: X_t, Y_t, Y_pred, coef, intercept, fit_result
    """

    ## Get X and Y lagged

    X_t, Y_t = get_XY_lagged(X_matr, Y_matr, lag, replace_row=replace_row,
                             has_reps=has_reps, bootstrap=bootstrap, seed=seed)

    t = time.time()
    Y_pred, coef, intercept, fit_result = fit_method(X_t, Y_t, **kwargs)
    fit_result["time"] = round(time.time() - t, 3)

    return X_t, Y_t, Y_pred, coef, intercept, fit_result


def perform_test_random(X_matr, rand_X_matr, Y_matr, lag, fit_method, replace_row,
                 has_reps=False, verbose=False, bootstrap=False, seed=None, **kwargs):
    """
    Perform a single fit.

    X_matr: n x T (x r) matrix of input genes
    rand_X_matr: n x T (x r) matrix of randomized input genes, where randomized across time.
    Y_matr: 1 x T (x r) matrix of output gene
    lag: lag
    fit_method: one of the fit_ methods, e.g. fit_lasso
    hyper: hyperparams for calling test
    replace_rows: which row Y_matr is in in X. Set to None, otherwise. If it is inside, set the coefficients to zero.

    return: coef, intercept, fit_result. Note the 0, i * lag indices of the coef are for the output genes.
    """

    ## Get X and Y lagged

    # iterate through all possible predictors

    n = X_matr.shape[0]
    T = X_matr.shape[1]

    coef = np.zeros(n * lag)

    coef_temps = np.zeros( (n*lag, n))
    intercept_temps = np.zeros((1, n))

    # There will be m different fits. fit_result should just average all of the individual fit_results
    fit_result_temps = []

    X_t_orig, Y_t_orig = get_XY_lagged(X_matr, Y_matr, lag, replace_row=replace_row,
                                     has_reps=has_reps, bootstrap=bootstrap, seed=seed)
    Y_pred_orig, coef_orig, intercept_orig, fit_result_orig = fit_method(X_t_orig, Y_t_orig, **kwargs)



    for p in range(n):

        if p != replace_row:

            X_matr_temp = X_matr.copy()

            # replace the predictor row with the randomized row
            X_matr_temp[p] = rand_X_matr[p]



            # X_t_temp is a matrix of T - lag replace_rows where for each row j, i = j + lag-1
            #[A_i, B_i, .... Z_i, A_{i-1}, B_{i-1}, .... Z_{i-1}......  A_{i-lag+1}, B_{i-lag+1}, .... Z_{i-lag+1}

            X_t_temp, Y_t_temp = get_XY_lagged(X_matr_temp, Y_matr, lag, replace_row=replace_row,
                                     has_reps=has_reps, bootstrap=bootstrap, seed=seed)




            # coef_temp is (A_i, B_i, .... Z_i, A_{i-1}, B_{i-1}, .... Z_{i-1}......  A_{i-lag+1}, B_{i-lag+1}, .... Z_{i-lag+1}, 1)

            # we only want the p, p + n,... p + (lag - 1) * n indices.

            t = time.time()
            Y_pred_temp, coef_temp, intercept_temp, fit_result_temp = fit_method(X_t_temp, Y_t_temp, **kwargs)
            
            fit_result_temp["time"] = round(time.time() - t, 3)



            # Not
            coef[p: p + (lag - 1) * n + 1: n] = coef_temp[p : p + (lag - 1) * n + 1: n].flatten()

            # Store in the list of all coefs
            coef_temps[:, p] = coef_temp.flatten()

            intercept_temps[:, p] = intercept_temp


            fit_result_temps.append(fit_result_temp)

            if verbose:

                print "Original TS: ", X_matr[p]
                print "Randomized TS: ", X_matr_temp[p]

                print "Original X_t: ", X_t_orig[:T - lag, p:-1:n]
                print "Randomized X_t", X_t_temp[:T - lag, p:-1:n]

                print "Same Y_t?", (Y_t_orig == Y_t_temp).all()

                print "Orig coefs: ", coef_orig[p: lag * n + 1: n]
                print "Right around: "
                print coef_orig[p - 1: lag * n + 1: n]
                print coef_orig[p + 1: lag * n + 1: n]

                print "Updated coefs: ", coef_temp[p: lag * n + 1: n]
                print "Right around: "
                print coef_temp[p - 1: lag * n + 1: n]
                print coef_temp[p + 1: lag * n + 1: n]

                print "Updated coefs:", coef
                print "Y_pred_temp", Y_pred_temp

    coef = coef.reshape(n * lag, 1)

    fit_result_temps_df = pd.DataFrame(fit_result_temps)

    fit_result_std_dict = fit_result_temps_df.std()

    fit_result = fit_result_temps_df.mean().to_dict()

    keys = fit_result.keys()

    for key in keys:
        fit_result[key + "_mean"] = fit_result[key] # LEFT OFF HERE 1/25
        del fit_result[key]
        fit_result[key + "_std"] = fit_result_std_dict[key]

    return coef, fit_result, coef_temp, coef_temps, intercept_temps



def perform_loto_cv(X_matr, Y_matr, lag, fit_method, replace_row, verbose=False, has_reps=False,
                    **kwargs):
    """
    Perform leave-one-timepoint-out cross-validation.

    X_matr: n x T (x r) matrix of input genes, where r is # reps
    Y_matr: 1 x T (x r) matrix of output gene
    lag: lag
    fit_method: one of the fit_ methods, e.g. fit_lasso
    hyper: hyperparams for calling test
    replace_row: which row Y_matr is in in X

    return: fit_result
    """

    X_t, Y_t = get_XY_lagged(X_matr, Y_matr, lag, replace_row=replace_row, has_reps=has_reps)


    T_test = X_t.shape[0]

    loo = LeaveOneOut().split(X_t)


    Y_tests = np.zeros(T_test)
    Y_preds = np.zeros(T_test)
    dfs = np.zeros(T_test)

    for train_index, test_index in loo:
        X_train = X_t[train_index]
        Y_train = Y_t[train_index]
        X_test = X_t[test_index]
        Y_test = Y_t[test_index]


        _, coef, intercept, _ = fit_method(X_train, Y_train, **kwargs)

        Y_pred = predict(X_test, coef, intercept)

        Y_tests[test_index] = Y_test
        Y_preds[test_index] = Y_pred
        dfs[test_index] = len(np.nonzero(coef)[0])



    mse = mean_squared_error(Y_tests, Y_preds)
    sse = np.sum((Y_tests - Y_preds)**2)
    avg_df = np.average(dfs)
    r2 = r2_score(Y_tests, Y_preds)
    # formula for adjusted R^2 : R^2_a = (n-1) R^2 / (n- k - 1) - k/(n - k -1)
    # note: this is a little hacky, since the dfs vary in between
    #r2_adj =
    # TODO - Jonathan

    fit_result = collections.OrderedDict()
    fit_result["n"] = T_test
    fit_result["lag"] = lag
    fit_result["mse"] = mse
    fit_result["sse"] = sse
    fit_result["avg_df"] = avg_df
    fit_result["r2"] = r2
    #fit_result["r2_adj"] = r2_adj



    if verbose:
        print "Y_tests: ", Y_tests
        print "Y_preds: ", Y_preds
        print fit_result

    return fit_result