__author__ = 'jlu96'

import pandas as pd
import numpy as np
import scipy.stats as stats
import pickle
import sys
from sklearn.linear_model import LinearRegression
import geneTSmunging as gtm


def fit_all_pairwise_conditional(geneTS, lag, rows, coeflag_options, has_reps=True,
                                                               only_array=False):
    """
    rows: the effect genes to use
    """

    assert isinstance(lag, int)

    if rows == None:
        rows = range(geneTS.shape[0])


    if coeflag_options == None:
        coeflag_options = range(1, lag + 1)

    assert hasattr(rows, '__iter__')
    assert hasattr(coeflag_options, '__iter__')

    if rows == None:
        rows = range(geneTS.shape[0])


    if coeflag_options == None:
        coeflag_options = range(1, lag + 1)


    n = geneTS.shape[0]
    T = geneTS.shape[1]


    coefs = np.zeros((n * lag, n))

#     print "Gene TS shape:" , geneTS.shape


    # you just ignore the ones that are self-on-self
    for effect_index in rows:

        print "Effect # ", effect_index

        if has_reps:
            effect_ts = geneTS[effect_index, :, :]
        else:
            effect_ts = geneTS[effect_index, :].flatten()

        for cause_index in range(n):

            if cause_index != effect_index:
                if has_reps:
                    cause_ts = geneTS[cause_index, :, :]
                else:
                    cause_ts = geneTS[cause_index, :].flatten()


                for coeflag in coeflag_options:

                    slope, _, _, _, _, _, _, _, _, _ = fit_pairwise_conditional(cause_ts=cause_ts,
                                                                   effect_ts=effect_ts,
                                                                   lag=lag,
                                                                   coeflag=coeflag,
                                                                   T=T,
                                                                   has_reps=has_reps)

                    coefs[(cause_index + (coeflag - 1) * n), effect_index] = slope

    return coefs







def fit_pairwise_conditional(cause_ts, effect_ts, lag, coeflag, T, has_reps):
    """
    Simplest form: cause_ts and effect_ts are literally just arrays.

    cause_ts: T (x r)
    effect_ts: T (x r)


    For edge C-> D, I want the coefficient c from:

    D_t = c C_{t-coeflag} + sum_{i=1}^lag d_i D_{t-i}  + intercept

    Then we fit ResidualCD_t = c C_{t-1}

    """

    assert cause_ts.shape[0] == effect_ts.shape[0]
    assert cause_ts.shape[0] == T
    assert isinstance(lag, int) and (lag > 0)



    if has_reps:
        # T x r

        cause_pts = cause_ts[(lag - coeflag): (T - coeflag), :]

        effect_pts = effect_ts[lag:T, :]

        prev_effect_pts = np.array([effect_ts[(lag - x): (T - x), :] for x in range(1, lag + 1)])

    else:

        cause_pts = cause_ts[(lag - coeflag): (T - coeflag)]

        effect_pts = effect_ts[lag:T]

        # Take the lagged versions starting from lag 1

        prev_effect_pts = np.array([effect_ts[(lag - x): (T - x)] for x in range(1, lag + 1)])


    # Y is T-lag x r

    Y = effect_pts.T
    n = Y.shape[0]

    # flatten the replicates as independent samples. Y was of form r x T
    if has_reps:
        Y = np.reshape(Y, (Y.shape[1] * Y.shape[0],))


#     print "Cause TS shape is ", cause_ts.shape
#     print "Cause PTS shape is ", cause_pts.shape
#     print "Effect TS shape is ", effect_ts.shape
#     print "Effect PTS shape is ", effect_pts.shape
#     print "Prev Effect PTS shape is ", prev_effect_pts.shape


#     print "Y shape is ", Y.shape

    X = np.append(np.array([cause_pts]), prev_effect_pts, axis=0).T

    # flatten the replicates as independent samples. X was of form r x T x n
    if has_reps:
        X = np.reshape(X, (X.shape[1] * X.shape[0], X.shape[2]))

    p_all = X.shape[1]

#     print "X shape: ", X.shape



    X_nocause = prev_effect_pts.T

    # flatten the replicates as independent samples. X was of form r x T x n
    if has_reps:
        X_nocause = np.reshape(X_nocause, (X_nocause.shape[1] * X_nocause.shape[0], X_nocause.shape[2]))

    p_nocause = X_nocause.shape[1]

#     print "X_nocause shape: ", X_nocause.shape

    lmodel = LinearRegression()

    lmodel.fit(X, Y)
    Y_pred = lmodel.predict(X)
    RSS_all = (np.power(Y_pred - Y, 2)).sum()

    lmodel_nocause = LinearRegression()
    lmodel_nocause.fit(X_nocause, Y)
    Y_pred_nocause = lmodel_nocause.predict(X_nocause)
    RSS_nocause = (np.power(Y_pred_nocause - Y, 2)).sum()


    p_value, F = F_test(RSS_1=RSS_nocause, RSS_2=RSS_all, n=n,
                       p_1=p_nocause, p_2=p_all)


    ## F = [(RSS_1 - RSS_2)/(p_2 - p_1)]/[RSS_2/(n - p_2)]



    #print p_value

    slope = lmodel.coef_[0]


    intercept = lmodel.intercept_
    r2 = lmodel.score(X, Y)



    #print "Coef: ", lmodel.coef_

    #slope_pv =
    # get the p-value for the slope (don't use the F-test??)

    #Y_pred, fit_result = fm.compute_fit(X, Y, lmodel.coef_, intercept)

    return slope, intercept, p_value, F, r2, X, Y, cause_pts, effect_pts, prev_effect_pts


def F_test(RSS_1, RSS_2, n, p_1, p_2, silent_RSS_error=False):
    """
    RSS_1: Residual Sum of Squares from smaller model
    RSS_2: RSS from larger model (more parameters)
    n: sample size
    p_1: smaller number of parameters
    p_2: larger number of parameters
    """
    assert p_2 > p_1
    if RSS_1 < RSS_2:
        if silent_RSS_error:
            print "Smaller model RSS1 = ", RSS_1
            print "Smaller model p = ", p_1
            print "Larger model RSS2 = ", RSS_2
            print "Larger model p = ", p_2
            print "Larger model has higher RSS??"
            return None, None
        else:
            raise ValueError("Smaller model RSS " + str(RSS_1) + " less than larger model RSS" + str(RSS_2))

    F = ((RSS_1 - RSS_2)/(p_2 - p_1))/(RSS_2/(n - p_2))
    p_value = 1.0 - stats.f.cdf(F, p_2 - p_1, n - p_2)

    return p_value, F


def get_parser():
    # Parse arguments
    import argparse

    description = 'Perform all the fits.'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-d', '--data_file', required=True)

    parser.add_argument('-lr', '--load_reps', required=True, type=int)

    parser.add_argument('-o', '--out_prefix', required=True, help="Prefix for saving results")

    parser.add_argument('-l', '--lag',  type=int, required=True)

    parser.add_argument('-rl', '--row_file', default=None)

    return parser




def load_and_run(args):

    data_file = args.data_file

    lag = args.lag

    if args.row_file != None:
        rows = pickle.load(open(args.row_file, 'rB'))
    else:
        rows = None


    # Load data file    # Load data file
    if args.load_reps:
        genes, geneTS = gtm.load_basic_rep_file_list(data_file)

    else:
        df = pd.read_csv(data_file, sep="\t")
        genes, geneTS = gtm.get_gene_TS(df)


    coefs = fit_all_pairwise_conditional(geneTS=geneTS, lag=lag, rows=rows, coeflag_options=None, has_reps=args.load_reps)


    outfile = args.out_prefix + "_coefs.p"
    pickle.dump(coefs, open(outfile, 'w'))
    print "Coefs saved to ", outfile


def main():
    load_and_run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()