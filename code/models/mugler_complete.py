model = """
basal_comK:
    $pool > mRNAk
    k_1

bind_comK:
    $pool > mRNAk
    k_2 * comK**h / (k_k**h + comK**h)

trans_comK:
    mRNAk > mRNAk + comK
    k_3

basal_comS:
    $pool > mRNAs
    k_4 bind_comS: $pool > mRNAs k_5/ (1 + (comK/k_s)**p)
    trans_comS: mRNAs > mRNAs + comS
    k_6

deg_mRNAk:
    mRNAk > $pool
    k_7

deg_comK:
    comK > $pool
    k_8

deg_mRNAs:
    mRNAs > $pool
    k_9

deg_comS:
    comS > $pool
    k_10

bind_mecAk:
    mecA + comK > mecAk
    k_11/omega

unbind_mecAk:
    mecAk > mecA + comK
    k_11un

deg_mecAk:
    mecAk > mecA
    k_12

bind_mecAs:
    mecA + comS > mecAs
    k_13/omega

unbind_mecAs:
    mecAs > mecA + comS
    k_13un

deg_mecAs:
    mecAs > mecA
    k_14


#initialization
mRNAk = 0
mRNAs = 0
comK = 0
comS = 0
mecA = 500
mecAs = 0
mecAk = 0

#parameters
k_1 = 0.025
k_2 = 0.1 9
k_3 = 0.2
k_4 = 0
k_5 = 5000
k_6 = 0.2
k_7 = 0.005
k_8 = 0.0001
k_9 = 0.005
k_10 = 0.0001
k_11 = 0.000002
k_11un = 0.0005
k_12 = 0.05
k_13 = 0.0000045
k_13un = 0.00005
K_14 = 0.00004
omega = 1.66
h = 2
p = 5
k_k = 5000
k_s = 8333"""
