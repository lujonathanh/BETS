# BETS
Causal Network Inference and Analysis from Time Series

BETS is a Python package that infers causal networks from time series data. 

# Is BETS right for you?

BETS is a good method for your problem if:
* you have a dataset of gene expression values (or other assay type) over time
  * you can have multiple replicates 
  * the time points are *about* equally spaced

The method requires that all input  TODO

# Requirements

BETS requires installation of the following libraries:

* numpy (1.13.1 or above)
* scipy (0.19.1 or above)
* pandas (0.20.3 or above)
* matplotlib (1.4.3 or above)
* sklearn (0.19.0 or above)

It has been tested in Python 2.7 on MacOSX v.10.11.6.

# How To Run

See `BETS_tutorial.md` for a step-by-step walk through of BETS.

# Questions?

Please post them at [our google group](https://groups.google.com/forum/#!forum/bets-support)

# About the Method 

BETS is short for "Bootstrap Elastic net regression from Time Series", a statistical framework based on Granger causality for the recovery of a directed gene network from transcriptional time series data. applies regularized vector autoregression along with a permutation-based 
null and False Discovery control to infer causal networks. It was designed
for a high-dimensional gene-expression time series data. 

# Authors

BETS was developed by Jonathan Lu, Bianca Dumitrascu, and Professor Barbara Engelhardt in the [Engelhardt Group](beehive.cs.princeton.edu) at the Department of Computer Science at Princeton University over 2016-2019.

# Citation
This work is in submission.