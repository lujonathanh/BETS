import numpy as np
import pandas as pd
import sys
import os
import collections
import matplotlib.pyplot as plt
import fdr_control as fc
import time
import pickle

fdr = 0.2

orig_run_folder = "run_l-fdr"

# MODIFY THE BELOW, REPLACING URANDFOLDERNAME WITH THE NAME OF THE FOLDER WITH THE PERMUTED RUN
# e.g. replace URANDFOLDERNAME with insilico_size100_1_urand-0mean-reps-enet-2-g
urand_run_folder = "../URANDFOLDERNAME/run_l-fdr"




################################################################################


outdir = "analysis"
signetworkname = os.path.join(outdir, "bootstrap-fdr-" + str(fdr) + "-network" + ".txt")


# get the exact run_l-d
orig_bootstrap_filename = [os.path.join(orig_run_folder, "bootstrap-results-fdr-0.05-effect", f) for f in os.listdir(os.path.join(orig_run_folder, "bootstrap-results-fdr-0.05-effect")) if f.endswith("-union-bootstrap-network.txt")][-1]
urand_bootstrap_filename = [os.path.join(urand_run_folder, "bootstrap-results-fdr-0.05-effect", f) for f in os.listdir(os.path.join(urand_run_folder, "bootstrap-results-fdr-0.05-effect")) if f.endswith("-union-bootstrap-network.txt")][-1]


if not os.path.exists(outdir):
    os.makedirs(outdir)

# # WARNING: this assumes that the networks and corresponding files are already ordered and will correspond


orig_bnet = pd.read_csv(orig_bootstrap_filename, sep="\t")
urand_bnet = pd.read_csv(urand_bootstrap_filename, sep="\t")
        
print("Original")
print(orig_bootstrap_filename)
print(orig_bnet.shape)
print(orig_bnet.head())
print()

print("Umbrella null")
print(urand_bootstrap_filename)
print(urand_bnet.shape)
print(urand_bnet.head())

bootstrap_freq_key = "Bootstrap:Freq"

orig_freqs = orig_bnet[bootstrap_freq_key].values
urand_freqs = urand_bnet[bootstrap_freq_key].values


orig_naboves = []
urand_naboves = []
threshes = []

print(fdr)

thresh = fc.FDR_above_threshold(orig_freqs, urand_freqs, fdr)
print("Thresh: ", thresh)
threshes.append(thresh)

orig_naboves.append(len(np.where(orig_freqs >= thresh)[0]))
urand_naboves.append(len(np.where(urand_freqs >= thresh)[0]))

orig_thresh_bnet = orig_bnet[orig_bnet[bootstrap_freq_key] >= thresh]

orig_thresh_bnet.to_csv(signetworkname, sep="\t", index=0)
print("Written to ", signetworkname)





umbrella_result_df = pd.DataFrame()
umbrella_result_df["FDR"] = [fdr]
umbrella_result_df["Bootstrap Frequency Threshold"] = threshes
umbrella_result_df["Orig Bootstrap Network"] = orig_naboves
umbrella_result_df["Null Bootstrap Network"] = urand_naboves

umbrella_result_df


# In[21]:





out_filename = signetworkname[:-4] + "-umbrella-results.csv"
umbrella_result_df.to_csv(out_filename, index = False)
print("Written to ", out_filename)


# In[22]:

#
# linestyles = ["-", "--", "-.", ":"]
# linestyles = linestyles + linestyles
#
# plt.figure(figsize=(12,8))
# plt.title("Edge Frequency in Bootstrap networks", fontsize=18)
# plt.hist(orig_freqs, alpha=0.5, label="Original", color='r', bins=100)
# plt.hist(urand_freqs, alpha=0.5, label="Null", color='b', bins=100)
# for i, row in umbrella_result_df.iterrows():
#     if row["FDR"] < 1.0:
#         plt.axvline(row["Bootstrap Frequency Threshold"], label= "FDR = " + str(row["FDR"]),
#                    linestyle = linestyles[i], linewidth=5, color='g')
# plt.legend(loc="best")
# plt.tick_params(labelsize=14)
# plt.xlabel("Frequency among Bootstrap FDR-thresholded networks", fontsize=18)
# plt.ylabel("Count", fontsize=18)
# plt.show()
# plt.close()
#
#
# plt.figure(figsize=(12,8))
# plt.ylim((0,100))
# plt.title("Edge Frequency in Bootstrap networks", fontsize=18)
# plt.hist(orig_freqs, alpha=0.5, label="Original", color='r', bins=100)
# plt.hist(urand_freqs, alpha=0.5, label="Null", color='b', bins=100)
# for i, row in umbrella_result_df.iterrows():
#     if row["FDR"] < 1.0:
#         plt.axvline(row["Bootstrap Frequency Threshold"], label= "FDR = " + str(row["FDR"]),
#                    linestyle = linestyles[i], linewidth=5, color='g')
# plt.legend(loc="best")
# plt.tick_params(labelsize=14)
# plt.xlabel("Frequency among Bootstrap FDR-thresholded networks", fontsize=18)
# plt.ylabel("Count", fontsize=18)
# plt.show()
# plt.close()



print(fdr)

thresh = fc.FDR_above_threshold(orig_freqs, urand_freqs, fdr)
orig_thresh_bnet = orig_bnet[orig_bnet[bootstrap_freq_key] >= thresh]
print("Edges: ", orig_thresh_bnet.shape[0])


# In[24]:


print("Causal Genes", len(np.unique(orig_thresh_bnet["Cause"].values)))

print("Effect Genes", len(np.unique(orig_thresh_bnet["Effect"].values)))

print("Genes in network ", len(np.unique(np.concatenate((orig_thresh_bnet["Cause"].values,
                                             orig_thresh_bnet["Effect"].values)))))

print("# Edges ", orig_thresh_bnet.shape[0])


# # Plot the Thresholds

# In[34]:


# fdr = 0.2

# bins = np.linspace(0, 1.0, 40)
# fig = plt.figure(figsize=(12,8))
# plt.ylim((0,500))
# #plt.title("Edge Frequency in Bootstrap networks", fontsize=18)
# plt.hist(orig_freqs, alpha=0.5, label="Original", color='r', bins=bins)
# plt.hist(urand_freqs, alpha=0.5, label="Null", color='b', bins=bins)

# plt.axvline(umbrella_result_df[umbrella_result_df["FDR"] == fdr]["Bootstrap Frequency Threshold"].values[0], label= "FDR = " + str(fdr),
#                    linestyle = "dashed", linewidth=8, color='k')
# plt.legend(loc="best", fontsize=20)
# plt.tick_params(labelsize=20)
# plt.xlabel("Edge Bootstrap Frequency", fontsize=20)
# plt.ylabel("Count", fontsize=20)


# outfile = "analysis/plots/joint-unperturbed-enet-2_bootstrap-fdr-0.2/bootstrap-frequency-distribution.pdf"
# fig.savefig(outfile)
# print("Written to ", outfile)

# plt.show()
# plt.close()


# In[ ]:

