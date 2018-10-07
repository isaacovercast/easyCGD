import glob
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from collections import Counter
from itertools import combinations
from scipy.stats import entropy,linregress
from scipy.spatial.distance import pdist, squareform


#######################################################
## Stats functions
#######################################################


## Get one hill humber from a list of abundances (a column vector from the OTU table)
def hill_number(abunds, order):
    if order == 0:
        return len(np.nonzero(abunds)[0])
    if order == 1:
        h1 = np.exp(entropy(abunds))
        return h1
    tot = float(np.sum(abunds))
    proportions = np.array(abunds[abunds > 0])/tot
    prop_order = proportions**order
    h2 = np.sum(prop_order)**(1/(1-order))
    return h2


## Get all hill numbers from 0 to 'orders' from a column vector of abundances
def hill_numbers(abunds, orders, granularity=None, do_negative=False):
    ret = []
    min_order = 0
    if not granularity: granularity = orders + 1
    if do_negative:
        min_order = -orders
        granularity *= 2
    for order in np.linspace(min_order, orders, granularity):
        ret.append(hill_number(abunds, order))
    return np.array(ret)


def pi(fasta_file, quiet=True):
    ## Calculate pi from a fasta file
    pi = 0
    len_seq = 0
    try:
        f = open(fasta_file).readlines()
        ## Get just the sequences
        dat = [list(x.strip()) for x in f if ">" not in x]
        len_seq = len(dat[0])

        ## Transpose, so now we have a list of lists of all bases at each
        ## position.
        dat = np.transpose(np.array(dat))

        ## for each position
        for d in dat:
            ## If the position is _not_ monomorphic
            if len(Counter(d)) > 1:
                ## Enumerate the possible comparisons and for each
                ## comparison calculate the number of pairwise differences,
                ## summing over all sites in the sequence.
                base_count = Counter(d)
                ## ignore indels
                del base_count["-"]
                del base_count["N"]
                for c in combinations(base_count.values(), 2):
                    #print(c)
                    n = c[0] + c[1]
                    n_comparisons = float(n) * (n - 1) / 2
                    pi += float(c[0]) * (n-c[0]) / n_comparisons
        pi = pi/len_seq
    except ZeroDivisionError as inst:
        if not quiet: print("No sequences in file {}".format(fasta_file))
        pass
    except Exception as inst:
        print("Error in file - {}\n{}".format(fasta_file, inst))
        pi = 0
    ## Average over the length of the whole sequence.
    return pi


#######################################################
## Plotting functions
#######################################################


def plot_hill_numbers(in_df, order=5, do_negative=False, title="Hill numbers", quiet=True, axis_labels=True, ax=None, color=None):
    hill_results = {}
    for name, vals in in_df.items():
        hill_results[name] = hill_numbers(vals, orders=order, granularity=None, do_negative=do_negative)
    if not quiet: print(hill_results)
    if not ax:
        fig, ax = plt.subplots(figsize=(10, 10))
    if do_negative:
        X = np.linspace(-order, order, num=order*2+2)
        ax.axvline(0, color="k", alpha=0.2, linewidth=3, linestyle="--")

    else:
        X = range(order+1)
    for region, h in hill_results.items():
        ## If not specifying the color then we want to just use random different colors for each line
        if color:
            ax.semilogy(X, h, label=region, linewidth=4, alpha=0.65, color=color)
        else:
            ax.semilogy(X, h, label=region, linewidth=4, alpha=0.65)

    if axis_labels:
        plt.xlabel("Hill Number", fontsize=25)
        plt.ylabel("# Effective Species", fontsize=25)

    plt.xticks(range(order))
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.title(title, fontsize=30)
    plt.legend(fontsize=15)
    return ax


def plot_RACs(in_df, trim_zeros=True, title="Rank abundance", label=None, axis_labels=True, ax=None, color=None):
    if not ax:
        fig, ax = plt.subplots(figsize=(10, 10))
    if label == None:
       pass 
    for c in in_df.columns:
        all_abunds = in_df[c].fillna(0).values
        if trim_zeros:
            all_abunds = np.trim_zeros(np.sort(all_abunds)[::-1])
        else:
            all_abunds = np.sort(all_abunds)[::-1]
        X = np.arange(0, len(all_abunds))
        if color:
            ax.semilogy(X, all_abunds, label=c, linewidth=4, alpha=0.65, color=color)
        else:
            ax.semilogy(X, all_abunds, label=c, linewidth=4, alpha=0.65)
    if axis_labels:
        plt.xlabel("Rank", fontsize=25)
        plt.ylabel("log(Abundance)", fontsize=25)
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.title(title, fontsize=30)
    plt.legend(fontsize=15)
    return ax


def plot_distances(abunds, labels, metric="braycurtis", title="distance"):
    def mksq(abunds, metric):
        dist = pdist(abunds, metric=metric)
        return squareform(dist)    
    if isinstance(metric, list):
        sq1 = mksq(abunds, metric[0])
        sq2 = mksq(abunds, metric[1])

        idx = np.triu_indices(len(labels))
        sq1[idx] = 0
        sq2[idx] = 0
        sq1 += sq2.T
        ## Reformat as a string for the plot title
        metric = "{}/{}".format(metric[0], metric[1])
    else:
        sq1 = mksq(abunds, metric)
    
    fix, ax = plt.subplots(figsize=(8,7))
    heatmap = ax.pcolor(sq1, cmap="magma")
    cbar = plt.colorbar(heatmap)
    # Set ticks in center of cells
    _ = ax.set_xticks(np.arange(len(labels)) + 0.5, minor=False)
    _ = ax.set_yticks(np.arange(len(labels)) + 0.5, minor=False)

    _ = ax.set_xticklabels(labels, rotation=45, ha="right", rotation_mode="anchor")
    _ = ax.set_yticklabels(labels)
    ax.set_title("{} - {} Distance(s)".format(title, metric))
    return sq1


## Plot correlation between abundances and genetic diversities. Input should be dataframes
## with index indicating OTU for both frames.
def plot_abundance_diversity_correlation(a_df, g_df, ax=None, drop_zeros=False, color=None, log_transform=False,\
                                        do_proportions=False, title="Abundance Genetic Diversity Correlation",\
                                        quiet=True):
    if ax == None:
        _, ax = plt.subplots(figsize=(10, 10))
    if isinstance(color, str):
        pass
    elif np.any(color):
        pass
    else:
        color = np.random.rand(3,)

    ## if doing proportions scale all values to sum to 1
    if do_proportions:
        tmp_pis_df = g_df/g_df.sum(axis=0)
        tmp_abund_df = a_df/a_df.sum(axis=0)
    else:
        tmp_pis_df = g_df
        tmp_abund_df = a_df

    if log_transform:
        tmp_pis_df = g_df
        tmp_abund_df = np.log(a_df)

    ## Get global max x and y values for each ABC across all regions
    xmax = tmp_abund_df.max().values[0]
    ymax = tmp_pis_df.max().values[0]
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)

    ## Merge abundances and genetic diversities preserving OTU index, and drop any with missing data
    tmp_df = pd.concat([tmp_abund_df, tmp_pis_df], axis=1, sort=True).dropna()
    if len(tmp_df) == 0:
        print("abundance and genetic diversity inputs do not have the same indices, so this function is unavailable.")
        return
    tmp_df.columns = ["abund", "genet"]

    ## drop all rows with no genetic information
    if drop_zeros:
        tmp_df = tmp_df[tmp_df["genet"] > 0]

    ## Draw the linear regression line
    linreg = linregress(tmp_df["abund"], tmp_df["genet"])
    X = np.linspace(0, max([xmax, ymax]), 100)
    Y = X*linreg.slope + linreg.intercept
    ax.plot(X, Y, color=color, alpha=.5)

    ax.scatter(tmp_df["abund"], tmp_df["genet"], label="{:.4f} - {}".format(linreg.rvalue**2, "MESS Simulations"), color=color)

    if not quiet: print(linreg)
    _ = plt.xlabel("Abundance", fontsize=15)
    _ = plt.ylabel("Genetic Diversity", fontsize=15)
    _ = plt.legend(loc="upper right")
    _ = plt.suptitle(title, fontsize=20)
    return ax

#######################################################
## Loading functions
#######################################################

def get_pis_from_fastas(fasta_dir, outfile=None, colname=None, quiet=True):
    try:
        files = glob.glob(os.path.join(fasta_dir, "*.fasta"))
        if not colname:
            colname = [fasta_dir.strip("/").split("/")[-1]]
        if not isinstance(colname, list):
            colname = [colname]
        pis = {}
        for f in files:
            ## Just use the OTU/species name for the pis file
            pis[f.split("/")[-1].split(".")[0]] = pi(f, quiet)
            pis_df = pd.DataFrame.from_dict(pis, orient="index", columns=colname)
        if outfile:
            with open(outfile, 'w') as out:
                pis_df.to_csv(outfile)
    except Exception as inst:
        print("Error in getting pis {}".format(inst))
        raise
    return pis_df


## Read in abundances and/or pi data. Accepts either a flat file with data as a row vector
## or a dataframe oriented as a column
def load_data(pis_file, labels=None):
    ## Allow to read in both one line style pi vectors as well as pd.dataframe style
    ## which is column oriented.
    pis_df = pd.DataFrame()
    with open(pis_file) as infile:
        dat = infile.readlines()
    if len(dat) > 1:
        pis_df = pd.concat([pis_df, pd.read_csv(pis_file, index_col=0)], ignore_index=True, axis=1)
    else:
        pi_df = pd.read_csv(pis_file, header=None).T
        pis_df = pd.concat([pis_df, pi_df], ignore_index=True, axis=1)
    ## relabel columns if labels are passed in
    if labels:
        if not isinstance(labels, list):
            labels = [labels]
    else:
        ## If no labels then just "guess" based on the file name
        labels = [pis_file.split("/")[-1].split(".")[0]]
    pis_df.columns = labels 
    return pis_df


## Load multiple files at once into a concatenated dataframe. Convenience function
## since this is something that happens regularly and the housekeeping is annoying.
def load_dfs(df_files, labels=None):
    dfs = {x.split("/")[-1]:load_data(x) for x in df_files}
    df = pd.concat(dfs.values(), ignore_index=False, axis=1, sort=True)
    return df


