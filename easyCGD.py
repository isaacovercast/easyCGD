import glob
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from collections import Counter
from itertools import combinations
from scipy.stats import entropy


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


def plot_hill_numbers(combined_df, order=5, granularity=None, do_negative=False, title="Hill numbers for ABC combined", quiet=True):
    hill_results = {}
    for region in regions:
        all_abunds = np.concatenate([combined_df[region][ABC].values for ABC in ABCs])
        hill_results[region] = hill_numbers(np.array(all_abunds), orders=order, granularity=granularity, do_negative=do_negative)
    if not quiet: print(hill_results)
    fig, ax = plt.subplots(figsize=(10, 10))
    if do_negative:
        X = np.linspace(-order, order, num=order*2+2)
        ax.axvline(0, color="k", alpha=0.2, linewidth=3, linestyle="--")

    else:
        X = range(order+1)
    for region, h in hill_results.items():
        ax.semilogy(X, h, label=region, linewidth=4, alpha=0.65)

    plt.title(title)
    plt.legend()
    return hill_results


def plot_hill_numbers(in_df, order=5, do_negative=False, title="Hill numbers", quiet=True, ax=None, color=None):
    hill_results = {}
    all_abunds = in_df[in_df.columns[0]].values
    hill_results[in_df.columns.values[0]] = hill_numbers(all_abunds, orders=order, granularity=None, do_negative=do_negative)
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

    plt.xticks(range(order))
    plt.title(title, fontsize=20)
    plt.legend(fontsize=15)
    return ax


def plot_RACs(in_df, title="Rank abundance", label=None, ax=None, color=None):
    all_abunds = in_df[in_df.columns[0]].values
    all_abunds = np.trim_zeros(np.sort(all_abunds)[::-1])
    #all_abunds = np.sort(all_abunds)[::-1]
    X = np.arange(0, len(all_abunds))
    if not ax:
        fig, ax = plt.subplots(figsize=(10, 10))
    if label == None:
       pass 
    if color:
        ax.semilogy(X, all_abunds, label=in_df.columns.values[0], linewidth=4, alpha=0.65, color=color)
    else:
        ax.semilogy(X, all_abunds, label=in_df.columns.values[0], linewidth=4, alpha=0.65)
    plt.title(title, fontsize=20)
    plt.legend(fontsize=15)
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
