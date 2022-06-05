import itertools
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns


def bin_proteome(df_curr, arg, bins):
    """For a single proteome: Performs binning and returns bin values."""

    x = np.asarray(df_curr[arg], float)
    if x.ndim > 1:
        x = x.squeeze()

    bin_kws = dict(bins=bins)
    hist_kws = {"density": True}

    hist, bin_edges = np.histogram(
        x,
        **bin_kws,
        **hist_kws,
    )
    hist = hist.astype(float) / hist.sum()

    return hist


def get_cis(df, p, bins, arg):
    """For a single proteome p, performs SE/CI calculation via bootstrapping for bins."""

    m = 1000
    n = len(bins) - 1

    df_curr = df[df["proteome"] == p]
    proteome_size = len(df_curr)
    artificials = np.ndarray(shape=(n, m))

    # Create m(=1000)x artificial/sampled proteome data
    for j in range(0, m):
        # Collect bin values for each (artificial) proteome of p
        extra_df = df_curr.sample(n=proteome_size, replace=True)
        bin_values = bin_proteome(extra_df, arg, bins)
        artificials[:, j] = bin_values

    # Calculate CIs for each bin
    cis = np.ndarray(n)
    for i in range(0, n):
        # SE bootstrapping formula: SE=SD(bin) over all artificial proteomes
        # CI formula: mean(/sum?) +- t_0.025 * SE
        # For t, use Student's t-distriution for 1000-1 = 999 = ~inf DOF and conf level = 0.95 -> 1.960
        t = 1.960
        cis[i] = artificials[i].std() * t

    return cis


def ci_per_bin(df, arg, bins):
    """Performs SE/CI calculation via bootstrapping for bins. Returns dictionary of CIs per proteome."""

    proteomes = pd.unique(df["proteome"])
    all_cis = {p: get_cis(df, p, bins, arg) for p in proteomes}

    return all_cis


def val_per_bin(df, arg, bins):
    """Performs binning of original data to get positions for error bars."""
    all_ys = {}

    # Get positions for error bars
    proteomes = pd.unique(df["proteome"])
    for p in proteomes:
        df_curr = df[df["proteome"] == p]
        all_ys[p] = bin_proteome(df_curr, arg, bins)

    # # Fill up for proteomes with less bins
    # max_num_regions = len(max(all_ys.values(), key=lambda v: len(v)))
    # for p in proteomes:
    #     all_ys[p] = np.resize(all_ys[p], max_num_regions)

    return all_ys


def plot_errorbars(all_cis, bins, all_y, palette, ax):
    x = bins[: (len(bins) - 1)]
    for p in all_cis.keys():
        ax.errorbar(
            x=x,
            y=all_y[p],
            yerr=all_cis[p],
            ecolor=(0, 0, 0, 0.4),
            elinewidth=1.5,
            fmt="o",
            color=palette[p],
            ms=1,
            zorder=10,
        )


def kl_via_binning(df, arg, bins):
    """Calculates the KL divergence between the resulting distributions."""

    proteomes = pd.unique(df["proteome"])
    bin_values = {}
    for p in proteomes:
        # Get bin values and add pseudo counts
        df_curr = df[df["proteome"] == p]
        distribution = bin_proteome(df_curr, arg, bins)
        bin_values[p] = distribution + math.exp(-12)

    # Calculate kl for every pair
    pairs = list(itertools.permutations(proteomes, 2))
    df_kl = pd.DataFrame({
        "first": [pair[0] for pair in pairs],
        "second": [pair[1] for pair in pairs],
        "value": [round(stats.entropy(bin_values[pair[0]], bin_values[pair[1]]), 2) for pair in pairs],
    })

    return df_kl


def plot_kl(df_kl, name_mapping, ax):
    """Plots a pairwise metric (KL) of multiple proteomes as specified in the given matrix in form of a heatmap."""

    df_matrix = df_kl.pivot(index="first", columns="second", values="value").fillna(0)
    names = list(name_mapping[p] for p in df_matrix.columns)
    df_matrix.index.name = ""
    df_matrix.columns.name = ""

    sns.heatmap(
        df_matrix,
        square=True,
        xticklabels=names,
        yticklabels=names,
        cbar_kws={
            "shrink": 0.5,
            "label": "KL",
        },
        annot=True,
        annot_kws={"fontsize": 9},
        ax=ax,
        cmap=plt.cm.binary,
    )
    ax.set_title("KL Between Whole Distributions", y=1.1)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=9)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=9)


