from abc import ABC, abstractmethod

import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import pandas as pd
import seaborn as sns

from ppprint.visualization.plot import Plot
from ppprint.visualization.plot_extras import (
    val_per_bin,
    ci_per_bin,
    plot_errorbars,
    kl_via_binning,
    plot_kl,
)


class PLengthDistributionPlot(Plot):
    PLOT_NAME = "Protein Length Distribution"
    SOURCE_TYPE = "mdisorder pbased"
    FILE_NAME = "p_length_hist"

    def _run(self, df: pd.DataFrame):
        # Add KL-heatmap to histogram
        gs = gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[2, 1])
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])
        fig = plt.gcf()
        fig.set_size_inches(9.5, 4.5)

        arg = "protein length"
        bins = np.arange(start=0, stop=2540, step=40)
        ax1 = sns.histplot(
            data=df,
            x=arg,
            hue="proteome",
            palette=self.get_color_scheme(),
            kde=True,
            kde_kws={
                "bw_adjust": 0.3,
                "gridsize": 2000,
            },
            line_kws={"lw": 1.5},
            bins=bins,
            common_norm=False,
            stat="proportion",
            alpha=0.5,
            edgecolor=(1.0, 1.0, 1.0, 0.5),
            linewidth=1.0,
            ax=ax1,
        )
        self.add_mean_to_legend(df, "protein length", ax1)
        ax1.set_xlim(-10, 2500)
        ax1.set_ylim(0.0, 0.12)

        # Calculate CIs and their positions
        all_cis = ci_per_bin(df, arg, bins)
        all_ys = val_per_bin(df, arg, bins)

        # Calculate centers of bins for plotting
        half = float(bins[1] - bins[0]) / 2
        bin_centers = [i + half for i in bins]

        # Plot the errorbars
        plot_errorbars(all_cis, bin_centers, all_ys, self.get_color_scheme(), ax1)

        # Calculate and plot KL
        df_kl = kl_via_binning(df, arg, bins)
        plot_kl(df_kl, self.get_proteome_names(), ax2)

        # Set true plot title
        ax1.set_title(self.PLOT_NAME)

        fig.tight_layout()

    def set_title(self):
        """Overwrites title function for plots including a KL part."""
        pass


class RLengthDistributionPlot(Plot):
    RELATIVE: bool = True
    MINLENGTH: float = -0.02
    MAXLENGTH: float = 1.0
    STEPSIZE: int = 0.02
    BW_ADJUST: float = 0

    # Maybe move down to child classes
    def _run(self, df: pd.DataFrame):
        # Add KL-heatmap to histogram
        gs = gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[2, 1])
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])
        fig = plt.gcf()
        fig.set_size_inches(9.5, 4.5)

        arg = "rel reg length" if self.RELATIVE else "reg length"
        bins = np.arange(
            start=0, stop=self.MAXLENGTH + (0.5 * self.STEPSIZE), step=self.STEPSIZE
        )
        df = df.loc[df[arg] <= self.MAXLENGTH]

        ax1 = sns.histplot(
            data=df,
            x=arg,
            hue="proteome",
            palette=self.get_color_scheme(),
            bins=bins,
            kde=True,
            kde_kws={
                "bw_adjust": self.BW_ADJUST,
                "gridsize": 100,
                "clip": (self.MINLENGTH, self.MAXLENGTH),
            },
            common_norm=False,
            stat="proportion",
            edgecolor=(1.0, 1.0, 1.0, 0.5),
            linewidth=1.0,
            alpha=0.4,
            discrete=False,
            ax=ax1,
        )

        self.add_mean_to_legend(df, arg, ax1)
        ylim = ax1.get_ylim()
        ax1.set_xlabel("Relative Region Length" if self.RELATIVE else "Region Length")
        ax1.set_xlim(0, self.MAXLENGTH)
        ax1.set_ylim(bottom=0.0, top=(ylim[1] + 0.08))

        # Calculate CIs and their positions
        all_cis = ci_per_bin(df, arg, bins)
        all_ys = val_per_bin(df, arg, bins)

        # Calculate centers of bins for plotting
        half = float(bins[1] - bins[0]) / 2
        bin_centers = [i + half for i in bins]

        # Plot the errorbars
        plot_errorbars(all_cis, bin_centers, all_ys, self.get_color_scheme(), ax1)

        # Calculate and plot KL
        df_kl = kl_via_binning(df, arg, bins)
        plot_kl(df_kl, self.get_proteome_names(), ax2)

        # Set true plot title
        ax1.set_title(self.PLOT_NAME)

        fig.tight_layout()

    def set_title(self):
        """Overwrites title function for plots including a KL part."""
        pass


class RLengthDistributionPlotAbsMdisorder(RLengthDistributionPlot):
    PLOT_NAME = "Disordered Region Length Distribution"
    SOURCE_TYPE = "mdisorder rbased"
    FILE_NAME = "mdisorder_r_length_hist_abs"
    RELATIVE = False
    MINLENGTH = 30
    MAXLENGTH = 250
    STEPSIZE = 7
    BW_ADJUST: float = 0.3


class RLengthDistributionPlotAbsTmseg(RLengthDistributionPlot):
    PLOT_NAME = "TMH Length Distribution"
    SOURCE_TYPE = "tmseg rbased"
    FILE_NAME = "tmseg_r_length_hist_abs"
    RELATIVE = False
    MINLENGTH = 12
    MAXLENGTH = 35
    STEPSIZE = 1
    BW_ADJUST: float = 1.0


class RLengthDistributionPlotAbsProna(RLengthDistributionPlot):
    PLOT_NAME = "Binding Region Length Distribution"
    SOURCE_TYPE = "prona rbased"
    FILE_NAME = "prona_r_length_hist_abs"
    RELATIVE = False
    MINLENGTH = 6
    MAXLENGTH = 50
    STEPSIZE = 1
    BW_ADJUST: float = 0.7


class RLengthDistributionPlotRelMdisorder(RLengthDistributionPlot):
    PLOT_NAME = "Relative Disordered Region Length Distribution"
    SOURCE_TYPE = "mdisorder rbased"
    FILE_NAME = "mdisorder_r_length_hist_rel"
    BW_ADJUST: float = 0.6


class RLengthDistributionPlotRelTmseg(RLengthDistributionPlot):
    PLOT_NAME = "Relative TMH Length Distribution"
    SOURCE_TYPE = "tmseg rbased"
    FILE_NAME = "tmseg_r_length_hist_rel"
    BW_ADJUST: float = 0.8


class RLengthDistributionPlotRelProna(RLengthDistributionPlot):
    PLOT_NAME = "Relative Binding Region Length Distribution"
    SOURCE_TYPE = "prona rbased"
    FILE_NAME = "prona_r_length_hist_rel"
    BW_ADJUST: float = 0.8


class RLengthsPlotAbsMdisorderProna(Plot):
    PLOT_NAME = "Distribution of Lengths of Disordered PBRs"
    SOURCE_TYPE = "mdisorder rbased"
    FILE_NAME = "mixed_mdis_prona_r_dpbr_lengths"
    BW_ADJUST: float = 0.2
    RELATIVE = False
    MINLENGTH = 6
    MAXLENGTH = 51
    STEPSIZE = 5

    def overlap_for_pbr(self, df: pd.DataFrame):
        # One protein in this df -> step through/check for each PBR
        # Note: (1) this far, we used PBR length as number of residues to count (here we have RI)
        #       (2) any overlap counts
        df_pbr = df[df["feature"] == "prona"]
        df_dr = df[df["feature"] == "mdisorder"]
        num_res_in_drs = 0
        res_in_drs_list = []
        for i in df_pbr.index:
            start = df_pbr.at[i, "start"]
            end = df_pbr.at[i, "end"]
            num_overlap = len(df_dr[~((df_dr["start"] > end) | (df_dr["end"] < start))])
            if num_overlap > 0:
                valid_length = end - start + 1
                num_res_in_drs += valid_length
                res_in_drs_list.append(valid_length)
        return num_res_in_drs

    def _run(self, df_mdisorder: pd.DataFrame):
        fig, ax1 = plt.subplots()

        df_prona = self.dataframes["prona rbased"]
        df_sizes = (
            self.dataframes["mdisorder pbased"][["proteome", "protein length"]]
            .groupby("proteome")
            .sum()
        )

        # Fraction of residues in DRs

        df_mdisorder_counts = (
            df_mdisorder[["proteome", "reg length"]].groupby(["proteome"]).sum()
        )
        df_mdisorder_counts = df_mdisorder_counts.join(
            df_sizes, on="proteome", how="left"
        )
        df_mdisorder_counts["mdisorder content"] = (
            df_mdisorder_counts["reg length"] / df_mdisorder_counts["protein length"]
        )

        # Fraction of disordered residues used in protein binding

        df_mdisorder["feature"] = ["mdisorder"] * len(df_mdisorder)
        df_prona["feature"] = ["prona"] * len(df_prona)
        df_all = pd.concat([df_mdisorder, df_prona])
        df_all["start"], df_all["end"] = zip(*df_all["region"])

        # Calculate number of overlapping residues for each proteome
        grouped = df_all.groupby(["proteome", "protein"])

        overlap_series = grouped.apply(func=self.overlap_for_pbr)
        overlap_series.name = "overlap"
        hue = "proteome"
        dpbrs = overlap_series[overlap_series.notnull()]
        dpbrs = (
            dpbrs.apply(pd.Series)
            .reset_index()
            .melt(id_vars=["proteome", "protein"])
            .dropna()[["proteome", "protein", "value"]]
        )

        # bins = [6, 11, 16, 21, 26, 31, 41, 51, 101, max(111, max(dpbrs["value"]))]
        bins = np.arange(
            start=0, stop=self.MAXLENGTH + (0.5 * self.STEPSIZE), step=self.STEPSIZE
        )
        dpbrs = dpbrs.loc[dpbrs["value"] <= self.MAXLENGTH]

        ax1 = sns.histplot(
            data=dpbrs,
            x="value",
            hue=hue,
            palette=self.get_color_scheme(),
            bins=bins,
            # kde=True,
            kde_kws={
                "bw_adjust": self.BW_ADJUST,
                "bw_method": "scott",
                "gridsize": 3000,
                "clip": (self.MINLENGTH, self.MAXLENGTH),
                # "common_grid": True, # not working
            },
            # common_grid=True, # not working
            common_norm=False,
            stat="proportion",
            # edgecolor=(1.0, 1.0, 1.0, 0.5),
            alpha=0.8,
            fill=False,
            edgecolor="k",
            linewidth=1.5,
            # cumulative=True,
            ax=ax1,
        )

        ax1.set_title("Distribution of Lengths of Disordered PBRs")
        ax1.set_xlabel("Binding Region Length")
        ax1.set_xlim(right=self.MAXLENGTH)  # absolute lengths
        ylim = ax1.get_ylim()
        ax1.set_ylim(bottom=0.0, top=(ylim[1] + 0.01))
        ax1.set_xticks(bins)
        # self.rename_legend(ax1)
        self.add_mean_to_legend(dpbrs, "value", ax1)
