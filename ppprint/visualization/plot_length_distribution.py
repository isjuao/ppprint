from abc import ABC, abstractmethod

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from ppprint.visualization.plot import Plot


class PLengthDistributionPlot(Plot):
    PLOT_NAME = "Protein Length Distribution"
    SOURCE_TYPE = "mdisorder pbased"
    FILE_NAME = "p_length_hist"

    def _run(self, df: pd.DataFrame):
        ax1 = plt.subplot()

        bins = np.arange(start=0, stop=2540, step=40)
        ax1 = sns.histplot(
            data=df,
            x="protein length",
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
        self.rename_legend(ax1)
        ax1.set_xlim(-10, 2500)
        ax1.set_ylim(0.0, 0.12)


class RLengthDistributionPlot(Plot):
    RELATIVE: bool = True
    MINLENGTH: float = -0.02
    MAXLENGTH: float = 1.0
    STEPSIZE: int = 0.02
    BW_ADJUST: float = 0

    # Maybe move down to child classes
    def _run(self, df: pd.DataFrame):
        ax1 = plt.subplot()

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

        # plt.figure().set_size_inches(6.4, 5.5)

        self.rename_legend(ax1)
        ylim = ax1.get_ylim()
        ax1.set_xlabel("Relative Region Length" if self.RELATIVE else "Region Length")
        ax1.set_xlim(0, self.MAXLENGTH)
        ax1.set_ylim(bottom=0.0, top=(ylim[1] + 0.08))


class RLengthDistributionPlotAbsMdisorder(RLengthDistributionPlot):
    PLOT_NAME = "Disordered Region Length Distribution"
    SOURCE_TYPE = "mdisorder rbased"
    FILE_NAME = "mdisorder_r_length_hist_abs"
    RELATIVE = False
    MINLENGTH = 30
    MAXLENGTH = 250
    STEPSIZE = 2
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
        self.rename_legend(ax1)
