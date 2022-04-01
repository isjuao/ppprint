from abc import ABC, abstractmethod

import matplotlib.pyplot as plt  # TODO: how to get plt object (new one for every plot or clear!)
import numpy as np
import pandas as pd
import seaborn as sns

from ppprint.visualization.plot import Plot

# TODO: what about KL -> show optionally? always?


class PLengthDistributionPlot(Plot):
    PLOT_NAME = "Protein Length Distribution"
    SOURCE_TYPE = "mdisorder pbased"
    FILE_NAME = "p_length_hist"

    def _run(self, df: pd.DataFrame):
        ax1 = plt.subplot()

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
    MAXLENGTH = 50  # TODO search best cutoff
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
