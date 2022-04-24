from abc import abstractmethod
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from ppprint.visualization.plot import Plot


class PNumberOfRegions(Plot):
    MAXLENGTH: int

    def _run(self, df: pd.DataFrame):
        ax1 = plt.subplot()
        colors = self.get_color_scheme()

        start, stop = df["number of regions"].min(), df["number of regions"].max()
        bin_edges = np.arange(start - 0.5, stop + 1.5)
        bins = bin_edges
        sns.histplot(
            data=df,
            x="number of regions",
            hue="proteome",
            palette=colors,
            bins=bins,
            common_norm=False,
            stat="proportion",
            fill=False,
            edgecolor="k",
            linewidth=1.5,
            shrink=0.5,
            discrete=True,
            alpha=0.7,
            ax=ax1,
        )

        half = float(bins[1] - bins[0]) / 2
        bin_centers = [int(i + half) for i in bins]
        ax1.set_xticks(bin_centers)
        ax1.set_xticklabels(bin_centers)
        ax1.set_xlabel("Number of Regions")
        ax1.set_xlim(-0.5, self.MAXLENGTH)
        ylim = ax1.get_ylim()
        ax1.set_ylim(bottom=0.0, top=(ylim[1] + 0.025))

        self.rename_legend(ax1)


class PNumberOfRegionsMdisorder(PNumberOfRegions):
    PLOT_NAME = "Distribution of Number of DRs Per Protein"
    SOURCE_TYPE = "mdisorder pbased"
    FILE_NAME = "mdisorder_p_num_regions"
    MAXLENGTH = 11


class PNumberOfRegionsTmseg(PNumberOfRegions):
    PLOT_NAME = "Distribution of Number of TMHs Per Protein"
    SOURCE_TYPE = "tmseg pbased"
    FILE_NAME = "tmseg_p_num_regions"
    MAXLENGTH = 16


class PNumberOfRegionsProna(PNumberOfRegions):
    PLOT_NAME = "Distribution of Number of PBRs Per Protein"
    SOURCE_TYPE = "prona pbased"
    FILE_NAME = "prona_p_num_regions"
    MAXLENGTH = 5


class PNumberOfRegTopoTmseg(Plot):
    PLOT_NAME = "Distribution of Number and Orientation of TMHs Per TMP"
    SOURCE_TYPE = "tmseg pbased"
    FILE_NAME = "tmseg_p_num_regions_topo"
    colors: Dict[int, Tuple]
    bins: np.ndarray
    ax1: plt.Axes

    def plot_half(self, df_half: pd.DataFrame):
        sns.histplot(
            data=df_half,
            x="number of regions",
            hue="proteome",
            palette=self.colors,
            bins=self.bins,
            common_norm=False,
            stat="proportion",
            ax=self.ax1,
            shrink=0.5,
            discrete=True,
            fill=False,
            edgecolor="k",
            linewidth=1.5,
            alpha=0.8,
        )

    def _run(self, df: pd.DataFrame):
        self.ax1 = plt.subplot()
        self.colors = self.get_color_scheme()

        start, stop = df["number of regions"].min(), df["number of regions"].max()
        self.bins = np.arange(start - 0.5, stop + 1.5)
        x_rlim = 16
        # Drop proteins without TMH
        df = df[df["orientation"] != "No region"]

        # Cytoplasmic half

        df_half = df[df["orientation"] == "Cytoplasmic"]
        self.plot_half(df_half)
        y_blim = -self.ax1.get_ylim()[1]

        # Turn the histogram upside down
        for p in self.ax1.patches:
            p.set_height(-p.get_height())
        # for line in self.ax1.lines:
        #     # Turn the KDE curve upside down
        #     line.set_ydata(-line.get_ydata())

        # Extracellular half

        df_half = df[df["orientation"] == "Extracellular"]
        self.plot_half(df_half)
        y_ulim = self.ax1.get_ylim()[1]

        self.ax1.set_yticks(np.arange(0.0, (y_ulim + 0.1), 0.1))
        pos_ticks = np.array([t for t in self.ax1.get_yticks() if t > 0])
        ticks = np.concatenate([-pos_ticks[::-1], [0], pos_ticks])
        self.ax1.set_yticks(ticks)
        self.ax1.set_yticklabels([f"{abs(t):.1f}" for t in ticks])
        self.ax1.set_ylim(bottom=(y_blim - 0.05), top=(y_ulim + 0.05))

        self.ax1.text(x=0.2, y=y_ulim, s="N-terminus extracellular", fontsize="x-small")
        self.ax1.text(
            x=0.2, y=(y_blim - 0.03), s="N-terminus cytoplasmic", fontsize="x-small"
        )
        self.ax1.set_xlim(0.0, x_rlim)
        self.ax1.axhline(y=0.0, color="black", lw=0.7, zorder=4)
        self.ax1.set_xlabel("Number of Regions")
        self.ax1.set_xticks(np.arange(0, x_rlim, 1))

        self.rename_legend(self.ax1)
