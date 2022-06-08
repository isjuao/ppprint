import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import pandas as pd
import seaborn as sns

from ppprint.visualization.plot import Plot
from ppprint.visualization.plot_extras import (
    ci_per_bin,
    val_per_bin,
    plot_errorbars,
    kl_via_binning,
    plot_kl,
)


class PContentPerProteinPlot(Plot):
    def _run(self, df: pd.DataFrame):
        # Add KL-heatmap to histogram
        gs = gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[2, 1])
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])
        fig = plt.gcf()
        fig.set_size_inches(9.5, 4.5)

        arg = "region content"
        # Filter to region-containing proteins only
        df = df.loc[df["region content"] > 0]
        bins = np.arange(start=0, stop=1.01, step=0.02)
        ax1 = sns.histplot(
            data=df,
            x=arg,
            hue="proteome",
            palette=self.get_color_scheme(),
            kde=True,
            kde_kws={
                "bw_adjust": 0.5,
            },
            bins=bins,
            # cumulative=True,
            common_norm=False,
            stat="proportion",
            # multiple="stack"
            ax=ax1,
            edgecolor=(1.0, 1.0, 1.0, 0.5),
            linewidth=1.0,
            alpha=0.4,
            hue_order=self.get_proteome_names(),
        )
        ax1.set_xlabel("Region Content")
        ax1.set_xlim(-0.02, 1.0)
        ylim = ax1.get_ylim()
        ax1.set_ylim(bottom=0.0, top=(ylim[1] + 0.025))

        self.add_mean_to_legend(df, arg, ax1)

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


class PContentPerProteinPlotMdisorder(PContentPerProteinPlot):
    PLOT_NAME = "Distribution of DR Content Per Disordered Protein"
    SOURCE_TYPE = "mdisorder pbased"
    FILE_NAME = "mdisorder_p_content_protein"


class PContentPerProteinPlotTmseg(PContentPerProteinPlot):
    PLOT_NAME = "Distribution of TM Content Per TM Protein"
    SOURCE_TYPE = "tmseg pbased"
    FILE_NAME = "tmseg_p_content_protein"


class PContentPerProteinPlotProna(PContentPerProteinPlot):
    PLOT_NAME = "Distribution of PBR Content Per Protein-Binding Protein"
    SOURCE_TYPE = "prona pbased"
    FILE_NAME = "prona_p_content_protein"
