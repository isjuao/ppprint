import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from ppprint.visualization.plot import Plot


class PContentPerProteinPlot(Plot):
    def _run(self, df: pd.DataFrame):
        ax1 = plt.subplot()

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
        self.rename_legend(ax1)


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
