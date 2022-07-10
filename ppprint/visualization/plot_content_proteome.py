from abc import ABC
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from ppprint.visualization.plot import Plot


class RContentPerProteomePlot(Plot, ABC):
    """Requires r-based dataframe."""

    def _run_old(self, df: pd.DataFrame):
        """
        Defines proteome content differently: counts protein lengths as often as regions occur:
        Multiple times for multiple regions, 0 times for proteins without regions.
        """
        ax1 = plt.subplot()
        names = self.get_proteome_names()
        colors = self.get_color_scheme()

        df_counts = df.groupby(["proteome"]).sum().reset_index()
        df_counts["DR content"] = df_counts["reg length"] / df_counts["protein length"]

        sns.barplot(
            x="proteome", y="DR content", data=df_counts, ax=ax1, palette=colors
        )

        ax1.set_xlabel("Proteome")
        ax1.set_ylabel("Fraction of Region Residues")
        ax1.set_xticklabels(labels=names.values(), fontsize=7)


class PContentPerProteomePlot(Plot):
    def _run(self, df: pd.DataFrame):
        ax1 = plt.subplot()
        names = self.get_proteome_names()
        colors = self.get_color_scheme()

        gb = pd.DataFrame(df.groupby("proteome")["protein length"].sum()).rename(
            columns={"protein length": "sum protein length"}
        )
        df = df.join(gb, on="proteome", how="left")
        # Regenerate region lengths
        df["region content"] = df["region content"] * df["protein length"]
        # Calculate content
        df["proteome content"] = df["region content"] / df["sum protein length"]

        sns.barplot(
            x="proteome",
            y="proteome content",
            data=df,
            estimator=sum,
            ci=95,
            n_boot=1000,
            errwidth=1,
            capsize=0.1,
            palette=colors,
            ax=ax1,
        )
        ax1.set_xlabel("Proteome")
        ax1.set_ylabel("Fraction of Region Residues")
        ax1.set_xticklabels(labels=names.values(), fontsize=7)


class PContentPerProteomePlotMdisorder(PContentPerProteomePlot):
    PLOT_NAME = "Disorder Content Per Proteome"
    SOURCE_TYPE = "mdisorder pbased"
    FILE_NAME = "mdisorder_p_content_proteome"


class PContentPerProteomePlotTmseg(PContentPerProteomePlot):
    PLOT_NAME = "TMH Content Per Proteome"
    SOURCE_TYPE = "tmseg pbased"
    FILE_NAME = "tmseg_p_content_proteome"


class PContentPerProteomePlotProna(PContentPerProteomePlot):
    PLOT_NAME = "Binding Content Per Proteome"
    SOURCE_TYPE = "prona pbased"
    FILE_NAME = "prona_p_content_proteome"


class PCompositionPerProteomePlotMdisorder(Plot):
    PLOT_NAME = "Disorder Composition"
    SOURCE_TYPE = "mdisorder pbased"
    FILE_NAME = "mdisorder_p_composition"

    def _run(self, df: pd.DataFrame):
        ax1 = plt.subplot()
        names = self.get_proteome_names()
        colors = self.get_color_scheme()

        df["composition"] = (df["number of regions"] >= 1).replace({True: 1, False: 0})
        df_sizes = pd.DataFrame(df.groupby("proteome").size()).rename(
            columns={0: "proteome size"}
        )
        df = df.join(df_sizes, how="left", on=["proteome"])
        df["composition"] = df["composition"] / df["proteome size"]

        sns.barplot(
            x="proteome",
            y="composition",
            data=df,
            palette=colors,
            ax=ax1,
            estimator=sum,
            ci=95,
            n_boot=1000,
            errwidth=1,
            capsize=0.1,
        )

        ax1.set_xlabel("Proteome")
        ax1.set_ylabel("Fraction")
        ax1.set_xticklabels(labels=names.values(), fontsize=8, rotation=20)
