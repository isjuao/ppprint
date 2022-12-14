from abc import abstractmethod
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from ppprint.visualization.plot import Plot


class PResidueFractions(Plot):
    @abstractmethod
    def configure(self, df: pd.DataFrame):
        pass

    @abstractmethod
    def plotting(self, df_long: pd.DataFrame, ax1: plt.Axes):
        pass

    def _run(self, df: pd.DataFrame):
        ax1 = plt.subplot()

        # Get correct df and classes
        df, value_vars = self.configure(df)
        df = df.drop(columns=["protein length", "number of regions"], axis=1)
        df_long = df.melt(
            id_vars=["proteome"], value_vars=value_vars, var_name="element"
        )

        self.plotting(df_long, ax1)

        ax1.set_ylabel("Frequency")
        ax1.set_xlabel("Proteome")
        names = self.get_proteome_names()
        ax1.set_xticklabels(labels=names.values(), fontsize=7)
        ax1.legend(
            bbox_to_anchor=(1.01, 1),
            borderaxespad=0.0,
            frameon=True,
            fancybox=True,
            framealpha=0.8,
            facecolor="white",
        )


class PResidueFractionsTmseg(PResidueFractions):
    PLOT_NAME = "Fraction of TMP Residues (I)nside/(M)embrane/(O)utside"
    SOURCE_TYPE = "tmseg pbased"

    def configure(self, df: pd.DataFrame):
        df = df[df["number of regions"] > 0]
        return df, ["I", "M", "O"]


class PResidueFractionsBarsTmseg(PResidueFractionsTmseg):
    FILE_NAME = "tmseg_p_res_fractions_bars"

    def plotting(self, df_long: pd.DataFrame, ax1: plt.Axes):
        sns.barplot(
            data=df_long,
            x="proteome",
            y="value",
            hue="element",
            palette="Greys_r",
            errwidth=1,
            capsize=0.1,
            estimator=np.mean,
            ci=95,
            n_boot=1000,
            edgecolor=(1.0, 1.0, 1.0, 1.0),
            linewidth=1,
            ax=ax1,
        )


class PResidueFractionsViolinsTmseg(PResidueFractionsTmseg):
    FILE_NAME = "tmseg_p_res_fractions_violins"

    def plotting(self, df_long: pd.DataFrame, ax1: plt.Axes):
        sns.violinplot(
            x="proteome",
            y="value",
            hue="element",
            data=df_long,
            ax=ax1,
            dodge=True,
            alpha=0.1,
            size=4,
            palette="Greys_r",
            cut=0,
        )


class PProtClassFractionsProna(Plot):
    PLOT_NAME = "Fraction of DNA/RNA/Protein Binding Proteins"
    SOURCE_TYPE = "prona pbased"
    FILE_NAME = "prona_p_prot_fractions"

    def _run(self, df: pd.DataFrame):
        ax1 = plt.subplot()

        # Calculate number of X-binding proteins and melt dataframe
        molecule_contents = ["DBR content", "RBR content", "PBR content"]
        for cat in molecule_contents:
            df[cat] = (df[cat] > 0.0).replace({True: 1, False: 0})
        df_long = df.melt(
            id_vars=["proteome"], value_vars=molecule_contents, var_name="molecule"
        )

        # Calculate sizes for relativity
        df_sizes = pd.DataFrame(df.groupby("proteome").size()).rename(
            columns={0: "proteome size"}
        )
        df_long = df_long.join(df_sizes, how="left", on=["proteome"])
        df_long["value"] = df_long["value"] / df_long["proteome size"]

        sns.barplot(
            data=df_long,
            x="proteome",
            y="value",
            hue="molecule",
            palette="rocket",
            errwidth=1,
            capsize=0.1,
            estimator=sum,
            ci=95,
            n_boot=1000,
            edgecolor=(1.0, 1.0, 1.0, 1.0),
            linewidth=1,
            ax=ax1,
        )

        ax1.set_ylabel("Frequency")
        ax1.set_xlabel("Proteome")
        names = self.get_proteome_names()
        ax1.set_xticklabels(labels=names.values(), fontsize=8)
        mol_categories = ["DNA binding", "RNA binding", "Protein binding"]
        leg = ax1.legend()
        ax1.legend(
            handles=leg.legendHandles,
            labels=mol_categories,
            frameon=True,
            fancybox=True,
            framealpha=0.8,
            facecolor="white",
        )


class PResidueFractionsReprof(PResidueFractions):
    PLOT_NAME = "Fraction of Residues H(Helix)/E(Strand)/O(Other)"
    SOURCE_TYPE = "reprof pbased"
    FILE_NAME = "reprof_p_res_fractions_bars"

    def configure(self, df: pd.DataFrame):
        return df, ["H", "E", "O"]

    def plotting(self, df_long: pd.DataFrame, ax1: plt.Axes):
        sns.barplot(
            data=df_long,
            x="proteome",
            y="value",
            hue="element",
            palette="Greys_r",
            errwidth=1,
            capsize=0.1,
            estimator=np.mean,
            ci=95,
            n_boot=1000,
            edgecolor=(1.0, 1.0, 1.0, 1.0),
            linewidth=1,
            ax=ax1,
        )
