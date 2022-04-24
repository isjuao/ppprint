from abc import ABC, abstractmethod

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import figure, gridspec

from ppprint.visualization.plot import Plot


class PiePlotTmseg(Plot, ABC):
    SOURCE_TYPE = "tmseg pbased"

    @abstractmethod
    def get_pie(self, df):
        """Generates pie values and labels for a pie chart."""
        pass

    def set_title(self):
        pass

    @abstractmethod
    def set_suptitle(self, fig):
        """Pie chart dashboards of multiple subplots require setting a suptitle."""
        pass

    def _run(self, df: pd.DataFrame):
        proteomes = pd.unique(df["proteome"])
        n = len(proteomes)
        names = self.get_proteome_names()

        gs = gridspec.GridSpec(nrows=int((n + 1) / 2), ncols=2)
        fig = plt.figure()
        fig.set_size_inches(11, int((n + 1) / 2) * 4)
        fig.tight_layout()

        for i, p in enumerate(proteomes):
            df_curr = df[df["proteome"] == p]

            pie, pie_labels = self.get_pie(df_curr)

            current_ax = plt.subplot(gs[i])
            _, _, autotexts = current_ax.pie(
                pie,
                labels=pie_labels,
                autopct="%.2f%%",
                colors=sns.color_palette("pastel"),
                textprops={"size": 18},
                labeldistance=1.08,
                wedgeprops={"alpha": 0.8},
            )
            plt.setp(autotexts, color="white")
            current_ax.set_title(
                names[p],
                size=18,
                pad=0,
                y=-0.05,
                bbox=dict(facecolor="grey", alpha=0.2),
            )
        fig.subplots_adjust(top=0.92)
        self.set_suptitle(fig)


class POrientationsPlotTmseg(PiePlotTmseg):
    PLOT_NAME = "Orientation (Location of N-Terminus) of all TMPs"
    FILE_NAME = "tmseg_p_orientations"

    def set_suptitle(self, fig: figure):
        fig.suptitle(
            "Orientation (Location of N-Terminus) of all TMPs", fontsize=23, y=0.99
        )

    def get_pie(self, df_curr: pd.DataFrame):
        # Drop globular porteins
        abs_pie = df_curr["orientation"].value_counts().drop(0)
        # Get pie fractions
        num_TMPs = sum(abs_pie.values)
        pie = abs_pie / num_TMPs
        # Merge two categories
        pie["Membrane/SP"] = pie["Signal Peptide"] + pie["Membrane"]
        pie = pie.drop(["Signal Peptide", "Membrane"])

        pie_labels = pie.index

        return pie, pie_labels


class PProtClassPlotTmseg(PiePlotTmseg):
    PLOT_NAME = "Protein Classes"
    FILE_NAME = "tmseg_p_prot_classes"

    def set_suptitle(self, fig: figure):
        fig.suptitle("(Transmembrane-) Protein Classes", fontsize=23, y=0.99)

    def get_pie(self, df_curr: pd.DataFrame):
        df_count = pd.DataFrame(df_curr.groupby(["number of regions"]).size())

        num_tmps = sum(df_count.iloc[1:][0])
        num_glob = sum(df_count.iloc[0])
        num_single = sum(df_count.iloc[1])
        num_multi = sum(df_count.iloc[2:][0])

        pie = np.array([num_glob, num_single, num_multi])
        pie = pie / sum(pie)
        pie_labels = ["Globular", "Single-pass", "Multi-pass"]

        return pie, pie_labels
