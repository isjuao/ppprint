from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from ppprint.visualization.plot import Plot


class PProteomeSizes(Plot):
    PLOT_NAME = "Proteome Sizes"
    SOURCE_TYPE = "mdisorder pbased"
    FILE_NAME = "p_proteome_sizes"

    def _run(self, df: pd.DataFrame):
        ax1 = plt.subplot()

        sizes = df.groupby(["proteome"]).size().to_frame("size")

        names = self.get_proteome_names()
        sizes = sizes.join(
            pd.DataFrame.from_dict(names, orient="index", columns=["name"]), how="left"
        ).sort_values("size", ascending=True)

        sns.barplot(
            data=sizes,
            x="name",
            y="size",
            color="black",
            dodge=False,
            linewidth=2,
            edgecolor="black",
            alpha=0.6,
            ax=ax1,
        )

        ax1.set_title("Proteome Sizes")
        ax1.set_ylabel("Number of Proteins")
        ax1.set_xlabel("Organism")
        ax1.set_xticklabels(labels=ax1.get_xticklabels(), rotation=30)

        # Annotate actual numbers
        for x in ax1.patches:
            ax1.annotate(
                format(x.get_height(), ".0f"),
                (0.4 + x.get_x() + x.get_width() / 2.0, x.get_height()),
                ha="center",
                va="center",
                xytext=(0, 9),
                textcoords="offset points",
            )
        plt.subplots_adjust(top=0.92, bottom=0.25, right=0.95)

    # def store_plot(self):
    #     """Storing routine for thesis document. ^^"""
    #     fig = plt.gcf()
    #     fig.set_size_inches(6.4, 3.5)
    #     out = Path("/home/isabell/work/python/thesis/ppprint/plot_pdfs/")
    #     plt.savefig(out / "proteome_sizes.pdf", bbox_inches="tight")
