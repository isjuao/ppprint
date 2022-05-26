import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import gridspec
from matplotlib_venn import venn2

from ppprint.visualization.plot import Plot


class POverlapPlotTmsegMdisorder(Plot):
    SOURCE_TYPE = "tmseg pbased"
    PLOT_NAME = "Relative Overlap of TMPs and Disordered Proteins"
    FILE_NAME = "mixed_tmseg_mdis_p_overlap"

    def set_title(self):
        pass

    def _run(self, df_tmseg: pd.DataFrame):
        df_mdisorder = self.dataframes["mdisorder pbased"]

        proteomes = pd.unique(df_tmseg["proteome"])
        n = len(proteomes)
        proteomes_m = pd.unique(df_mdisorder["proteome"])
        if len(proteomes_m) > n:
            n = len(proteomes_m)
            proteomes = proteomes_m
        names = self.get_proteome_names()
        relative = True

        fig = plt.figure()
        # 3 columns
        # gs = gridspec.GridSpec(nrows=int((n + 2) / 3), ncols=3)
        # fig.set_size_inches(10, int((n + 2) / 3) * 3.3)
        # 2 columns
        gs = gridspec.GridSpec(nrows=int((n + 1) / 2), ncols=2)
        fig.set_size_inches(8, int((n + 1) / 2) * 2.2)
        # fig.set_size_inches(7, int((n + 1) / 2) * 3)

        fig.subplots_adjust(top=0.82)
        plt.tight_layout()
        # plt.tight_layout(pad=0, w_pad=0.5, h_pad=0.5)
        # fig.tight_layout(rect=[0, 0, 1, 0.95])

        for i, p in enumerate(proteomes):
            current_ax = plt.subplot(gs[i])
            df_t = df_tmseg[df_tmseg["proteome"] == p]
            df_m = df_mdisorder[df_mdisorder["proteome"] == p]

            df_m["mdis"] = df_m["number of regions"] > 0
            df_m = df_m[["mdis"]]
            df_t["tm"] = df_t["number of regions"] > 0
            df_t = df_t.join(df_m, how="left")
            df_t["both"] = df_t["tm"] & df_t["mdis"]

            intersect = df_t["both"].value_counts()[True]
            tm = df_t["tm"].value_counts()[True]
            mdis = df_t["mdis"].value_counts()[True]
            total = len(df_t)

            if relative:
                subsets = [
                    round(tm / total, 2),
                    round(mdis / total, 2),
                    round(intersect / total, 2),
                ]
            else:
                subsets = [tm, mdis, intersect]
            venn2(
                subsets=subsets,
                set_labels=["TM", "Disordered", "Both"],
                set_colors=("b", "y"),
            )

            current_ax.set_title(
                names[p], size=11, pad=0, y=-0.2, bbox=dict(facecolor="grey", alpha=0.2)
            )
        if relative:
            fig.suptitle(
                "Relative Overlap of TMPs and Disordered Proteins",
                y=0.9,  # fontsize=15,
            )
        else:
            fig.suptitle(
                "Absolute Overlap of TMPs and Disordered Proteins",
                y=0.9,  # fontsize=15,
            )
