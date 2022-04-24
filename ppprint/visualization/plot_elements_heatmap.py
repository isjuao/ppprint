import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import gridspec
from matplotlib.collections import QuadMesh
from matplotlib.text import Text

from ppprint.visualization.plot import Plot


class PBindingElementsPlotProna(Plot):
    SOURCE_TYPE = "prona pbased"
    PLOT_NAME = "Fraction of Proteins with Binding Element"
    FILE_NAME = "prona_p_elements"

    def set_title(self):
        pass

    def _run(self, df: pd.DataFrame):
        proteomes = pd.unique(df["proteome"])
        n = len(proteomes)
        names = self.get_proteome_names()

        fig = plt.figure()
        # 3 columns
        # gs = gridspec.GridSpec(nrows=int((n + 2) / 3), ncols=3)
        # fig.set_size_inches(10, int((n + 2) / 3) * 3.3)
        # 2 columns
        gs = gridspec.GridSpec(nrows=int((n + 1) / 2), ncols=2)
        fig.set_size_inches(7, int((n + 1) / 2) * 3)

        fig.subplots_adjust(top=0.82)
        plt.tight_layout()
        # plt.tight_layout(pad=0, w_pad=0.5, h_pad=0.5)
        # fig.tight_layout(rect=[0, 0, 1, 0.95])

        for i, p in enumerate(proteomes):
            df_curr = df[df["proteome"] == p]

            matrix_abs = pd.DataFrame(
                columns=["Prot bind", "No PBR"], index=["DNA bind", "No DBR"]
            )
            matrix_abs.at["DNA bind", "Prot bind"] = len(
                df_curr[(df_curr["DBR content"] > 0.0) & (df_curr["PBR content"] > 0.0)]
            )
            matrix_abs.at["No DBR", "Prot bind"] = len(
                df_curr[
                    (df_curr["DBR content"] == 0.0) & (df_curr["PBR content"] > 0.0)
                ]
            )
            matrix_abs.at["DNA bind", "No PBR"] = len(
                df_curr[
                    (df_curr["DBR content"] > 0.0) & (df_curr["PBR content"] == 0.0)
                ]
            )
            matrix_abs.at["No DBR", "No PBR"] = len(
                df_curr[
                    (df_curr["DBR content"] == 0.0) & (df_curr["PBR content"] == 0.0)
                ]
            )

            matrix_abs["Total"] = matrix_abs["Prot bind"] + matrix_abs["No PBR"]
            matrix_abs = matrix_abs.append(
                pd.Series(
                    {
                        "Prot bind": sum(matrix_abs["Prot bind"]),
                        "No PBR": sum(matrix_abs["No PBR"]),
                    },
                    name="Total",
                )
            )

            # Sanity checking totals
            total_column = (
                matrix_abs.at["DNA bind", "Total"] + matrix_abs.at["No DBR", "Total"]
            )
            total_row = (
                matrix_abs.at["Total", "Prot bind"] + matrix_abs.at["Total", "No PBR"]
            )
            assert total_row == total_column
            matrix_abs.at["Total", "Total"] = total_column

            matrix_abs = matrix_abs.astype(int)
            matrix_rel = matrix_abs / len(df_curr)
            matrix_rel = matrix_rel.astype(float)

            current_ax = plt.subplot(gs[i])
            sns.heatmap(
                matrix_rel,
                square=True,
                cbar_kws={
                    "shrink": 0.5,
                    "label": "Fraction",
                },
                annot=True,
                annot_kws={"fontsize": 10},
                ax=current_ax,
                cmap="viridis_r",
            )
            current_ax.set_title(
                names[p],
                size=10,
                pad=0,
                y=-0.02,
                bbox=dict(facecolor="grey", alpha=0.2),
            )

            # Make colors of the last column white
            quadmesh = current_ax.findobj(QuadMesh)[0]
            facecolors = quadmesh.get_facecolors()
            facecolors[np.array([2, 5, 6, 7, 8])] = np.array([1.0, 1.0, 1.0, 1.0])
            quadmesh.set_facecolors(facecolors)

            for j in current_ax.findobj(Text):
                pos = j.get_position()
                if pos[0] > 2 or pos[1] > 2:
                    # Black text for total columns
                    j.set_color("black")

            current_ax.xaxis.tick_top()
            current_ax.xaxis.set_label_position("top")

        fig.suptitle("Fraction of Proteins with Binding Element", fontsize=15, y=0.99)
