import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from ppprint.visualization.plot import Plot


class RPointLinePlot(Plot):
    def split_point_regions(self, df: pd.DataFrame):
        """Extracts all points covered by the given regions in point format."""

        # To ensure a consecutive index
        df = df.reset_index(drop=True)

        touched_points = []
        points_proteomes = []
        for i, pr in enumerate(df["point region"]):
            # Append each touched point (round to 2 decimals)
            t = [
                round(x, 2)
                for x in list(
                    np.arange(round(pr[0], 2), round(pr[1], 2) + 0.001, step=0.01)
                )
            ]
            touched_points.extend(t)
            # Add proteome information
            p = df.at[i, "proteome"]
            points_proteomes.extend([p for _ in range(0, len(t))])

        df_result = pd.DataFrame(
            {"touched point": touched_points, "proteome": points_proteomes}
        )
        return df_result

    def _run(self, df: pd.DataFrame):
        ax1 = plt.subplot()

        df_splitreg = self.split_point_regions(df)

        sizes = pd.DataFrame(df_splitreg.groupby(["proteome"]).size()).rename(
            {0: "proteome size"}, axis=1
        )
        df_splitreg = df_splitreg.join(sizes, how="left", on="proteome")
        df_splitreg["value"] = 1 / df_splitreg["proteome size"]
        sns.lineplot(
            data=df_splitreg,
            x="touched point",
            y="value",
            estimator=sum,
            hue="proteome",
            palette=self.get_color_scheme(),
            err_style="bars",
            ci=95,
            n_boot=5,
            ax=ax1,
        )

        # plt.figure().set_size_inches(7.5, 4.8)

        ax1.set_xlabel("Point in Protein")
        ax1.set_xlim(0.0, 1.0)
        ax1.set_ylabel("Frequency of Disorder at Point")
        ax1.set_ylim(bottom=0.0)
        self.rename_legend(ax1)


class RPointLinePlotMdisorder(RPointLinePlot):
    PLOT_NAME = "Spread of all Disordered Regions"
    SOURCE_TYPE = "mdisorder rbased"
    FILE_NAME = "mdisorder_r_points"


class RPointLinePlotTmseg(RPointLinePlot):
    PLOT_NAME = "Spread of all TMH Regions"
    SOURCE_TYPE = "tmseg rbased"
    FILE_NAME = "tmseg_r_points"


class RPointLinePlotProna(RPointLinePlot):
    PLOT_NAME = "Spread of all Binding Regions"
    SOURCE_TYPE = "prona rbased"
    FILE_NAME = "prona_r_points"
