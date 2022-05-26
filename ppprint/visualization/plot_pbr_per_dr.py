import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from ppprint.visualization.plot import Plot


class RRegionPlotMdisorderProna(Plot):
    SOURCE_TYPE = "mdisorder rbased"
    PLOT_NAME = "Distribution of Number of PBRs Per DR"
    FILE_NAME = "mixed_mdis_prona_r_pbr_per_dr"

    def overlap_for_dr(self, df: pd.DataFrame):
        # One protein in this df -> step through/check for each DR
        # Note: (1) counting number of PBRs per DR
        #       (2) any overlap counts
        df_pbr = df[df["feature"] == "prona"]
        df_dr = df[df["feature"] == "mdisorder"]
        pbrs_per_dr = []
        for i in df_dr.index:
            start = df_dr.at[i, "start"]
            end = df_dr.at[i, "end"]
            num_overlap = len(
                df_pbr[~((df_pbr["start"] > end) | (df_pbr["end"] < start))]
            )
            pbrs_per_dr.append(num_overlap)
        if not pbrs_per_dr:
            # No DR, do not include
            return None
        return pbrs_per_dr

    def _run(self, df_mdisorder: pd.DataFrame):
        fig, ax1 = plt.subplots()

        df_prona = self.dataframes["prona rbased"]
        df_sizes = (
            self.dataframes["mdisorder pbased"][["proteome", "protein length"]]
            .groupby("proteome")
            .sum()
        )

        # Fraction of residues in DRs

        df_mdisorder_counts = (
            df_mdisorder[["proteome", "reg length"]].groupby(["proteome"]).sum()
        )
        df_mdisorder_counts = df_mdisorder_counts.join(
            df_sizes, on="proteome", how="left"
        )
        df_mdisorder_counts["mdisorder content"] = (
            df_mdisorder_counts["reg length"] / df_mdisorder_counts["protein length"]
        )

        # Fraction of disordered residues used in protein binding

        df_mdisorder["feature"] = ["mdisorder"] * len(df_mdisorder)
        df_prona["feature"] = ["prona"] * len(df_prona)
        df_all = pd.concat([df_mdisorder, df_prona])
        df_all["start"], df_all["end"] = zip(*df_all["region"])

        # Calculate number of overlapping residues for each proteome
        grouped = df_all.groupby(["proteome", "protein"])

        overlap_series = grouped.apply(func=self.overlap_for_dr)
        overlap_series.name = "overlap"
        proteins_with_drs = overlap_series[overlap_series.notnull()]
        drs = (
            proteins_with_drs.apply(pd.Series)
            .reset_index()
            .melt(id_vars=["proteome", "protein"])
            .dropna()[["proteome", "protein", "value"]]
        )

        discrete = True
        # Bins are overwritten by seaborn when plotting hist
        start, stop = drs["value"].min(), drs["value"].max()
        x_rlim = min(6, stop)
        bins = np.arange(start - 0.5, stop + 1.5)
        sns.histplot(
            data=drs,
            x="value",
            hue="proteome",
            palette=self.get_color_scheme(),
            bins=bins,
            fill=False,
            common_norm=False,
            stat="proportion",
            ax=ax1,
            shrink=0.5,
            discrete=discrete,
            alpha=0.6,
        )
        ax1.set_title("Distribution of Number of PBRs Per DR")
        ax1.set_xlabel("Number of Binding Regions Per DR")
        half = float(bins[1] - bins[0]) / 2
        bin_centers = [int(i + half) for i in bins]
        ax1.set_xticks(bin_centers)
        ax1.set_xticklabels(bin_centers)
        ax1.set_xlim(-0.5, x_rlim)
        ylim = ax1.get_ylim()
        ax1.set_ylim(bottom=0.0, top=(ylim[1] + 0.025))
        self.rename_legend(ax1)

        ax1.set_xlim(left=-0.5)
        ax1.set_ylim(top=1.05)
