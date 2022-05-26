import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from ppprint.visualization.plot import Plot


class RScatterPlotMdisorderProna(Plot):
    SOURCE_TYPE = "mdisorder rbased"
    PLOT_NAME = "Relative Overlap of TMPs and Disordered Proteins"
    FILE_NAME = "mixed_mdis_prona_r_scatter"

    def overlap_for_pbr(self, df: pd.DataFrame):
        # One protein in this df -> step through/check for each PBR
        # Note: (1) this far, we used PBR length as number of residues to count (here we have RI)
        #       (2) any overlap counts
        df_pbr = df[df["feature"] == "prona"]
        df_dr = df[df["feature"] == "mdisorder"]
        num_res_in_drs = 0
        res_in_drs_list = []
        for i in df_pbr.index:
            start = df_pbr.at[i, "start"]
            end = df_pbr.at[i, "end"]
            num_overlap = len(df_dr[~((df_dr["start"] > end) | (df_dr["end"] < start))])
            if num_overlap > 0:
                valid_length = end - start + 1
                num_res_in_drs += valid_length
                res_in_drs_list.append(valid_length)
        return num_res_in_drs

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

        overlap_series = grouped.apply(func=self.overlap_for_pbr)
        overlap_series.name = "overlap"
        grouped = (
            grouped.count().join(overlap_series, how="right").groupby("proteome").sum()
        )
        # Calculate relative
        # (1) Normalized by number of residues in proteome
        # df_prona_counts = grouped[["overlap"]].join(df_sizes)
        # df_prona_counts["overlap content"] = df_prona_counts["overlap"] / df_prona_counts["protein length"]
        # (2) Normalized by number of disordered residues
        df_prona_counts = grouped[["overlap"]].join(
            df_mdisorder_counts[["reg length"]], how="left"
        )
        df_prona_counts["overlap content"] = (
            df_prona_counts["overlap"] / df_prona_counts["reg length"]
        )
        df_counts = df_mdisorder_counts[["mdisorder content"]].join(
            df_prona_counts[["overlap content"]]
        )

        # plot
        sns.scatterplot(
            data=df_counts,
            x="mdisorder content",
            y="overlap content",
            hue="proteome",
            alpha=1,
            ax=ax1,
            palette=self.get_color_scheme(),
            s=100,
        )
        self.rename_legend(ax1)
        ax1.set_ylim(bottom=0.0)
        ax1.set_ylabel("Fraction of DR Residues Used in Binding")
        ax1.set_xlabel("Fraction of Residues in DRs")
        ax1.set_title("Relationship of Disorder/Protein Binding", fontsize=11, y=0.99)
