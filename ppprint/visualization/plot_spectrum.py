import itertools
import random

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt  # TODO: how to get plt object (new one for every plot or clear!)
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import gridspec
from scipy import stats

from ppprint.visualization.plot import Plot


class RSpectrumPlotMdisorder(Plot):
    PLOT_NAME = "Spectrum of Disordered Regions"  # FIXME: Title for whole plot
    SOURCE_TYPE = "mdisorder rbased"
    FILE_NAME = "mdisorder_r_spectrum"

    def collect_lists(self, df: pd.DataFrame):
        """Determines center and distance from start to center for each region and adds info to dataframe."""

        # TODO: Flip x and y
        df["start"], df["end"] = zip(*df["point region"])
        df["y"] = df[["start", "end"]].mean(axis=1).round(decimals=2)
        df["x"] = -(df["y"] - df["start"])
        # Select only columns of interest
        df = df[["x", "y", "proteome"]]

        return df

    def group_and_metrics(self, df_points: pd.DataFrame):
        """Calculates metrics for a spectrum plot. Returns a dataframe with data ('binned' per center), means, SEs."""

        # Calculate stats
        df_grouped = (
            df_points.groupby(["proteome", "y"])
            .agg(
                **{
                    "mean": pd.NamedAgg(column="x", aggfunc="mean"),
                    "se": pd.NamedAgg(column="x", aggfunc=stats.sem),
                }
            )
            .reset_index()
        )

        # Fill up 0 values for missing center values
        multi_index = pd.MultiIndex.from_product(
            [df_grouped["proteome"].unique(), np.arange(0, 1.01, 0.01)],
            names=["proteome", "y"],
        )
        df_grouped = (
            df_grouped.set_index(["proteome", "y"])
            .reindex(multi_index)
            .reset_index()
            .fillna(0)
        )

        return df_grouped

    def plot_spectrum(self, ax1, df_grouped, colors, names, used_proteomes):
        """Plots the spectrum part."""

        patches = []
        for p in used_proteomes:
            # Color
            current_color = colors[p]

            # Plot means and CIs as shadow
            y = -df_grouped[df_grouped["proteome"] == p][
                "mean"
            ]  # TODO: Adjust when flipped
            y_none = [None if val == 0 else val for val in y]
            x = df_grouped[df_grouped["proteome"] == p]["y"]
            ax1.plot(
                x,
                y_none,
                color=current_color,
                marker="o",
                ms=1,
                zorder=10,
                linestyle="-",
                linewidth=1,
            )
            se = df_grouped[df_grouped["proteome"] == p]["se"]
            where_not = [(val is not None) for val in y_none]
            ax1.fill_between(
                x, y - se, y + se, where=where_not, color=current_color, alpha=0.3
            )

            patches.append(mpatches.Patch(color=current_color, label=names[p]))

        ax1.spines["left"].set_position("zero")
        ax1.spines["bottom"].set_position("zero")

        ax1.set_xlim(-0.0, 1.0)
        ax1.set_ylim(-0.0, 0.5)
        xticklabels = np.arange(0.0, 1.05, 0.1)
        ax1.set_xticks(xticklabels)

        ax1.set_xlabel("Position of Center of DR in Protein", loc="center")
        ax1.set_ylabel("Start/end of DR in Protein", loc="center")
        main_legend = ax1.legend(
            handles=patches,
            fontsize="small",
            frameon=True,
        )
        main_legend.get_frame().set_edgecolor("grey")
        main_legend.get_frame().set_linewidth(0.5)
        ax1.add_artist(main_legend)
        ax1.set_title(
            "Width of DRs at Center Positions", pad=7, backgroundcolor=(0, 0, 0, 0.1)
        )

    def set_title(self):
        pass

    def plot_cc(self, names, pairs, paired_cc, ax):
        """Plots the foreground CC plot."""

        chosen_palette = plt.get_cmap("Dark2")
        colors = random.sample(chosen_palette.colors, len(pairs))
        x_end = 1.0
        x = np.arange(
            0 - (x_end / 2), 0 + (x_end / 2) + 0.01, 0.01
        )  # unsure but apparently one value too much
        # x = np.arange(0 - (x_end / 2), 0 + (x_end / 2), 0.01)
        patches = []

        for i, pair in enumerate(pairs):
            cross_corr = paired_cc[i]
            sns.lineplot(y=cross_corr, x=x, ax=ax, color=colors[i], alpha=0.8)

            peak = max(cross_corr)
            peak_index = (np.argmax(cross_corr) - 50) / 100

            patches.append(
                mpatches.Patch(
                    color=colors[i],
                    label=f"{[names[pair[0]], names[pair[1]]]} [max: {round(peak, 2)} at {peak_index}]",
                )
            )

        x_ticks = np.arange(0 - (x_end / 2), 0 + (x_end / 2 + 0.01), 0.1)
        ax.xaxis.set_ticks(x_ticks)
        ax.set_ylim(
            bottom=0.0,
        )  # top=(peak + 0.1))
        ax.set_title(
            f"Cross-Correlation of\n {list(names.values())}",
            fontsize=10,
            backgroundcolor=(0, 0, 0, 0.1),
        )
        cc_legend = ax.legend(
            handles=patches,
            loc="lower center",
            fontsize="small",
            frameon=True,
        )
        cc_legend.get_frame().set_edgecolor("grey")
        cc_legend.get_frame().set_linewidth(0.5)
        ax.add_artist(cc_legend)

    def plot_background(self, ax):
        """Plots the background zones for the CC part."""

        x = np.linspace(-0.5, 0.5, 1000)
        fill_color = (0, 0, 0, 0.05)

        # Between mixed and euk/euk
        y = -2 * abs(x) + 1.5
        ax.fill_between(x, 0, y, color=fill_color, zorder=0)
        # Between prok/prok and mixed
        y = -3 * abs(x) + 2.5
        ax.fill_between(x, 0, y, color=fill_color, zorder=0)

        extra_legend = ax.legend(
            handles=[
                mpatches.Patch(color="white", label="Prokaryotic Pairings"),
                mpatches.Patch(color=fill_color, label="Mixed Zone"),
                mpatches.Patch(
                    color=tuple(np.add(fill_color, fill_color)),
                    label="Eukaryotic Pairings",
                ),
            ],
            loc="upper right",
            fontsize="small",
            frameon=True,
        )
        extra_legend.get_frame().set_edgecolor("grey")
        extra_legend.get_frame().set_linewidth(0.5)
        ax.add_artist(extra_legend)

    def _run(self, df: pd.DataFrame):
        """Creates a spectrum plot."""

        # -- Preprocessing --

        proteomes = pd.unique(df["proteome"])
        n = len(proteomes)
        assert n <= 4  # TODO: Ensure this when creating a VisualizationJob

        # Get points df
        df_points = self.collect_lists(df)

        # Calculate SE
        df_grouped = self.group_and_metrics(df_points)

        # # Add pseudo counts
        # df_grouped["mean"] -= math.exp(-12)

        def get_means(df, proteome):
            """For a given proteome, gets the mean column as a series with a reset index."""
            return df[df["proteome"] == proteome]["mean"].reset_index(drop=True)

        proteome_pairs = list(itertools.combinations(proteomes, 2))
        all_cross_corr = []
        for pair in proteome_pairs:
            # Calculate cross-correlation
            cross_corr = get_means(df_grouped, pair[0]).corr(
                get_means(df_grouped, pair[1]),
                method=(lambda a, b: np.correlate(a, b, "same")),
            )
            all_cross_corr.append(cross_corr)

        # -- Plotting --

        fig = plt.figure()
        gs = gridspec.GridSpec(2, 1, height_ratios=[6, 11])
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])

        # Plot paired CC
        names = self.get_proteome_names()
        self.plot_cc(names, proteome_pairs, all_cross_corr, ax2)

        # Insert background into CC plot
        self.plot_background(ax2)

        # TODO: Maybe assign colors based on kingdom
        # kingdom_colors = get_kingdom_colors()
        # colors = [kingdom_colors[v] for v in quadruple]
        # alphas = self.calc_alphas(df, combination)

        # Plot spectrum
        fig.set_size_inches(6, 9)
        fig.tight_layout(pad=2)
        colors = self.get_color_scheme()
        self.plot_spectrum(ax1, df_grouped, colors, names, proteomes)
