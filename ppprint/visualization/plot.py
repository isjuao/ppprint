"""
Provides all basic required functionalities for plotting.
All plots to be used by ppprint have to correspond to subclasses of `Plot`.
"""

import logging
import math
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, Set, Tuple

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)


class Plot(ABC):
    SOURCE_TYPE: str
    PLOT_NAME: str
    FILE_NAME: str

    base_folder: Path
    dataframes: Dict[str, pd.DataFrame]
    proteome_mapping: Dict[int, Tuple[str, Tuple[float, float, float]]]

    def __init__(
        self,
        dataframes: Dict[str, pd.DataFrame],
        proteome_mapping: Dict[int, Tuple[str, Tuple[float, float, float]]],
        base_folder: Path,
    ):
        logger.debug("Initialized Plot object!")
        self.dataframes = dataframes
        self.proteome_mapping = proteome_mapping
        self.base_folder = base_folder

    @abstractmethod
    def _run(self, df: pd.DataFrame):
        pass

    def get_df(self):
        """Returns correct dataframe based on source type."""
        return self.dataframes[self.SOURCE_TYPE]

    def run(self):
        # clear plot TODO object oriented?
        plt.style.use("seaborn-whitegrid")
        plt.clf()
        plt.figure().clf()
        plt.cla()
        self._run(self.get_df())
        self.set_title()
        self.store_plot()
        plt.close(plt.gcf())

    def set_title(self):
        plt.title(self.PLOT_NAME)

    def store_plot(self):
        out_path = self.base_folder / f"{self.FILE_NAME}.png"
        plt.savefig(out_path, dpi=200, bbox_inches="tight")

    def get_color_scheme(self):
        return {key: value[1] for key, value in self.proteome_mapping.items()}

    def get_proteome_names(self):
        return {key: value[0] for key, value in self.proteome_mapping.items()}

    def rename_legend(self, ax):
        leg = ax.get_legend()
        # proteome pk order in legend = proteome pk order in df = proteome pk order in mapping
        ax.legend(
            handles=leg.legendHandles,
            labels=self.get_proteome_names().values(),
            frameon=True,
            fancybox=True,
            framealpha=0.8,
            facecolor="white",
        )

    def add_mean_to_legend(self, df, arg, ax):
        """Calculates the mean and SEM of a column for each proteome and adds information to plot legend."""

        df_means = df[["proteome", arg]].groupby("proteome").agg(["mean", stats.sem])
        df_means.columns = df_means.columns.droplevel()

        patches = []
        colors = self.get_color_scheme()
        names = self.get_proteome_names()

        for p in pd.unique(df["proteome"]):
            mean = df_means.at[p, "mean"]
            sem = df_means.at[p, "sem"]
            conf = -int(math.log10(abs(sem + math.exp(-12)))) + 1
            patches.append(
                mpatches.Patch(
                    color=colors[p],
                    label=f"{names[p]} [mean: {round(mean, conf)} +/- {round(sem, conf)}]",
                )
            )
        ax.legend(
            handles=patches,
            frameon=True,
            fancybox=True,
            framealpha=0.8,
            facecolor="white",
            # fontsize=8,
        )
