from abc import ABC, abstractmethod
from typing import Dict
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


class Plot(ABC):
    SOURCE_TYPE: str
    PLOT_NAME: str
    FILE_NAME: str

    base_folder: Path
    dataframes: Dict[str, pd.DataFrame]
    proteome_mapping: Dict[int, str]

    def __init__(self, dataframes: Dict[str, pd.DataFrame], proteome_mapping: Dict[int, str], base_folder: Path):
        self.dataframes = dataframes
        self.proteome_mapping = proteome_mapping
        self.base_folder = base_folder

    @abstractmethod
    def _run(self, df: pd.DataFrame):
        pass

    def get_df(self):
        return self.dataframes[self.SOURCE_TYPE]

    def run(self):
        self._run(self.get_df())
        self.set_title()
        self.store_plot()

    def set_title(self):
        plt.title = self.PLOT_NAME

    def store_plot(self):
        out_path = self.base_folder / f"{self.FILE_NAME}.png"
        plt.savefig(out_path, dpi=200)

    def get_color_scheme(self):
        colors = [(0.8901960784313725, 0.4666666666666667, 0.7607843137254902),
                  (0.5333333333333333, 0.37254901960784315, 0.6392156862745098),
                  (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
                  (0.23529411764705882, 0.5058823529411764, 0.7411764705882353),
                  (0.5490196078431373, 0.10980392156862745, 0.13333333333333333),
                  (0.9098039215686274, 0.47058823529411764, 0.2),
                  (0.7372549019607844, 0.7411764705882353, 0.13333333333333333),
                  (0.5176470588235295, 0.8784313725490196, 0.7764705882352941),
                  (0.06274509803921569, 0.1803921568627451, 0.3607843137254902)]
        palette = {key: colors[i] for i, key in enumerate(self.proteome_mapping.keys())}
        return palette

    def rename_legend(self, ax):
        leg = ax.get_legend()
        # proteome pk order in legend = proteome pk order in df = proteome pk order in mapping
        ax.legend(handles=leg.legendHandles, labels=self.proteome_mapping.values())


class LengthDistributionPlot(Plot):
    SOURCE_TYPE = "mdisorder pbased"
    PLOT_NAME = "Length Distribution"
    FILE_NAME = "mdisorder_length_dist"

    def _run(self, df: pd.DataFrame):
        ax1 = plt.subplot()

        arg = "protein length"
        discrete = True
        bins = np.arange(start=0, stop=2540, step=40)
        ax1 = sns.histplot(data=df, x=arg, hue="proteome", palette=self.get_color_scheme(),
                           kde=True,
                           kde_kws={"bw_adjust": .5, "gridsize": 2000
                                    },
                           bins=bins,
                           discrete=discrete,
                           # cumulative=True,
                           common_norm=False,
                           stat="proportion",
                           # multiple="stack"
                           ax=ax1,
                           )
        self.rename_legend(ax1)
        ax1.set_xlim(-10, 2500)
        ax1.set_ylim(0.0, 0.12)
