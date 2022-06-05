import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from ppprint.visualization.plot import Plot


class PContentRelatePlotReprof(Plot):
    PLOT_NAME = "Helix (H) and Sheet (E) Content Per Protein"
    SOURCE_TYPE = "reprof pbased"
    FILE_NAME = "reprof_p_content_relate"

    jp: sns.JointGrid

    def _run(self, df: pd.DataFrame):
        ax1 = plt.subplot()
        self.jp = sns.jointplot(data=df, x="E", y="H",
                           hue="proteome",
                           kind="kde",
                           joint_kws={'alpha': 0.2,
                                      'thresh': .2,     # 'shade_lowest': True == 'thresh': 0
                                      'levels': [0.05, 1.0], # standard: [0.2, 1.0]
                                      'bw_adjust': .8,
                                      'fill': True,
                                      'kwargs': {'linewidth': 2,},
                                      },
                           xlim=(-0.05, 1.05),
                           ylim=(-0.05, 1.05),
                           marginal_kws={
                               "common_norm": False,
                               "bw_adjust": .5,
                               "fill": True,
                           },
                           # ax=ax1,
                           palette=self.get_color_scheme(),
                           hue_order=self.get_proteome_names(),
                           label="xkcd",
                           )
        self.jp.plot_joint(sns.kdeplot, hue="proteome", zorder=10, fill=False, bw_adjust=0.8, levels=[0.05, 1.0], thresh=0.2)

        self.rename_legend(ax1)
        self.jp.fig.subplots_adjust(top=0.95)  # Reduce plot to make room

    def set_title(self):
        self.jp.fig.suptitle(self.PLOT_NAME, fontsize=12, y=0.97)

    def rename_legend(self, ax):
        handles = self.jp.ax_joint.legend_.legendHandles
        self.jp.ax_joint.legend_._visible = False
        self.jp.fig.legend(
            handles=handles,
            labels=self.get_proteome_names().values(),
            loc=1,
            bbox_to_anchor=(0.82, 0.75),
            frameon=True,
            fancybox=True,
            framealpha=0.8,
            facecolor="white",
        )
        leg = ax.legend()


