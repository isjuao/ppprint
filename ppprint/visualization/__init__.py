"""
Contains all `Plot` subclasses to be displayed in ppprint.
"""

from ppprint.visualization.plot_content_proteome import (
    PCompositionPerProteomePlotMdisorder,
    PContentPerProteomePlotMdisorder,
    PContentPerProteomePlotProna,
    PContentPerProteomePlotTmseg,
)
from ppprint.visualization.plot_elements_heatmap import (
    PBindingElementsPlotProna,
    PSecStrElementsPlotReprof,
)
from ppprint.visualization.plot_fractions import (
    PProtClassFractionsProna,
    PResidueFractionsBarsTmseg,
    PResidueFractionsViolinsTmseg,
    PResidueFractionsReprof,
)
from ppprint.visualization.plot_length_distribution import (
    PLengthDistributionPlot,
    RLengthDistributionPlotAbsMdisorder,
    RLengthDistributionPlotAbsProna,
    RLengthDistributionPlotAbsTmseg,
    RLengthDistributionPlotRelMdisorder,
    RLengthDistributionPlotRelProna,
    RLengthDistributionPlotRelTmseg,
    RLengthsPlotAbsMdisorderProna,
)
from ppprint.visualization.plot_pie_charts import (
    POrientationsPlotTmseg,
    PProtClassPlotTmseg,
)
from ppprint.visualization.plot_points import (
    RPointLinePlotMdisorder,
    RPointLinePlotProna,
    RPointLinePlotTmseg,
    RPointLinePlotReprof,
)
from ppprint.visualization.plot_spectrum import RSpectrumPlotMdisorder
from ppprint.visualization.plot_proteome_sizes import PProteomeSizes
from ppprint.visualization.plot_num_regions import (
    PNumberOfRegionsMdisorder,
    PNumberOfRegTopoTmseg,
    PNumberOfRegionsTmseg,
    PNumberOfRegionsProna,
)
from ppprint.visualization.plot_content_protein import (
    PContentPerProteinPlotMdisorder,
    PContentPerProteinPlotTmseg,
    PContentPerProteinPlotProna,
)
from ppprint.visualization.plot_pbr_per_dr import RRegionPlotMdisorderProna
from ppprint.visualization.plot_venn_overlap import POverlapPlotTmsegMdisorder
from ppprint.visualization.plot_scatter_mixed import RScatterPlotMdisorderProna
from ppprint.visualization.plot_content_relate import PContentRelatePlotReprof


MDISORDER = [
    RLengthDistributionPlotRelMdisorder,
    RLengthDistributionPlotAbsMdisorder,
    RPointLinePlotMdisorder,
    RSpectrumPlotMdisorder,  # first set max num choices validator in ppprint/forms.py!
    PContentPerProteomePlotMdisorder,
    PCompositionPerProteomePlotMdisorder,
    PNumberOfRegionsMdisorder,
    PContentPerProteinPlotMdisorder,
]
TMSEG = [
    RLengthDistributionPlotRelTmseg,
    RLengthDistributionPlotAbsTmseg,
    RPointLinePlotTmseg,
    POrientationsPlotTmseg,
    PProtClassPlotTmseg,
    PContentPerProteomePlotTmseg,
    PResidueFractionsBarsTmseg,
    PResidueFractionsViolinsTmseg,
    PNumberOfRegionsTmseg,
    PNumberOfRegTopoTmseg,
    PContentPerProteinPlotTmseg,
]
PRONA = [
    RLengthDistributionPlotRelProna,
    RLengthDistributionPlotAbsProna,
    RPointLinePlotProna,
    PBindingElementsPlotProna,
    PContentPerProteomePlotProna,
    PProtClassFractionsProna,
    PNumberOfRegionsProna,
    PContentPerProteinPlotProna,
]
REPROF = [
    RPointLinePlotReprof,
    PSecStrElementsPlotReprof,
    PResidueFractionsReprof,
    PContentRelatePlotReprof,
]
ALL = [
    PLengthDistributionPlot,
    PProteomeSizes,
]
COMBINED = [
    POverlapPlotTmsegMdisorder,
    RLengthsPlotAbsMdisorderProna,
    RRegionPlotMdisorderProna,
    RScatterPlotMdisorderProna,
]

PLOTS = MDISORDER + TMSEG + PRONA + REPROF + COMBINED + ALL
