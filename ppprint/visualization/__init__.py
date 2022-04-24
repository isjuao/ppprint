"""
Contains all `Plot` subclasses to be displayed in ppprint.
"""

from ppprint.visualization.plot_content_proteome import (
    PCompositionPerProteomePlotMdisorder,
    PContentPerProteomePlotMdisorder,
    PContentPerProteomePlotProna,
    PContentPerProteomePlotTmseg,
)
from ppprint.visualization.plot_elements_heatmap import PBindingElementsPlotProna
from ppprint.visualization.plot_fractions import (
    PProtClassFractionsProna,
    PResidueFractionsBarsTmseg,
    PResidueFractionsViolinsTmseg,
)
from ppprint.visualization.plot_length_distribution import (
    PLengthDistributionPlot,
    RLengthDistributionPlotAbsMdisorder,
    RLengthDistributionPlotAbsProna,
    RLengthDistributionPlotAbsTmseg,
    RLengthDistributionPlotRelMdisorder,
    RLengthDistributionPlotRelProna,
    RLengthDistributionPlotRelTmseg,
)
from ppprint.visualization.plot_pie_charts import (
    POrientationsPlotTmseg,
    PProtClassPlotTmseg,
)
from ppprint.visualization.plot_points import (
    RPointLinePlotMdisorder,
    RPointLinePlotProna,
    RPointLinePlotTmseg,
)
from ppprint.visualization.plot_spectrum import RSpectrumPlotMdisorder
from ppprint.visualization.plot_proteome_sizes import PProteomeSizes
from ppprint.visualization.plot_num_regions import (
    PNumberOfRegionsMdisorder,
    PNumberOfRegTopoTmseg,
    PNumberOfRegionsTmseg,
    PNumberOfRegionsProna,
)


MDISORDER = [
    RLengthDistributionPlotRelMdisorder,
    RLengthDistributionPlotAbsMdisorder,
    RPointLinePlotMdisorder,
    RSpectrumPlotMdisorder,   # first set max num choices validator in ppprint/forms.py!
    PContentPerProteomePlotMdisorder,
    PCompositionPerProteomePlotMdisorder,
    PNumberOfRegionsMdisorder,
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
]
PRONA = [
    RLengthDistributionPlotRelProna,
    RLengthDistributionPlotAbsProna,
    RPointLinePlotProna,
    PBindingElementsPlotProna,
    PContentPerProteomePlotProna,
    PProtClassFractionsProna,
    PNumberOfRegionsProna,
]
ALL = [
    PLengthDistributionPlot,
    PProteomeSizes,
]

PLOTS = MDISORDER + TMSEG + PRONA + ALL
