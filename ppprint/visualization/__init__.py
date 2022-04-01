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

MDISORDER = [
    RLengthDistributionPlotRelMdisorder,
    RLengthDistributionPlotAbsMdisorder,
    RPointLinePlotMdisorder,
    RSpectrumPlotMdisorder,
    PContentPerProteomePlotMdisorder,
    PCompositionPerProteomePlotMdisorder,
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
]
PRONA = [
    RLengthDistributionPlotRelProna,
    RLengthDistributionPlotAbsProna,
    RPointLinePlotProna,
    PBindingElementsPlotProna,
    PContentPerProteomePlotProna,
    PProtClassFractionsProna,
]
COMBINED = [
    PLengthDistributionPlot,
]

PLOTS = MDISORDER + TMSEG + PRONA + COMBINED
