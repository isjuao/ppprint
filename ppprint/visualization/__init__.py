from ppprint.visualization.length_distribution import (
    PLengthDistributionPlot,
    RLengthDistributionPlotRelMdisorder,
    RLengthDistributionPlotRelTmseg,
    RLengthDistributionPlotRelProna,
    RLengthDistributionPlotAbsMdisorder,
    RLengthDistributionPlotAbsTmseg,
    RLengthDistributionPlotAbsProna,
)


MDISORDER = [
    RLengthDistributionPlotRelMdisorder,
    RLengthDistributionPlotAbsMdisorder,
]
TMSEG = [
    RLengthDistributionPlotRelTmseg,
    RLengthDistributionPlotAbsTmseg,
]
PRONA = [
    RLengthDistributionPlotRelProna,
    RLengthDistributionPlotAbsProna,
]
COMBINED = [
    PLengthDistributionPlot,
]

PLOTS = MDISORDER + TMSEG + PRONA + COMBINED
