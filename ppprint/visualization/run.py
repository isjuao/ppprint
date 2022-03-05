from typing import Dict
import pandas as pd
from collections import defaultdict
from ppprint.visualization.plot import LengthDistributionPlot
from ppprint.models import VisualizationJob
from django.conf import settings
from pathlib import Path


def run(visualization_job_pk: int, data: Dict[int, Dict[str, pd.DataFrame]]):
    result = defaultdict(pd.DataFrame)
    vj = VisualizationJob.objects.get(pk=visualization_job_pk)

    mapping = {}
    for pk, name in vj.sources.all().values_list("pk", "name"):
        mapping[pk] = name

    # Add proteome IDs (pks) to dataframe and build big source type df with all requested proteomes
    for proteome_id, dataset in data.items():
        for source_type, df in dataset.items():
            df["proteome"] = proteome_id
            result[source_type] = pd.concat([result[source_type], df], ignore_index=True)

    base_folder = (
        Path(settings.BASE_DIR)
        / settings.MEDIA_ROOT
        / "visualization_job"
        / str(visualization_job_pk)
    )
    base_folder.mkdir(exist_ok=True, parents=True)

    run_plotting(dict(result), mapping, base_folder)


def run_plotting(dataframes: Dict[str, pd.DataFrame], mapping: Dict[int, str], base_folder: Path):
    PLOTS = [LengthDistributionPlot,]
    for plot_cls in PLOTS:
        plot_cls(dataframes, mapping, base_folder).run()