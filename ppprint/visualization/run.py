from collections import defaultdict
from itertools import chain
from pathlib import Path
from typing import Dict, Tuple

import pandas as pd
from django.conf import settings

from ppprint.models import VisualizationJob
from ppprint.visualization import PLOTS


def run(visualization_job_pk: int, data: Dict[int, Dict[str, pd.DataFrame]]):
    result_dict, mapping, base_folder = prepare(visualization_job_pk, data)
    run_plotting(result_dict, mapping, base_folder)


def prepare(visualization_job_pk: int, data: Dict[int, Dict[str, pd.DataFrame]]):
    vj = VisualizationJob.objects.get(pk=visualization_job_pk)

    mapping = {}
    colors = {
        "#e377c2": (0.8901960784313725, 0.4666666666666667, 0.7607843137254902),
        "#885fa3": (0.5333333333333333, 0.37254901960784315, 0.6392156862745098),
        "#2ca02c": (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
        "#3c81bd": (0.23529411764705882, 0.5058823529411764, 0.7411764705882353),
        "#8c1c22": (0.5490196078431373, 0.10980392156862745, 0.13333333333333333),
        "#e87833": (0.9098039215686274, 0.47058823529411764, 0.2),
        "#bcbd22": (0.7372549019607844, 0.7411764705882353, 0.13333333333333333),
        "#84e0c6": (0.5176470588235295, 0.8784313725490196, 0.7764705882352941),
        "#102e5c": (0.06274509803921569, 0.1803921568627451, 0.3607843137254902),
    }
    used_colors = set()
    new_color_generator = iter(colors.items())
    duplicated = []

    # Add all proteomes with set colors to mapping
    for source in vj.sources.exclude(color=""):
        if source.color in used_colors:
            # User chose color twice or chose color of selected sample proteome
            duplicated.append(source)
        else:
            used_colors.add(source.color)
            mapping[source.pk] = source.name, source.get_rgb_parts()

    # Find available colors for proteomes without or duplicated colors
    color = next(new_color_generator)
    for source in chain(vj.sources.filter(color=""), duplicated):
        while color[0] in used_colors:
            color = next(new_color_generator)
        mapping[source.pk] = source.name, color[1]
        # Add chosen color to used_colors in order to find a new one for the next proteome
        used_colors.add(color[0])

    result = concat_proteomes(data)

    base_folder = (
        Path(settings.BASE_DIR)
        / settings.MEDIA_ROOT
        / "visualization_job"
        / str(visualization_job_pk)
    )
    base_folder.mkdir(exist_ok=True, parents=True)

    return dict(result), mapping, base_folder
    # run_plotting(dict(result), mapping, base_folder)


def concat_proteomes(data):
    """Adds proteome IDs (pks) to dataframe and builds big source type df with all requested proteomes."""

    result = defaultdict(pd.DataFrame)
    for proteome_id, dataset in data.items():
        for source_type, df in dataset.items():
            df["proteome"] = proteome_id
            result[source_type] = pd.concat(
                [result[source_type], df], ignore_index=True
            )
    return result


def run_plotting(
    dataframes: Dict[str, pd.DataFrame],
    mapping: Dict[int, Tuple[str, Tuple[float, float, float]]],
    base_folder: Path,
):
    for plot_cls in PLOTS:
        plot_cls(dataframes, mapping, base_folder).run()
