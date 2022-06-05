"""
Provides functionalities for testing in local Django shell.
For testing, import respective `Plot` subclass and `test_plots`.
Call `test_plots` with desired `VisualizationJob` pk and plot class.
"""

from pathlib import Path

import pandas as pd
from django.conf import settings

from ppprint.models import VisualizationJob
from ppprint.preprocessing.run import run_info
from ppprint.visualization.run import prepare


def test_plots(visualization_job_pk, plot_cls):
    job = VisualizationJob.objects.get(pk=visualization_job_pk)
    results = {}
    for source in job.sources.all():  # sources are ImportJobs
        json_path = (
            Path(settings.BASE_DIR)
            / settings.MEDIA_ROOT
            / "import_job"
            / str(source.pk)
            / "data.json"
        )
        results[source.pk] = run_info(json_path)

    result_dict, mapping, base_folder = prepare(visualization_job_pk, results)

    plot_cls(result_dict, mapping, base_folder).run()
