from pathlib import Path
from typing import Type

from django.conf import settings

from ppprint.celery import app
from ppprint.models import ImportJob, Job, StatusChoices, VisualizationJob
from ppprint.preprocessing.run import (
    get_base_folder,
    load,
    run_extract,
    run_info,
    store,
)
from ppprint.visualization.run import run


def watchdog(cls: Type[Job]):
    def outer(f):
        def inner(self, pk, *args, **kwargs):
            cls.objects.filter(pk=pk).update(status=StatusChoices.RUNNING)

            try:
                result = f(self, pk, *args, **kwargs)
            except Exception as exc:
                cls.objects.filter(pk=pk).update(status=StatusChoices.FAILURE)
                raise exc

            cls.objects.filter(pk=pk).update(status=StatusChoices.SUCCESS)
            return result

        return inner

    return outer


@app.task(bind=True, name="run_import_job")
@watchdog(ImportJob)
def run_import_job(self, import_job_pk: int):
    json_path = run_extract(import_job_pk)
    results = run_info(json_path)
    result_file = get_base_folder(import_job_pk) / "results.pickle"
    store(results, result_file)


@app.task(bind=True, name="run_visualization_job")
@watchdog(VisualizationJob)
def run_visualization_job(self, visualization_job_pk: int):
    job = VisualizationJob.objects.get(pk=visualization_job_pk)

    results = {}
    for source in job.sources.all():  # sources are ImportJobs
        pickle_path = get_base_folder(source.pk) / "results.pickle"
        results[source.pk] = load(pickle_path)

    run(visualization_job_pk, results)
