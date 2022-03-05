from typing import Type
from django.conf import settings
from pathlib import Path

from ppprint.celery import app
from ppprint.models import ImportJob, Job, StatusChoices, VisualizationJob
from ppprint.preprocessing.run import run_extract, run_info
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
    results = run_info(import_job_pk, json_path)


@app.task(bind=True, name="run_visualization_job")
@watchdog(VisualizationJob)
def run_visualization_job(self, visualization_job_pk: int):
    # FIXME load data from pickle

    job = VisualizationJob.objects.get(pk=visualization_job_pk)
    results = {}
    for source in job.sources.all():    # sources are ImportJobs
        json_path = (
                Path(settings.BASE_DIR)
                / settings.MEDIA_ROOT
                / "import_job"
                / str(source.pk)
                / "data.json"
        )
        results[source.pk] = run_info(json_path)

    run(visualization_job_pk, results)

