from typing import Type

from ppprint.celery import app
from ppprint.models import ImportJob, Job
from ppprint.preprocessing.run import run_extract, run_info


def watchdog(cls: Type[Job]):
    def outer(f):
        def inner(self, pk, *args, **kwargs):
            cls.objects.filter(pk=pk).update(status=Job.StatusChoices.RUNNING)

            try:
                result = f(self, pk, *args, **kwargs)
            except Exception as exc:
                cls.objects.filter(pk=pk).update(status=Job.StatusChoices.FAILURE)
                raise exc

            cls.objects.filter(pk=pk).update(status=Job.StatusChoices.SUCCESS)
            return result

        return inner

    return outer


@app.task(bind=True, name="run_import_job")
@watchdog(ImportJob)
def run_import_job(self, import_job_pk: int):
    json_path = run_extract(import_job_pk)
    results = run_info(import_job_pk, json_path)
