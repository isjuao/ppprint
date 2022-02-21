from ppprint.celery import app
from typing import Type
from ppprint.models import Job, ImportJob

from ppprint.preprocessing import run

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
    run(import_job_pk)
