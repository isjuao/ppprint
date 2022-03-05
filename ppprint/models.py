from django.db import models
from django.utils.translation import gettext_lazy as _


class StatusChoices(models.TextChoices):
    CREATED = "CREATED", _("Created")
    RUNNING = "RUNNING", _("Running")
    SUCCESS = "SUCCESS", _("Success")
    FAILURE = "FAILURE", _("Failure")


class Job(models.Model):
    status = models.CharField(
        max_length=7, choices=StatusChoices.choices, default=StatusChoices.CREATED
    )

    class Meta:
        abstract = True


class ImportJob(Job):
    name = models.CharField(max_length=200, blank=False)

    def __str__(self):
        return f"Proteome: {self.name}"


class VisualizationJob(Job):
    sources = models.ManyToManyField("ImportJob")

