from django.db import models
from django.utils.translation import gettext_lazy as _


class Job(models.Model):
    class StatusChoices(models.TextChoices):
        CREATED = "CREATED", _("Created")
        RUNNING = "RUNNING", _("Running")
        SUCCESS = "SUCCESS", _("Success")
        FAILURE = "FAILURE", _("Failure")

    status = models.CharField(
        max_length=7, choices=StatusChoices.choices, default=StatusChoices.CREATED
    )

    class Meta:
        abstract = True


class ImportJob(Job):
    proteome = models.ForeignKey("Proteome", on_delete=models.CASCADE)


class DataSource(models.Model):
    proteome = models.ForeignKey("Proteome", on_delete=models.CASCADE)


class VisualizationJob(Job):
    sources = models.ManyToManyField("DataSource")


class Proteome(models.Model):
    name = models.CharField(max_length=200, blank=False)
    uniprot_id = models.CharField(max_length=11, blank=True, default="")
    description = models.TextField(blank=True, default="")
