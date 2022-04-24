import datetime
from typing import Optional, Tuple

from django.db import models
from django.utils.translation import gettext_lazy as _

from ppprint.validators import validate_color


class StatusChoices(models.TextChoices):
    CREATED = "CREATED", _("Created")
    RUNNING = "RUNNING", _("Running")
    SUCCESS = "SUCCESS", _("Success")
    FAILURE = "FAILURE", _("Failure")


class Job(models.Model):
    status = models.CharField(
        max_length=7, choices=StatusChoices.choices, default=StatusChoices.CREATED
    )
    created_at = models.DateTimeField(
        auto_now_add=True,
    )

    class Meta:
        abstract = True


class ImportJob(Job):
    name = models.CharField(max_length=200, blank=False)
    color = models.CharField(
        max_length=7, blank=True, default="", validators=[validate_color]
    )

    def __str__(self):
        return f"Proteome: {self.name}"

    def get_rgb_parts(self) -> Optional[Tuple[float, float, float]]:
        """Returns RGB tuple mapped between 0 and 1 for a user-selected hex color, else None."""

        get_part = lambda x: int(x, 16) / 255
        if color := self.color:
            return get_part(color[1:3]), get_part(color[3:5]), get_part(color[5:7])
        else:
            return None


class VisualizationJob(Job):
    sources = models.ManyToManyField("ImportJob")
