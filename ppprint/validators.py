import re
from typing import List

from django.core.exceptions import ValidationError
from django.utils.translation import gettext_lazy as _


def validate_color(value: str):
    if value and not re.match(r"^#(?:[0-9a-fA-F]{3}){1,2}$", value):
        raise ValidationError(
            _("Invalid color given."),
            params={"value": value},
        )


def limit_num_choices(limit: int):
    def inner(value: List[str]):
        if len(value) > limit:
            raise ValidationError(
                _(f"Selected more than {limit} proteomes."),
                params={"value": value},
            )

    return inner
