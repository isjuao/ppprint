# Generated by Django 4.0.2 on 2022-03-16 11:15

from django.db import migrations, models

import ppprint.validators


class Migration(migrations.Migration):

    dependencies = [
        ("ppprint", "0003_remove_importjob_proteome_importjob_name_and_more"),
    ]

    operations = [
        migrations.AddField(
            model_name="importjob",
            name="color",
            field=models.CharField(
                blank=True,
                default="",
                max_length=7,
                validators=[ppprint.validators.validate_color],
            ),
        ),
    ]
