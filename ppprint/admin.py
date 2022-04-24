from django.contrib import admin

from ppprint.models import ImportJob, VisualizationJob


class ImportJobAdmin(admin.ModelAdmin):
    pass


class VisualizationJobAdmin(admin.ModelAdmin):
    pass


admin.site.register(ImportJob, ImportJobAdmin)
admin.site.register(VisualizationJob, VisualizationJobAdmin)
