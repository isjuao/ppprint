from django.conf import settings
from django.conf.urls.static import static
from django.contrib import admin
from django.urls import path

from django.contrib.staticfiles.urls import staticfiles_urlpatterns

from ppprint.views import (
    compare_proteomes,
    create_import_job,
    detail_visualization_job,
    home,
    list_visualization_jobs,
    loading_screen,
    direct_visualization,
)

urlpatterns = [
    path("admin/", admin.site.urls),
    path("upload", create_import_job, name="create_import_job"),
    path("", home, name="home"),
    path("select", compare_proteomes, name="compare_proteomes"),
    path("list", list_visualization_jobs, name="list_visualization_jobs"),
    path("visualization/<int:pk>", detail_visualization_job, name="visualization"),
    path("load/<int:pk>", loading_screen, name="loading_screen"),
    path("finished/<int:pk>", direct_visualization, name="direct_visualization")
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)

urlpatterns += staticfiles_urlpatterns()
