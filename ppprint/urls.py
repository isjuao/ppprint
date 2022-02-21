from django.contrib import admin
from django.urls import path
from ppprint.views import create_import_job, home

urlpatterns = [
    path("admin/", admin.site.urls),
    path("upload", create_import_job, name="create_import_job"),
    path("", home, name="home")
]
