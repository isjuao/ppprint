from django.core.files.storage import default_storage
from django.shortcuts import redirect, render, reverse

from ppprint.forms import UploadForm
from ppprint.models import ImportJob, Proteome
from ppprint.tasks import run_import_job


def home(request):
    return render(request, "ppprint/home.html")


def create_import_job(request):
    if request.method == "POST":
        form = UploadForm(request.POST, request.FILES)
        if form.is_valid():
            proteome = Proteome.objects.create(name=form.cleaned_data["name"])
            import_job = ImportJob.objects.create(proteome=proteome)
            filename = request.FILES["file"].name
            default_storage.save(
                f"import_job/{import_job.pk}/{filename}", request.FILES["file"]
            )
            run_import_job.delay(import_job.pk)
            return redirect(reverse("home"))
    else:
        form = UploadForm()
    return render(request, "ppprint/upload.html", {"form": form})
