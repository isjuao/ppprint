from django.conf import settings
from django.core.files.storage import default_storage
from django.http import Http404
from django.shortcuts import redirect, render, reverse

from ppprint.forms import SelectionForm, UploadForm
from ppprint.models import ImportJob, VisualizationJob
from ppprint.tasks import run_import_job, run_visualization_job
from ppprint.visualization import ALL, MDISORDER, PLOTS, PRONA, TMSEG


def home(request):
    return render(request, "ppprint/home.html")


def create_import_job(request):
    if request.method == "POST":
        form = UploadForm(request.POST, request.FILES)
        if form.is_valid():
            import_job = ImportJob.objects.create(
                name=form.cleaned_data["name"], color=form.cleaned_data["color"]
            )
            filename = request.FILES["file"].name
            default_storage.save(
                f"import_job/{import_job.pk}/{filename}", request.FILES["file"]
            )
            run_import_job.delay(import_job.pk)
            return redirect(reverse("home"))
    else:
        form = UploadForm()
    return render(request, "ppprint/upload.html", {"form": form})


def compare_proteomes(request):
    if request.method == "POST":
        form = SelectionForm(request.POST)
        if form.is_valid():
            visualization_job = form.save()
            run_visualization_job.delay(visualization_job.pk)
            return redirect(reverse("home"))
    else:
        form = SelectionForm()
    return render(request, "ppprint/selection.html", {"form": form})


def list_visualization_jobs(request):
    jobs = VisualizationJob.objects.all()
    return render(request, "ppprint/list.html", {"jobs": jobs})


def detail_visualization_job(request, pk):
    view = request.GET.get("view", "overview")
    plots = {"overview": ALL, "mdisorder": MDISORDER, "tmseg": TMSEG, "prona": PRONA}
    try:
        plot_classes = plots[view]
    except KeyError:
        raise Http404("Feature view does not exist.")
    vj = VisualizationJob.objects.get(pk=pk)
    base_path = settings.MEDIA_URL + f"visualization_job/{pk}/"
    mapping = [
        (base_path + cls.FILE_NAME + ".png", cls.PLOT_NAME) for cls in plot_classes
    ]
    return render(request, "ppprint/plots.html", {"job": vj, "mapping": mapping, "view": view,})
