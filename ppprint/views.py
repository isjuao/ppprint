from pathlib import Path

from django.conf import settings
from django.core.files.storage import default_storage
from django.http import Http404, JsonResponse, HttpResponse
from django.shortcuts import redirect, render, reverse

from ppprint.forms import SelectionForm, UploadForm
from ppprint.models import ImportJob, VisualizationJob
from ppprint.tasks import run_import_job, run_visualization_job
from ppprint.visualization import ALL, MDISORDER, PRONA, TMSEG, REPROF, COMBINED


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
            # return redirect(reverse("home"))
            data = {"pk": import_job.pk}
            # return render(request, "ppprint/load.html", data)
            return redirect("loading_screen", pk=import_job.pk)
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
    plots = {
        "overview": ALL,
        "mdisorder": MDISORDER,
        "tmseg": TMSEG,
        "prona": PRONA,
        "reprof": REPROF,
        "combined": COMBINED,
    }
    try:
        plot_classes = plots[view]
    except KeyError:
        raise Http404("Feature view does not exist.")
    vj = VisualizationJob.objects.get(pk=pk)
    base_path = settings.MEDIA_URL + f"visualization_job/{pk}/"
    mapping = [
        (base_path + cls.FILE_NAME + ".png", cls.PLOT_NAME) for cls in plot_classes
    ]
    return render(
        request,
        "ppprint/plots.html",
        {
            "job": vj,
            "mapping": mapping,
            "view": view,
        },
    )


def loading_screen(request, pk):
    ij = ImportJob.objects.get(pk=pk)
    if ij.status == "SUCCESS" or ij.status == "FAILURE":
        return render(request, "ppprint/job_finished.html", {"job": ij})
    else:
        return render(request, "ppprint/load.html")


# def whatever(request):
#     data = {"foo": "hi"}
#     return JsonResponse(data, safe=False)


def direct_visualization(request, pk):
    ij = ImportJob.objects.filter(pk=pk)
    vj = VisualizationJob.objects.create()
    vj.sources.set(ij)
    run_visualization_job.delay(vj.pk)
    return redirect(reverse("list_visualization_jobs"))
