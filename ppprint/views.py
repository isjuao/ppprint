from pathlib import Path

from django.conf import settings
from django.core.files.storage import default_storage
from django.http import Http404, JsonResponse, HttpResponse
from django.shortcuts import redirect, render, reverse

from ppprint.forms import SelectionForm, UploadForm
from ppprint.models import ImportJob, VisualizationJob, StatusChoices
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
            return redirect("import_job_status_page", pk=import_job.pk)
    else:
        form = UploadForm()
    return render(request, "ppprint/upload.html", {"form": form})


def compare_proteomes(request):
    if request.method == "POST":
        form = SelectionForm(request.POST)
        if form.is_valid():
            visualization_job = form.save()
            run_visualization_job.delay(visualization_job.pk)
            return redirect("visualization_job_status_page", pk=visualization_job.pk)
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
    if view != "reprof":
        try:
            plot_classes = plots[view]
        except KeyError:
            raise Http404("Feature view does not exist.")
        vj = VisualizationJob.objects.get(pk=pk)
        base_path = settings.MEDIA_URL + f"visualization_job/{pk}/"
        mapping = [
            (base_path + cls.FILE_NAME + ".png", cls.PLOT_NAME, i)
            for i, cls in enumerate(plot_classes)
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
    # Testing purposes
    else:
        vj = VisualizationJob.objects.get(pk=pk)
        base_path = settings.MEDIA_URL + f"visualization_job/{pk}/"
        mapping_dict = {
            cls.FILE_NAME: (base_path + cls.FILE_NAME + ".png", cls.PLOT_NAME, i)
            for i, cls in enumerate(plots[view])
        }
        mapping = [
            (base_path + cls.FILE_NAME + ".png", cls.PLOT_NAME, i)
            for i, cls in enumerate(plots[view])
        ]
        return render(
            request,
            "ppprint/plots_reprof.html",
            {
                "job": vj,
                "view": view,
                "mapping": mapping_dict,
                "basepath": base_path,
            },
        )


def import_job_status_page(request, pk):
    ij = ImportJob.objects.get(pk=pk)
    if ij.status == StatusChoices.SUCCESS or ij.status == StatusChoices.FAILURE:
        return render(request, "ppprint/import_finished.html", {"job": ij})
    else:
        return render(request, "ppprint/import_load.html")


def direct_visualization(request, pk):
    ij = ImportJob.objects.filter(pk=pk)
    vj = VisualizationJob.objects.create()
    vj.sources.set(ij)
    run_visualization_job.delay(vj.pk)
    return redirect("visualization_job_status_page", pk=vj.pk)


def visualization_job_status_page(request, pk):
    vj = VisualizationJob.objects.get(pk=pk)
    if vj.status == StatusChoices.SUCCESS:
        return redirect(reverse("visualization", kwargs={"pk": pk}))
    elif vj.status == StatusChoices.FAILURE:
        return redirect("list_visualization_jobs", pk=pk)
    else:
        return render(request, "ppprint/vis_load.html")
