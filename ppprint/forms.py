from django import forms
from ppprint.models import VisualizationJob, ImportJob, StatusChoices


class UploadForm(forms.Form):
    name = forms.CharField(max_length=200)
    file = forms.FileField()


class SelectionForm(forms.ModelForm):
    sources = forms.ModelMultipleChoiceField(queryset=ImportJob.objects.filter(status=StatusChoices.SUCCESS))

    class Meta:
        model = VisualizationJob
        fields = ["sources"]
