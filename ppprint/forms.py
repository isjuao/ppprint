from django import forms
from ppprint.models import VisualizationJob, ImportJob, StatusChoices
from ppprint.validators import validate_color


class UploadForm(forms.Form):
    name = forms.CharField(max_length=200)
    file = forms.FileField()
    color = forms.CharField(
        max_length=7,
        validators=[validate_color],
        widget=forms.TextInput(attrs={"type": "color"}),
        initial="#ffffff",
    )

    def clean_color(self):
        color = self.cleaned_data["color"]
        return "" if color == "#ffffff" else color


class SelectionForm(forms.ModelForm):
    sources = forms.ModelMultipleChoiceField(
        queryset=ImportJob.objects.filter(status=StatusChoices.SUCCESS)
    )

    class Meta:
        model = VisualizationJob
        fields = ["sources"]
