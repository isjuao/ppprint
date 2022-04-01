from django import forms

from ppprint.models import ImportJob, StatusChoices, VisualizationJob
from ppprint.validators import limit_num_choices, validate_color


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
        queryset=ImportJob.objects.filter(status=StatusChoices.SUCCESS),
        validators=[limit_num_choices(4)],
    )

    class Meta:
        model = VisualizationJob
        fields = ["sources"]
