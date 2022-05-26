from crispy_forms.helper import FormHelper
from django import forms

from ppprint.models import ImportJob, StatusChoices, VisualizationJob
from ppprint.validators import limit_num_choices, validate_color
from django.forms import ModelMultipleChoiceField


class UploadForm(forms.Form):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.form_tag = False

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
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.form_tag = False

    class MyModelMultipleChoiceField(ModelMultipleChoiceField):
        def label_from_instance(self, obj):
            return f"({obj.pk}) Proteome {obj.name}"

    sources = MyModelMultipleChoiceField(
        queryset=ImportJob.objects.filter(status=StatusChoices.SUCCESS),
        validators=[limit_num_choices(4)],
    )

    class Meta:
        model = VisualizationJob
        fields = ["sources"]
