import shutil
from pathlib import Path
import pytest
from ppprint.models import ImportJob, StatusChoices


@pytest.fixture(autouse=True)
def custom_media_root(settings, request):
    """Fixture to set media folder settings for tests and cleanup extracted media after."""

    settings.MEDIA_ROOT = "tmp/test_media"

    def delete():
        try:
            shutil.rmtree(settings.MEDIA_ROOT)
        except Exception:
            pass

    request.addfinalizer(delete)


@pytest.fixture
def import_job_factory(settings):
    """Fixture to create a new ImportJob."""

    def get_job(test_data_path: Path):
        """Creates new ImportJob and copies files to respective folder in media."""

        ij = ImportJob.objects.create(name="example", status=StatusChoices.RUNNING)
        dst_path = Path(settings.BASE_DIR) / settings.MEDIA_ROOT / "import_job" / str(ij.pk)
        shutil.copytree(src=test_data_path, dst=dst_path)

        return ij

    return get_job