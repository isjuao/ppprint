import shutil
import os
from http import HTTPStatus
from unittest.mock import patch
import json
import pytest
from pathlib import Path
from django.conf import settings

from ppprint.preprocessing.run import extract_data, write_json, LoggedException
from ppprint.models import ImportJob, StatusChoices
from ppprint.tasks import run_import_job
from tests.steps.utils import build_true_segments_json, convert_mdisorder_to_latin1


def test_extract_tarfile():
    """
    Tests whether ppprint extracts job folders correctly upon upload of a correct archive
    with missing job folders (job_2), which should not influence ppprint behavior.
    """

    base_folder = Path(settings.BASE_DIR) / "tests" / "data" / "sarscov2"
    data_folder = base_folder / "data"

    extract_data(base_folder, data_folder)
    assert len(os.listdir(data_folder)) == 16

    # Delete data folder for next test
    shutil.rmtree(data_folder)


def test_not_tarfile():
    """Tests whether ppprint raises the correct exception upon upload of a non-archive type file."""

    base_folder = Path(settings.BASE_DIR) / "tests" / "data" / "novalidtar"

    with pytest.raises(LoggedException):
        extract_data(base_folder, (base_folder / "data"))


@pytest.mark.django_db()
def test_no_jobfolder(import_job_factory):
    """
    Tests whether ppprint raises the correct exception upon upload of an archive that does contain
    .fasta files, but a false structure (no job folders).
    """

    # Create ImportJob
    ij = import_job_factory(Path(settings.BASE_DIR) / "tests" / "data" / "nojobfolder")

    base_folder = Path(settings.BASE_DIR) / settings.MEDIA_ROOT / "import_job" / str(ij.pk)
    data_folder = base_folder / "data"
    json_path = base_folder / "data.json"

    extract_data(base_folder, data_folder)
    assert len(os.listdir(data_folder)) > 0

    with pytest.raises(LoggedException):
        write_json(data_folder, json_path, ij.pk)


@pytest.mark.django_db()
def test_import(client):
    """
    Tests whether ppprint extracts upload data into intermediate .pickle files
    and correctly imports upload data into the database.
    """

    # Patch celery method to run it explicitly instead of queueing task
    with patch("ppprint.tasks.run_import_job.delay") as mock_method:
        data_path = Path(settings.BASE_DIR) / "tests" / "data" / "sarscov2" / "sarscov2.tar.gz"
        with open(data_path, "rb") as f:
            response = client.post("/upload", {"name": "example", "file": f, "color": "#000000"})
        # Assert working redirect to loading screen
        assert response.status_code == HTTPStatus.FOUND
        mock_method.assert_called_once()

    pk = ImportJob.objects.first().pk
    run_import_job(pk)

    assert ImportJob.objects.get(pk=pk).status == StatusChoices.SUCCESS
    assert (Path(settings.BASE_DIR) / settings.MEDIA_ROOT / "import_job" / str(pk) / "results.pickle").exists()


@pytest.mark.django_db()
def test_failed_import(client):
    """
    Tests whether ppprint does not throw exception but stores error message when the user is at fault,
    e.g. uploading a file that is not a tar archive.
    """

    # Patch celery method to run it explicitly instead of queueing task
    with patch("ppprint.tasks.run_import_job.delay") as mock_method:
        data_path = Path(settings.BASE_DIR) / "tests" / "data" / "novalidtar" / "thisisnotatarfile.txt"
        with open(data_path, "rb") as f:
            response = client.post("/upload", {"name": "badexample", "file": f, "color": "#000000"})
        mock_method.assert_called_once()

    ij = ImportJob.objects.first()

    assert not ij.messages.exists()
    run_import_job(ij.pk)

    # No exception should be thrown, but there should be an error message in the database
    assert ij.messages.exists()


@pytest.mark.django_db()
def test_parse_segments(import_job_factory):
    """
    Tests whether ppprint correctly parses the region segments from PredictProtein output files for two E. coli
    proteins. Ground truth data is built manually.
    """

    # Create ImportJob
    ij = import_job_factory(Path(settings.BASE_DIR) / "tests" / "data" / "tinyproteome")

    base_folder = Path(settings.BASE_DIR) / settings.MEDIA_ROOT / "import_job" / str(ij.pk)
    data_folder = base_folder / "data"
    json_path = base_folder / "data.json"

    # Get test data
    extract_data(base_folder, data_folder)
    write_json(data_folder, json_path, ij.pk)

    # Get ground truth data
    truth_data = build_true_segments_json()

    with open(json_path, "r") as f:
        test_data = json.load(f)

    assert truth_data == test_data


@pytest.mark.django_db()
def test_corrupt_files(client):
    """
    Tests whether ppprint extracts upload data into intermediate .pickle files
    and correctly imports upload data into the database, while dealing with
    - missing (no .reprof file for P62524) or
    - corrupt (missing column in .mdisorder file for Q8XA85)
    PredictProtein output files.
    Warning messages should be written into the database entry and made available to the user,
    while still marking the ImportJob as SUCCESS.
    """

    # Patch celery method to run it explicitly instead of queueing task
    with patch("ppprint.tasks.run_import_job.delay") as mock_method:
        data_path = Path(settings.BASE_DIR) / "tests" / "data" / "tinyproteome_miss" / "tinyproteome_miss.tar.gz"
        with open(data_path, "rb") as f:
            response = client.post("/upload", {"name": "anotherexample", "file": f, "color": "#000000"})
        # Assert working redirect to loading screen
        assert response.status_code == HTTPStatus.FOUND
        mock_method.assert_called_once()

    ij = ImportJob.objects.first()
    run_import_job(ij.pk)
    ij.refresh_from_db()

    # No exception should be thrown, but there should be warnings in the database
    messages = ij.messages.all()
    assert len(messages) == 2
    assert messages[0].text == "Could not FIND P62524.reprof in job_1."
    assert messages[1].text == "Could not PARSE Q8XA85.mdisorder in job_1."
    # Parsing result file should still exist
    assert (Path(settings.BASE_DIR) / settings.MEDIA_ROOT / "import_job" / str(ij.pk) / "results.pickle").exists()
    # Job status should not be failure
    assert ij.status == StatusChoices.SUCCESS


@pytest.mark.django_db()
def test_latin1_files(client, import_job_factory):
    """
    Tests whether ppprint extracts and correctly imports upload data,
    for latin-1-encoded .mdisorder files.
    """
    # Create ImportJob
    ij = import_job_factory(Path(settings.BASE_DIR) / "tests" / "data" / "tinyproteome_latin-1")

    base_folder = Path(settings.BASE_DIR) / settings.MEDIA_ROOT / "import_job" / str(ij.pk)
    data_folder = base_folder / "data"
    json_path = base_folder / "data.json"

    extract_data(base_folder, data_folder)

    # Convert .mdisorder files to latin-1 encoding
    convert_mdisorder_to_latin1(data_folder)

    write_json(data_folder, json_path, ij.pk)

    # Assert job did not throw error messages
    assert ij.messages.count() == 0


