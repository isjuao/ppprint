import tarfile
from pathlib import Path

import pandas as pd
from django.conf import settings

from ppprint.models import ImportJob
from ppprint.preprocessing.parse import write_json
from ppprint.preprocessing.extract import read_json, extract_pbased, extract_rbased


def extract_data(base_folder: Path, data_folder: Path):
    """Unpacks .tar and .tar.gz files into job folders."""
    archive = next(base_folder.iterdir())

    with tarfile.open(archive, "r") as tf:
        tf.extractall(data_folder)


def run_extract(import_job_pk: int):
    """Extracts uploaded archives and preprocesses data into JSON format for a given proteome."""

    base_folder = (
        Path(settings.BASE_DIR)
        / settings.MEDIA_ROOT
        / "import_job"
        / str(import_job_pk)
    )
    data_folder = base_folder / "data"
    json_path = base_folder / "data.json"

    extract_data(base_folder, data_folder)
    write_json(data_folder, json_path)

    return json_path


def run_info(json_path: Path):
    """Preprocesses data from JSON into info-containing data frames for a given proteome."""

    df_source, df_seq = read_json(json_path)

    results = {}
    # will contain both pbased and rbased results for each feature, for the given proteome

    # Protein-level extraction
    # results.update(extract_pbased(df_source, df_seq))

    # Region-level extraction
    results.update(extract_rbased(df_source, df_seq))

    return results


if __name__ == '__main__':
    # comment out django imports first
    example_json_path = Path("/home/isabell/work/python/thesis/ppprint/media/import_job/5/data.json")
    run_info(example_json_path)