"""
Runs preprocessing of raw upload data.
"""

import os
import pickle
import tarfile
from pathlib import Path
from typing import Dict

import pandas as pd

from django.conf import settings

from ppprint.models import ImportJob
from ppprint.preprocessing.parse import write_json
from ppprint.preprocessing.utils import LoggedException
from ppprint.preprocessing.extract import extract_pbased, extract_rbased, read_json


def extract_data(base_folder: Path, data_folder: Path):
    """Unpacks .tar and .tar.gz files into job folders."""
    archive = next(base_folder.iterdir())

    # If .pickle already present, but need new DataFrames
    for item in base_folder.iterdir():
        if not item.is_dir() and tarfile.is_tarfile(item):
            archive = item

    try:
        with tarfile.open(archive, "r") as tf:
            tf.extractall(data_folder)
    except tarfile.ReadError:
        message = "Failed to read proteome file. No ImportJob was created!"
        raise LoggedException(message)


def run_extract(import_job_pk: int):
    """Extracts uploaded archives and preprocesses data into JSON format for a given proteome."""

    base_folder = get_base_folder(import_job_pk)
    data_folder = base_folder / "data"
    json_path = base_folder / "data.json"

    extract_data(base_folder, data_folder)
    write_json(data_folder, json_path, import_job_pk)

    return json_path


def get_base_folder(import_job_pk: int):
    base_folder = (
        Path(settings.BASE_DIR)
        / settings.MEDIA_ROOT
        / "import_job"
        / str(import_job_pk)
    )
    return base_folder


def run_info(json_path: Path) -> Dict[str, pd.DataFrame]:
    """Preprocesses data from JSON into info-containing data frames for a given proteome."""

    df_source, df_seq = read_json(json_path)

    results = {}
    # Will contain both pbased and rbased results for each feature, for the given proteome

    # Protein-level extraction
    results.update(extract_pbased(df_source, df_seq))

    # Region-level extraction
    results.update(extract_rbased(df_source, df_seq))

    return results


def store(results: Dict[str, pd.DataFrame], path: Path):
    with open(path, "wb") as f:
        pickle.dump(results, f)


def load(path: Path) -> Dict[str, pd.DataFrame]:
    with open(path, "rb") as f:
        return pickle.load(f)


# if __name__ == "__main__":
#     # For local testing only
#     # Comment out django imports first
#     example_json_path = Path(
#         "/home/isabell/work/python/thesis/ppprint/tests/data/oneprotein/data.json"
#     )
#     run_info(example_json_path)
