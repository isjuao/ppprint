from pathlib import Path
from typing import Dict
import pandas as pd
from django.conf import settings

from ppprint.preprocessing.run import run_info
from tests.steps.utils import build_true_dfs_pbased, build_true_dfs_rbased


def test_extract():
    """Tests whether ppprint extracts the correct protein-based and region-based info from a data.JSON file."""

    json_path = (
        Path(settings.BASE_DIR) / "tests" / "data" / "artificialproteome" / "data.json"
    )
    test_results = run_info(json_path)

    def compare_sources():
        for source in true_results.keys():
            test_df = test_results[source].reset_index(drop=True)
            true_df = true_results[source].reset_index(drop=True)

            # Assert equality of dtypes
            assert test_df.dtypes.equals(true_df.dtypes)
            # Assert equality of contents
            message = f"{source} is wrong!"
            pd.testing.assert_frame_equal(test_df, true_df, obj=message)

    # Protein-based
    true_results = build_true_dfs_pbased()
    compare_sources()

    # Region-based
    true_results = build_true_dfs_rbased()
    compare_sources()
