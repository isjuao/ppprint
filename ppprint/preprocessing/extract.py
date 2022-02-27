from pathlib import Path

import pandas as pd
import numpy as np

import json
from typing import List, Optional, Tuple


def read_json(path: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Reads the proteome-specific JSON with all features into a dataframe with rows corresponding to regions."""

    with open(path, "r") as f:
        data = json.load(f)
    sequences = []

    def get_data():
        """Fast way to build the dataframe on the fly while effectively looking at rows."""

        for i, p in enumerate(data):
            seq = p["sequence"]
            sequences.append(seq)
            for feature in ["tmseg", "mdisorder", "prona"]:
                for region in p[feature]:
                    yield (
                        i,
                        feature,
                        region["begin"],
                        region["end"],
                        region["description"],
                    )

    df = pd.DataFrame(
        get_data(), columns=["protein", "feature", "begin", "end", "description"]
    )

    # Store extracted sequences in a separate dataframe and extract lengths
    df_seq = pd.DataFrame(sequences, columns=["sequence"])
    df_seq["protein length"] = df_seq["sequence"].map(len)

    return df, df_seq


def extract_pbased_mdisorder(
    df_source,
    *args,
    **kwargs,
):
    """Extracts protein-based data for mdisorder."""

    # Filter (all) regions based on minlength requirement
    minlength = kwargs.pop("minlength")
    df_source = df_source[df_source["region length"] >= minlength]

    df_new = extract_pbased_all(df_source, *args, **kwargs)

    # Calculate region content as last step, when (correct) lengths of all proteins got collected
    df_new["region content"] = df_new["sum region lengths"] / df_new["protein length"]

    return df_new


def extract_pbased_tmseg(df_source: pd.DataFrame, *args, **kwargs):
    """Extracts protein-based data for tmseg.
    First, extracts proteins with signal peptides to apply index shifting of all contained regions.
    Then, performs common data extraction and clips protein length and calculates correct region content.
    """

    # Extract signal peptide containing proteins and shift the indices of their regions
    df_sig = df_source[df_source["description"] == "Signal Peptide"].drop_duplicates(
        "protein"
    )[["protein", "region length"]]
    df_sig = df_sig.rename(columns={"region length": "shift size"}).set_index("protein")
    df_source = df_source.join(df_sig, how="left", on="protein").fillna(0)
    df_source["begin"] = df_source["begin"] - df_source["shift size"]
    df_source["end"] = df_source["end"] - df_source["shift size"]

    # Exclude TMH regions that do not fulfill minlength requirement
    minlength = kwargs.pop("minlength")
    df_source = df_source[
        (df_source["description"] != "Transmembrane Helix")
        | (df_source["region length"] >= minlength)
    ]

    df_new = extract_pbased_all(df_source, *args, **kwargs)

    # Clip protein length (if signal peptide present, else clip == 0) and calculate correct region content
    df_new = df_new.join(df_sig, how="left").fillna(0)
    df_new["protein length"] = df_new["protein length"] - df_new["shift size"]
    df_new["region content"] = df_new["sum region lengths"] / df_new["protein length"]

    # Calculate additional region contents per region type
    region_type_counts = (
        df_source.groupby(["protein", "description"])["region length"]
        .sum()
        .unstack(level=-1)
    )
    df_new = df_new.join(region_type_counts, how="left").fillna(0)
    df_new["M"] = df_new.pop("Signal Peptide") + df_new.pop("Transmembrane Helix")
    df_new = df_new.rename(columns={"Cytoplasmic": "I", "Extracellular": "O"})
    df_new["I"] = df_new["I"] / df_new["protein length"]
    df_new["M"] = df_new["M"] / df_new["protein length"]
    df_new["O"] = df_new["O"] / df_new["protein length"]

    def orientation(df):
        """Determines the orientation for each (TM) protein."""

        df = df.reset_index()
        indices = df.index[df["description"] == "Transmembrane Helix"].tolist()
        if not indices:
            # All TMHs did not fulfill minlength requirement
            return ""
        first = indices[0]
        if first == 0:
            # The first region is the helix -> N-terminus of protein is in membrane
            return "Membrane"
        # TODO: check what happens if signal peptide is immediately before first TMH
        return str(df.at[first - 1, "description"])

    # TODO: check if possible that only Cytoplasmic as description (region but no Transmembrane Helix)
    # df_tmp = df_source[df_source["description"] == "Transmembrane Helix"].drop_duplicates("protein")["protein"]
    # df_source = df_source[df_source["protein"].isin(df_tmp)]

    # Compute the orientation for all TMPs
    orientation_series = df_source.groupby("protein").apply(
        func=orientation,
    )
    orientation_series.name = "orientation"
    df_new = df_new.join(orientation_series, how="left").fillna(0)

    return df_new


def extract_pbased_prona(df_source: pd.DataFrame, *args, **kwargs):
    """Extracts protein-based data for prona."""

    # Filter (all) regions based on minlength requirement
    minlength = kwargs.pop("minlength")
    df_source = df_source[df_source["region length"] >= minlength]

    df_new = extract_pbased_all(df_source, *args, **kwargs)

    # Calculate region content as last step, when (correct) lengths of all proteins got collected
    df_new["region content"] = df_new["sum region lengths"] / df_new["protein length"]

    # TODO: calculate additional region contents and num regions per (binding) region type

    return df_new


def extract_pbased_all(
    df_source: pd.DataFrame,
    df_seq: pd.DataFrame,
    region_type: Optional[List[str]] = None,
    additional_aggfuncs={},
):
    """Extracts protein-based information from source data frame."""

    if region_type:
        # Remove unwanted region types
        df_source = df_source[df_source["description"].isin(region_type)]

    # Calculate common information relevant for all features
    df_new = df_source.groupby("protein").agg(
        **{
            "number of regions": pd.NamedAgg(column="region length", aggfunc="size"),
            "median length": pd.NamedAgg(column="region length", aggfunc="median"),
            "sum region lengths": pd.NamedAgg(column="region length", aggfunc="sum"),
            **additional_aggfuncs,
        }
    )
    df_new = df_new.join(df_seq["protein length"], how="right").fillna(0)

    return df_new


def extract_pbased(df_source: pd.DataFrame, df_seq: pd.DataFrame):
    """"""

    # Store function and required params for each feature
    mapping = {
        "tmseg": (
            extract_pbased_tmseg,
            {"minlength": 12, "region_type": ["Transmembrane Helix"]},
        ),
        "mdisorder": (extract_pbased_mdisorder, {"minlength": 30}),
        "prona": (extract_pbased_prona, {"minlength": 6}),
    }
    results = {}

    for feature, (func, kwargs) in mapping.items():
        # Reduce source dataframe to feature and calculate region lengths
        df_curr = df_source[df_source["feature"] == feature]
        df_curr["region length"] = df_curr["end"] - df_curr["begin"] + 1

        # Perform feature-specific extraction and store pbased results
        df_p = func(df_curr, df_seq, **kwargs)
        results[f"{feature} pbased"] = df_p

    return results
