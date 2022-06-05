"""
Extracts info from proteome JSON into dataframes
split by feature and base (region/protein).
"""


import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


def read_json(path: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Reads the proteome-specific JSON with all features into a dataframe with rows corresponding to regions."""

    with open(path, "r") as f:
        try:
            data = json.load(f)
        except:
            print(f"Exception occured for {path}")
    sequences = []

    def get_data():
        """Fast way to build the dataframe on the fly while effectively looking at rows."""

        for i, p in enumerate(data):
            seq = p["sequence"]
            sequences.append(seq)
            # TODO: add features here!
            for feature in [
                "tmseg",
                "mdisorder",
                "prona",
                "reprof",
            ]:
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
    df_new = ensure_exists(
        ["Signal Peptide", "Cytoplasmic", "Extracellular", "Transmembrane Helix"],
        df_new,
    )
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
        return str(df.at[first - 1, "description"])

    # Filter to TMPs, because it is possible that regions of proteins without TMH are in df_source, need orientation 0
    df_tmp = df_source[
        df_source["description"] == "Transmembrane Helix"
    ].drop_duplicates("protein")["protein"]
    df_source = df_source[df_source["protein"].isin(df_tmp)]

    # Compute the orientation for all TMPs
    orientation_series = df_source.groupby("protein").apply(
        func=orientation,
    )
    orientation_series.name = "orientation"
    df_new = df_new.join(orientation_series, how="left").fillna(0)

    return df_new


def ensure_exists(columns: List[str], df: pd.DataFrame):
    for column in columns:
        if column not in df:
            df[column] = 0
    return df


def extract_pbased_prona(df_source: pd.DataFrame, *args, **kwargs):
    """Extracts protein-based data for prona."""

    region_classes = ["PBR", "DBR", "RBR"]

    # Filter all regions based on allowed categories (only highest RI class)
    allowed_ri = [
        "Protein Binding (RI: 67-100)",
        "DNA Binding (RI: 67-100)",
        "RNA Binding (RI: 67-100)",
    ]
    df_source = df_source[df_source["description"].isin(allowed_ri)]

    # Filter (all) regions based on minlength requirement
    minlength = kwargs.pop("minlength")
    df_source = df_source[df_source["region length"] >= minlength]

    df_new = extract_pbased_all(df_source, *args, **kwargs)

    # Calculate region content as last step, when (correct) lengths of all proteins got collected
    df_new["region content"] = df_new["sum region lengths"] / df_new["protein length"]

    # Calculate additional region contents per region type
    region_type_counts = (
        df_source.groupby(["protein", "description"])["region length"]
        .sum()
        .unstack(level=-1)
    )

    df_new = df_new.join(region_type_counts, how="left").fillna(0)
    df_new = ensure_exists(allowed_ri, df_new)
    df_new = df_new.rename(
        columns={
            "Protein Binding (RI: 67-100)": "PBR content",
            "DNA Binding (RI: 67-100)": "DBR content",
            "RNA Binding (RI: 67-100)": "RBR content",
        }
    )
    for region_class in region_classes:
        df_new[f"{region_class} content"] = (
            df_new[f"{region_class} content"] / df_new["protein length"]
        )

    # Calculate number of regions per (binding) region type
    num_regions_counts = (
        df_source.groupby(["protein", "description"]).size().unstack(level=-1)
    )

    df_new = df_new.join(num_regions_counts, how="left").fillna(0)
    df_new = ensure_exists(allowed_ri, df_new)
    df_new = df_new.rename(
        columns={
            "Protein Binding (RI: 67-100)": "num PBR",
            "DNA Binding (RI: 67-100)": "num DBR",
            "RNA Binding (RI: 67-100)": "num RBR",
        }
    )

    return df_new


def extract_pbased_reprof(df_source: pd.DataFrame, *args, **kwargs):
    region_classes = ["Helix", "Strand", "Other"]

    print(df_source.head)

    # Filter all regions based on allowed categories (should already be this way)
    df_source = df_source[df_source["description"].isin(region_classes)]

    # Filter (all) regions based on minlength requirement
    minlength = kwargs.pop("minlength")
    df_source = df_source[df_source["region length"] >= minlength]

    df_new = extract_pbased_all(df_source, *args, **kwargs)

    # Calculate region content as last step, when (correct) lengths of all proteins got collected
    df_new["region content"] = df_new["sum region lengths"] / df_new["protein length"]

    # Calculate additional region contents per region type
    region_type_counts = (
        df_source.groupby(["protein", "description"])["region length"]
        .sum()
        .unstack(level=-1)
    )

    df_new = df_new.join(region_type_counts, how="left").fillna(0)
    df_new = ensure_exists(region_classes, df_new)
    for region_class in region_classes:
        df_new[f"{region_class}"] = df_new[f"{region_class}"] / df_new["protein length"]

    df_new = df_new.rename(
        columns={
            "Helix": "H",
            "Strand": "E",
            "Other": "O",
        }
    )

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
    # Add protein length from other dataframe
    df_new = df_new.join(df_seq["protein length"], how="right").fillna(0)

    return df_new


def extract_pbased(df_source: pd.DataFrame, df_seq: pd.DataFrame):
    """Maps required extraction parameters to each feature and collects p-based extraction results."""

    # Store function and required params for each feature
    mapping = {
        "tmseg": (
            extract_pbased_tmseg,
            {"minlength": 12, "region_type": ["Transmembrane Helix"]},
        ),
        "mdisorder": (extract_pbased_mdisorder, {"minlength": 30}),
        "prona": (
            extract_pbased_prona,
            {"minlength": 6, "region_type": ["Protein Binding (RI: 67-100)"]},
        ),
        "reprof": (extract_pbased_reprof, {"minlength": 4, "region_type": ["Helix"]}),
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


def extract_rbased_mdisorder(df_source: pd.DataFrame, *args, **kwargs):
    """"""

    df_new = extract_rbased_all(df_source, *args, **kwargs)
    return df_new


def extract_rbased_tmseg(df_source: pd.DataFrame, *args, **kwargs):
    """"""

    # Extract signal peptide containing proteins and shift the indices of their regions
    df_sig = df_source[df_source["description"] == "Signal Peptide"].drop_duplicates(
        "protein"
    )[["protein", "region length"]]
    df_sig = df_sig.rename(columns={"region length": "shift size"}).set_index("protein")
    df_source = df_source.join(df_sig, how="left", on="protein").fillna(0)
    df_source["begin"] = df_source["begin"] - df_source["shift size"]
    df_source["end"] = df_source["end"] - df_source["shift size"]

    # Clip protein length (if signal peptide present, else clip == 0) in sequences dataframe
    df_seq = args[0]
    df_seq = df_seq.join(df_sig, how="left").fillna(0)
    df_seq["protein length"] = df_seq["protein length"] - df_seq["shift size"]

    df_new = extract_rbased_all(df_source, df_seq, **kwargs)

    return df_new


def extract_rbased_prona(df_source: pd.DataFrame, *args, **kwargs):
    """"""

    df_new = extract_rbased_all(df_source, *args, **kwargs)
    return df_new


def extract_rbased_reprof(df_source: pd.DataFrame, *args, **kwargs):

    df_new = extract_rbased_all(df_source, *args, **kwargs)
    return df_new


def extract_rbased_all(
    df_source: pd.DataFrame,
    df_seq: pd.DataFrame,
    minlength: int,
    region_type: Optional[List[str]] = None,
):
    if region_type:
        # Remove unwanted region types
        df_source = df_source[df_source["description"].isin(region_type)]
    # Filter based on minlength requirement
    df_source = df_source[df_source["region length"] >= minlength]

    # Build region and length from begin and end
    df_new = pd.DataFrame()
    df_new["protein"] = df_source["protein"]
    df_source["reg length"] = df_source["end"] - df_source["begin"] + 1
    df_source["region"] = df_source[["begin", "end"]].apply(tuple, axis=1)
    df_new = df_new.join(df_source[["region", "reg length"]], how="left").fillna(0)
    # Add description (relevant if multiple types allowed)
    df_new = df_new.join(df_source["description"], how="left").fillna(0)

    # Calculate point region
    df_source = df_source.join(
        df_seq["protein length"], on="protein", how="left"
    ).fillna(0)
    df_source["point region"] = (
        df_source[["begin", "end"]]
        .div(df_source["protein length"].values, axis=0)
        .apply(tuple, axis=1)
    )

    # Add protein length and calculate relative region length
    df_new = df_new.join(
        df_source[["point region", "protein length"]], how="left"
    ).fillna(0)
    df_new["rel reg length"] = df_new["reg length"] / df_new["protein length"]

    return df_new


def extract_rbased(df_source: pd.DataFrame, df_seq: pd.DataFrame):
    """Maps required extraction parameters to each feature and collects r-based extraction results."""

    # Store function and required params for each feature
    mapping = {
        "tmseg": (
            extract_rbased_tmseg,
            {"minlength": 12, "region_type": ["Transmembrane Helix"]},
        ),
        "mdisorder": (extract_rbased_mdisorder, {"minlength": 30}),
        "prona": (
            extract_rbased_prona,
            {"minlength": 6, "region_type": ["Protein Binding (RI: 67-100)"]},
        ),
        "reprof": (
            extract_rbased_reprof,
            {"minlength": 4, "region_type": ["Helix", "Strand"]},
        ),
    }
    results = {}

    for feature, (func, kwargs) in mapping.items():
        # Reduce source dataframe to feature and calculate region lengths
        df_curr = df_source[df_source["feature"] == feature]
        df_curr["region length"] = df_curr["end"] - df_curr["begin"] + 1

        # Perform feature-specific extraction and store pbased results
        df_r = func(df_curr, df_seq, **kwargs)
        results[f"{feature} rbased"] = df_r

    return results
