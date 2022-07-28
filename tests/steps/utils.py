import os
import pandas as pd

import numpy as np

def build_true_segments_json():
    """Manually build ground truth JSON from PredictProtein output files."""

    # Descriptions
    proteins = {
        "P62524": "MTALLRVISLVVISVVVIIIPPCGAALGRGKA",
        "Q8XA85": "MKRISTTITTTITTTITITITTGNGAG",
    }
    # Mdisorder
    dr = "Disordered Region"
    # Tmseg
    sp = "Signal Peptide"
    tmh = "Transmembrane Helix"
    extr = "Extracellular"
    cyto = "Cytoplasmic"
    # Prona
    pb33 = "Protein Binding (RI: 00-33)"
    pb66 = "Protein Binding (RI: 34-66)"
    pb100 = "Protein Binding (RI: 67-100)"
    # Reprof
    h = "Helix"
    e = "Strand"
    o = "Other"

    mdisorder = {
        "P62524": [
            (1, 7, dr),
            (9, 11, dr),
            (23, 23, dr),
            (25, 32, dr),
        ],
        "Q8XA85": [
            (1, 27, dr),
        ],
    }
    tmseg = {
        "P62524": [
            (1, 3, extr),
            (4, 27, tmh),
            (28, 32, cyto),
        ],
        "Q8XA85": [
            (1, 6, sp),
        ],
    }
    prona = {
        "P62524": [
            (18, 19, pb33),
            (20, 24, pb66),
            (25, 27, pb33),
        ],
        "Q8XA85": [
            (16, 19, pb33),
            (20, 20, pb66),
            (21, 22, pb33),
        ],
    }
    reprof = {
        "P62524": [
            (1, 1, o),
            (2, 28, h),
            (29, 32, o),
        ],
        "Q8XA85": [
            (1, 2, o),
            (3, 4, e),
            (5, 5, o),
            (6, 22, e),
            (23, 27, o),
        ],
    }

    # Features in the order the parser builds
    features = {
        "tmseg": tmseg,
        "prona": prona,
        "mdisorder": mdisorder,
        "reprof": reprof,
    }

    def build_segments(segments):
        segment_list = [
            {"begin": t[0], "end": t[1], "description": t[2]} for t in segments
        ]
        return segment_list

    data = []
    for p, seq in proteins.items():
        entry = dict()
        entry["sequence"] = seq
        for feature_name, feature_dict in features.items():
            entry[feature_name] = build_segments(feature_dict[p])

        data.append(entry)

    return data


def convert_mdisorder_to_latin1(data_folder):
    """Converts .mdisorder files to latin-1 encoding."""

    for folder in (folder for folder in data_folder.iterdir() if folder.is_dir()):
        for pp_file in os.listdir(folder):
            if pp_file.endswith(".mdisorder"):
                with open((data_folder / folder / pp_file), "r", encoding="utf-8") as f:
                    content = f.read()
                with open(
                    (data_folder / folder / pp_file), "w", encoding="latin-1"
                ) as f:
                    f.write(content)


def build_true_dfs_pbased():
    """Manually builds ground truth DataFrames on protein basis for extraction tests."""

    results = {}

    # Mdisorder pbased
    df = pd.DataFrame(
        columns=[
            "number of regions",
            "median length",
            "sum region lengths",
            "protein length",
            "region content",
        ]
    )
    # df = pd.DataFrame(
    #     {
    #         "number of regions": pd.Series(dtype=int),
    #         "median length": pd.Series(dtype=float),
    #         "sum region lengths": pd.Series(dtype=int),
    #         "protein length": pd.Series(dtype=int),
    #         "region content": pd.Series(dtype=float),
    #     }
    # )
    df.loc[0] = [1, 30, 30, 60, 0.5]
    df.loc[1] = [0, 0, 0, 70, 0]
    df = df.astype(
        {
            "number of regions": "int64",
            "median length": "float64",
            "sum region lengths": "int64",
            "protein length": "int64",
            "region content": "float64",
        }
    )
    results["mdisorder pbased"] = df

    # Tmseg pbased
    df = pd.DataFrame(
        columns=[
            "number of regions",
            "median length",
            "sum region lengths",
            "protein length",
            "shift size",
            "region content",
            "I",
            "O",
            "M",
            "orientation",
        ]
    )
    # df = pd.DataFrame(
    #     {
    #         "number of regions": pd.Series(dtype=int),
    #         "median length": pd.Series(dtype=float),
    #         "sum region lengths": pd.Series(dtype=int),
    #         "protein length": pd.Series(dtype=int),
    #         "shift size": pd.Series(dtype=int),
    #         "region content": pd.Series(dtype=float),
    #         "I": pd.Series(dtype=float),
    #         "O": pd.Series(dtype=float),
    #         "M": pd.Series(dtype=float),
    #         "orientation": pd.Series(dtype=str),
    #     }
    # )
    df.loc[0] = [
        1,
        12,
        12,
        60,
        0,
        0.2,
        8 / 60,
        33 / 60,
        12 / 60,
        "Extracellular",
    ]
    # FIXME: Is this correct behavior? Clipped protein length for content calculation, but counts SPs?
    df.loc[1] = [
        2,
        12,
        24,
        (70 - 10),
        10,
        0.4,
        32 / 60,
        4 / 60,
        34 / 60,
        "Signal Peptide",
    ]
    df = df.astype(
        {
            "number of regions": "int64",
            "median length": "float64",
            "sum region lengths": "int64",
            "protein length": "int64",
            "shift size": "int64",
            "region content": "float64",
            "I": "float64",
            "O": "float64",
            "M": "float64",
            "orientation": "str",
        }
    )
    results["tmseg pbased"] = df

    # Prona pbased
    df = pd.DataFrame(
        columns=[
            "number of regions",
            "median length",
            "sum region lengths",
            "protein length",
            "region content",
            "DBR content",
            "PBR content",
            "RBR content",
            "num DBR",
            "num PBR",
            "num RBR",
        ]
    )
    # df = pd.DataFrame(
    #     {
    #         "number of regions": pd.Series(dtype=int),
    #         "median length": pd.Series(dtype=float),
    #         "sum region lengths": pd.Series(dtype=int),
    #         "protein length": pd.Series(dtype=int),
    #         "region content": pd.Series(dtype=float),
    #         "DBR content": pd.Series(dtype=float),
    #         "PBR content": pd.Series(dtype=float),
    #         "RBR content": pd.Series(dtype=float),
    #         "num DBR": pd.Series(dtype=int),
    #         "num PBR": pd.Series(dtype=int),
    #         "num RBR": pd.Series(dtype=int),
    #     }
    # )
    df.loc[0] = [1, 6, 6, 60, 0.1, 0.1, 0.1, 0.1, 1, 1, 1]
    df.loc[1] = [0, 0, 0, 70, 0, 0, 0, 0, 0, 0, 0]
    df = df.astype(
        {
            "number of regions": "int64",
            "median length": "float64",
            "sum region lengths": "int64",
            "protein length": "int64",
            "region content": "float64",
            "DBR content": "float64",
            "PBR content": "float64",
            "RBR content": "float64",
            "num DBR": "int64",
            "num PBR": "int64",
            "num RBR": "int64",
        }
    )
    results["prona pbased"] = df

    # Reprof pbased
    df = pd.DataFrame(
        columns=[
            "number of regions",
            "median length",
            "sum region lengths",
            "protein length",
            "region content",
            "H",
            "O",
            "E",
        ]
    )
    # df = pd.DataFrame(
    #     {
    #         "number of regions": pd.Series(dtype=int),
    #         "median length": pd.Series(dtype=float),
    #         "sum region lengths": pd.Series(dtype=int),
    #         "protein length": pd.Series(dtype=int),
    #         "region content": pd.Series(dtype=float),
    #         "H": pd.Series(dtype=float),
    #         "O": pd.Series(dtype=float),
    #         "E": pd.Series(dtype=float),
    #     }
    # )
    # FIXME: Apparently, 4 is minimum length for all components, pbased and rbased. Sensible?
    df.loc[0] = [1, 10, 10, 60, 1 / 6, 1 / 6, 34 / 60, 1 / 6]
    df.loc[1] = [0, 0, 0, 70, 0, 0, 0, 0]
    df = df.astype(
        {
            "number of regions": "int64",
            "median length": "float64",
            "sum region lengths": "int64",
            "protein length": "int64",
            "region content": "float64",
            "H": "float64",
            "O": "float64",
            "E": "float64",
        }
    )
    results["reprof pbased"] = df

    return results


def build_true_dfs_rbased():
    """Manually builds ground truth DataFrames on region basis for extraction tests."""

    results = {}

    features = ["mdisorder", "tmseg", "prona", "reprof"]
    for feature in features:
        df = pd.DataFrame(
            columns=[
                "protein",
                "region",
                "reg length",
                "description",
                "point region",
                "protein length",
                "rel reg length",
            ]
        )

        if feature == "mdisorder":
            # Mdisorder
            df.loc[0] = [0, np.array((1, 30)), 30, "Disordered Region", np.array((1 / 60, 0.5)), 60, 0.5]
        elif feature == "tmseg":
            # Tmseg
            df.loc[0] = [
                0,
                np.array((13, 24)),
                12,
                "Transmembrane Helix",
                np.array((13 / 60, 24 / 60)),
                60,
                0.2,
            ]
            df.loc[1] = [
                1,
                np.array((1, 12)),
                12,
                "Transmembrane Helix",
                np.array((1 / 60, 12 / 60)),
                60,
                0.2,
            ]
            df.loc[2] = [
                1,
                np.array((17, 28)),
                12,
                "Transmembrane Helix",
                np.array((17 / 60, 28 / 60)),
                60,
                0.2,
            ]
        elif feature == "prona":
            # Prona
            df.loc[0] = [
                0,
                np.array((1, 6)),
                6,
                "Protein Binding (RI: 67-100)",
                np.array((1 / 60, 0.1)),
                60,
                0.1,
            ]
        else:
            # Reprof
            df.loc[0] = [0, np.array((1, 10)), 10, "Helix", np.array((1 / 60, 1 / 6)), 60, 1 / 6]
            df.loc[1] = [0, np.array((11, 20)), 10, "Strand", np.array((11 / 60, 1 / 3)), 60, 1 / 6]

        df = df.astype(
            {
                "protein": "int64",
                #"region": "float64",
                "reg length": "int64",
                "description": "str",
                #"point region": "float64",
                "protein length": "int64",
                "rel reg length": "float64",
            }
        )

        results[f"{feature} rbased"] = df

    return results
