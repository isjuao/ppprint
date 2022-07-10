import os


def build_true_segments_json():
    """Manually build ground truth JSON from PredictProtein output files."""

    # Descriptions
    proteins = {"P62524": "MTALLRVISLVVISVVVIIIPPCGAALGRGKA", "Q8XA85": "MKRISTTITTTITTTITITITTGNGAG", }
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
        "P62524": [(1, 7, dr), (9, 11, dr), (23, 23, dr), (25, 32, dr), ],
        "Q8XA85": [(1, 27, dr), ],
    }
    tmseg = {
        "P62524": [(1, 3, extr), (4, 27, tmh), (28, 32, cyto), ],
        "Q8XA85": [(1, 6, sp), ],
    }
    prona = {
        "P62524": [(18, 19, pb33), (20, 24, pb66), (25, 27, pb33), ],
        "Q8XA85": [(16, 19, pb33), (20, 20, pb66), (21, 22, pb33), ],
    }
    reprof = {
        "P62524": [(1, 1, o), (2, 28, h), (29, 32, o), ],
        "Q8XA85": [(1, 2, o), (3, 4, e), (5, 5, o), (6, 22, e), (23, 27, o), ],
    }

    # Features in the order the parser builds
    features = {
        "tmseg": tmseg,
        "prona": prona,
        "mdisorder": mdisorder,
        "reprof": reprof,
    }

    def build_segments(segments):
        segment_list = [{"begin": t[0], "end": t[1], "description": t[2]} for t in segments]
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
                with open((data_folder / folder / pp_file), "w", encoding="latin-1") as f:
                    f.write(content)