#!/usr/bin/env python3
import logging
import re
import json
import random
import argparse
import requests
import subprocess

from pathlib import Path
from itertools import groupby


HEX_COL1 = "#648FFF"
HEX_COL2 = "#785EF0"
HEX_COL3 = "#DC267F"
HEX_COL4 = "#FE6100"
HEX_COL5 = "#FFB000"


def __get_fasta(uniprot_id, temp_dir="./"):
    uniprot_url = "https://www.uniprot.org/uniprot/{}.fasta".format(uniprot_id)

    if not Path(temp_dir).exists():
        Path(temp_dir).mkdir(parents=True)

    fasta_file = Path(
        temp_dir, "{}_{:09}.fasta".format(uniprot_id, random.randrange(1000000000))
    )

    try:
        remote = requests.get(uniprot_url, stream=True, timeout=10)

        with fasta_file.open("wb") as ff:
            ff.write(remote.content)

        download_complete = True
    except:
        download_complete = False

    if download_complete:
        sequence = __parse_fasta(fasta_file)

    if fasta_file.exists():
        fasta_file.unlink()

    return sequence


def __parse_fasta(fasta_file):
    sequence = []

    with fasta_file.open("r") as ff:
        for line in ff:
            if line.startswith(">"):
                sequence = []
            else:
                sequence.append("".join(line.split()))

    sequence = "".join(sequence)

    return sequence


def __get_ppc_root_dir(sequence, ppc_hash):
    if sequence is not None:
        params = ["ppc_fetch", "--skip-hash-lock", "-d", "--seq", sequence]
    else:
        params = ["ppc_fetch", "--skip-hash-lock", "-d", "--hash", ppc_hash]

    try:
        ppc_root_dir = subprocess.check_output(params).decode("utf-8").strip()
    except subprocess.CalledProcessError:
        ppc_root_dir = None

    if ppc_root_dir is not None:
        ppc_root_dir = Path(ppc_root_dir)

        if ppc_root_dir.exists() and ppc_root_dir.is_dir():
            return ppc_root_dir

    return None


def __group_segments(annotation):
    position = 1
    segments = []

    for k, g in groupby(annotation):
        length = len(list(g))
        segment = dict()

        segment["type"] = k
        segment["category"] = "PREDICT_PROTEIN"
        segment["begin"] = position
        segment["end"] = position + length - 1

        position = position + length

        segments.append(segment)

    return segments


def __filter_segments(segments, type_dict):
    segments = [s for s in segments if s["type"] in type_dict]

    for segment in segments:
        seg_type, seg_desc, seg_color = type_dict[segment["type"]]

        segment["type"] = seg_type
        segment["description"] = seg_desc
        segment["color"] = seg_color

    return segments


def __parse_mdisorder(ppc_root_dir):
    mdisorder_file = Path(ppc_root_dir, "query.mdisorder")  # query.mdisorder

    if not mdisorder_file.exists():
        return []

    num_col = -1
    mdisorder = []

    with mdisorder_file.open("r") as input_file:
        for line in input_file:
            if num_col > 0:
                data = line.split()

                if len(data) == num_col:
                    mdisorder.append(data[10])
                else:
                    break
            elif line.lstrip().startswith("Number"):
                data = line.split()
                num_col = len(data)

    segments = __group_segments(mdisorder)

    type_dict = {
        "D": ("DISORDERED_REGION_(META-DISORDER)", "Disordered Region", HEX_COL3)
    }

    return __filter_segments(segments, type_dict)


def __parse_norsnet(ppc_root_dir):
    norsnet_file = Path(ppc_root_dir, "query.norsnet")

    if not norsnet_file.exists():
        return []

    num_col = -1
    norsnet = []

    with norsnet_file.open("r") as input_file:
        for line in input_file:
            if num_col > 0:
                data = line.split()

                if len(data) == num_col:
                    norsnet.append(data[6])
                else:
                    break
            elif line.lstrip().startswith("pos"):
                data = line.split()
                num_col = len(data)

    segments = __group_segments(norsnet)

    type_dict = {
        "N": (
            "NON-REGULAR_SECONDARY_STRUCTURE_(NORS-NET)",
            "Non-Regular Secondary Structure (NORS)",
            HEX_COL3,
        )
    }

    return __filter_segments(segments, type_dict)


def __parse_profbval(ppc_root_dir):
    profbval_file = Path(ppc_root_dir, "query.profbval")

    if not profbval_file.exists():
        return []

    profbval = []
    read_line = False

    with profbval_file.open("r") as input_file:
        for line in input_file:
            if read_line:
                data = line.split()

                if len(data) == 3:
                    profbval.append(int(data[1]))
                else:
                    break
            elif line.lstrip().startswith("* out vec: (I8,A1,100I4)"):
                read_line = True

    for i in range(len(profbval)):
        value = profbval[i]

        if value <= 30:
            profbval[i] = "0"
        elif value <= 70:
            profbval[i] = "1"
        else:
            profbval[i] = "2"

    segments = __group_segments(profbval)

    type_dict = {
        "0": (
            "RELATIVE_B-VALUE_(PROF-BVAL)",
            "Relative B-Value: 00-30 (Low)",
            HEX_COL1,
        ),
        "1": (
            "RELATIVE_B-VALUE_(PROF-BVAL)",
            "Relative B-Value: 31-70 (Intermediate)",
            HEX_COL3,
        ),
        "2": (
            "RELATIVE_B-VALUE_(PROF-BVAL)",
            "Relative B-Value: 71-99 (High)",
            HEX_COL5,
        ),
    }

    return __filter_segments(segments, type_dict)


def __parse_tmseg(ppc_root_dir):
    tmseg_file = Path(ppc_root_dir, "query.tmseg")

    if not tmseg_file.exists():
        return []

    header = None
    sequence = None
    annotation = None

    # logger = logging.getLogger(__name__)
    # logger.setLevel(logging.INFO)
    # tmseg_file_content = open(tmseg_file)
    # logger.warning(msg=str(tmseg_file_content.read()))
    with tmseg_file.open(
        "r", encoding="latin-1"
    ) as input_file:  # sometimes needs latin-1 (encoding="latin-1")
        for line in input_file:
            if line.startswith("#"):
                continue
            elif line.startswith(">"):
                header = line.strip()
            elif header and sequence is None:
                sequence = line.strip()
            elif sequence and annotation is None:
                annotation = line.strip()
            else:
                break

    if annotation is None or len(annotation) != len(sequence):
        return []

    segments = __group_segments(annotation)

    type_dict = {
        "S": ("TOPOLOGY_(TMSEG)", "Signal Peptide", HEX_COL3),
        "H": ("TOPOLOGY_(TMSEG)", "Transmembrane Helix", HEX_COL1),
        "1": ("TOPOLOGY_(TMSEG)", "Cytoplasmic", HEX_COL4),
        "2": ("TOPOLOGY_(TMSEG)", "Extracellular", HEX_COL5),
    }

    return __filter_segments(segments, type_dict)


def __parse_reprof(ppc_root_dir):
    reprof_file = Path(ppc_root_dir, "query.reprof")

    if not reprof_file.exists():
        return []

    num_col = -1
    structure = []
    accessibility = []

    with reprof_file.open("r") as input_file:  #  encoding="latin-1"
        for line in input_file:
            if num_col > 0:
                data = line.split()

                if len(data) == num_col:
                    structure.append(data[2])
                    # col 11 only has buried/exposed
                    accessibility.append(data[11])
                else:
                    break
            elif line.lstrip().startswith("#"):
                continue
            elif line.lstrip().startswith("No"):
                data = line.split()
                num_col = len(data)

    segments = []

    segments.extend(__group_segments(structure))
    segments.extend(__group_segments(accessibility))

    type_dict = {
        "H": ("SECONDARY_STRUCTURE_(REPROF)", "Helix", HEX_COL1),
        "E": ("SECONDARY_STRUCTURE_(REPROF)", "Strand", HEX_COL3),
        "L": ("SECONDARY_STRUCTURE_(REPROF)", "Other", HEX_COL5),
        "e": ("SOLVENT_ACCESSIBILITY_(REPROF)", "Exposed", HEX_COL1),
        "i": ("SOLVENT_ACCESSIBILITY_(REPROF)", "Intermediate", HEX_COL3),
        "b": ("SOLVENT_ACCESSIBILITY_(REPROF)", "Buried", HEX_COL5),
    }

    return __filter_segments(segments, type_dict)


def __parse_consurf(ppc_root_dir):
    consurf_file = Path(ppc_root_dir, "query.consurf.grades")

    if not consurf_file.exists():
        return []

    consurf = []
    num_col = -1
    skiped_lines = 0

    with consurf_file.open("r") as input_file:
        for line in input_file:
            if num_col > 0:
                data = line.split()

                if len(data) >= num_col:
                    consurf.append(int(data[3][0]))
                elif skiped_lines == 0:
                    skiped_lines = 1
                else:
                    break
            elif line.lstrip().startswith("POS"):
                data = line.split()
                num_col = len(data) - 5

    for i in range(len(consurf)):
        value = consurf[i]

        if value <= 3:
            consurf[i] = "0"
        elif value <= 6:
            consurf[i] = "1"
        else:
            consurf[i] = "2"

    consurf = "".join(consurf)

    segments = []

    segments.extend(__group_segments(consurf))

    type_dict = {
        "0": ("CONSERVATION_(CONSEQ)", "Conservation Score: 1-3 (Low)", HEX_COL1),
        "1": (
            "CONSERVATION_(CONSEQ)",
            "Conservation Score: 4-6 (Intermediate)",
            HEX_COL3,
        ),
        "2": ("CONSERVATION_(CONSEQ)", "Conservation Score: 7-9 (High)", HEX_COL5),
    }

    return __filter_segments(segments, type_dict)


def __parse_prona(ppc_root_dir):
    prona_file = Path(ppc_root_dir, "query.prona")

    if not prona_file.exists():
        return []

    prona_pro = []
    prona_dna = []
    prona_rna = []

    with prona_file.open("r") as input_file:
        for line in input_file:
            if line.lstrip().startswith("Res_"):
                data = line.split()

                if len(data) == 8:
                    pro = bool(int(data[3]))
                    dna = bool(int(data[5]))
                    rna = bool(int(data[7]))

                    pro_ri = int(data[2])
                    dna_ri = int(data[4])
                    rna_ri = int(data[6])

                    if pro:
                        prona_pro.append(pro_ri)
                    else:
                        prona_pro.append(None)

                    if dna:
                        prona_dna.append(dna_ri)
                    else:
                        prona_dna.append(None)

                    if rna:
                        prona_rna.append(rna_ri)
                    else:
                        prona_rna.append(None)
                else:
                    break

    for i in range(len(prona_pro)):
        value = prona_pro[i]

        if value is None:
            prona_pro[i] = "."
        elif value <= -1:
            prona_pro[i] = "X"
        elif value <= 33:
            prona_pro[i] = "P0"
        elif value <= 66:
            prona_pro[i] = "P1"
        else:
            prona_pro[i] = "P2"

    for i in range(len(prona_dna)):
        value = prona_dna[i]

        if value is None:
            prona_dna[i] = "."
        elif value <= -1:
            prona_dna[i] = "X"
        elif value <= 33:
            prona_dna[i] = "D0"
        elif value <= 66:
            prona_dna[i] = "D1"
        else:
            prona_dna[i] = "D2"

    for i in range(len(prona_rna)):
        value = prona_rna[i]

        if value is None:
            prona_rna[i] = "."
        elif value <= -1:
            prona_rna[i] = "X"
        elif value <= 33:
            prona_rna[i] = "R0"
        elif value <= 66:
            prona_rna[i] = "R1"
        else:
            prona_rna[i] = "R2"

    segments = []

    segments.extend(__group_segments(prona_pro))
    segments.extend(__group_segments(prona_dna))
    segments.extend(__group_segments(prona_rna))

    type_dict = {
        "P0": ("PROTEIN_BINDING_(PRONA)", "Protein Binding (RI: 00-33)", HEX_COL1),
        "P1": ("PROTEIN_BINDING_(PRONA)", "Protein Binding (RI: 34-66)", HEX_COL3),
        "P2": ("PROTEIN_BINDING_(PRONA)", "Protein Binding (RI: 67-100)", HEX_COL5),
        "D0": ("DNA_BINDING_(PRONA)", "DNA Binding (RI: 00-33)", HEX_COL1),
        "D1": ("DNA_BINDING_(PRONA)", "DNA Binding (RI: 34-66)", HEX_COL3),
        "D2": ("DNA_BINDING_(PRONA)", "DNA Binding (RI: 67-100)", HEX_COL5),
        "R0": ("RNA_BINDING_(PRONA)", "RNA Binding (RI: 00-33)", HEX_COL1),
        "R1": ("RNA_BINDING_(PRONA)", "RNA Binding (RI: 34-66)", HEX_COL3),
        "R2": ("RNA_BINDING_(PRONA)", "RNA Binding (RI: 67-100)", HEX_COL5),
    }

    return __filter_segments(segments, type_dict)


def __parse_disulfind(ppc_root_dir):
    disulfind_file = Path(ppc_root_dir, "query.disulfinder")

    if not disulfind_file.exists():
        return []

    dsf_bonds = []
    bond_pattern = re.compile(r"^.*>DB_bond<.*bond\((.*),(.*)\)")

    with disulfind_file.open("r") as input_file:
        for line in input_file:
            line = line.strip()
            match = bond_pattern.match(line)

            if match:
                p1 = int(match.group(1).strip())
                p2 = int(match.group(2).strip())

                dsf_bonds.append((p1, p2))

    segments = []

    for p1, p2 in dsf_bonds:
        segment = dict()

        segment["type"] = "B"
        segment["category"] = "PREDICT_PROTEIN"
        segment["begin"] = p1
        segment["end"] = p2

        segments.append(segment)

    type_dict = {"B": ("DISULFIDE_BOND_(DISULFIND)", "Disulfide Bond", HEX_COL1)}

    return __filter_segments(segments, type_dict)


def __add_evidences(pp_json, ppc_hash):
    source = {
        "name": "PredictProtein",
        "id": ppc_hash,
        "url": "https://open.predictprotein.org/visual_results?req_id={}".format(
            ppc_hash
        ),
    }
    evidence = {"code": "ECO:0007669", "source": source}

    for feature in pp_json["features"]:
        feature["evidences"] = [evidence]


def __raise_error(error_code, error_message):
    error_titles = {400: "Bad Request", 404: "Not Found", 500: "Internal Server Error"}

    error_title = error_titles[error_code]

    error_json = dict()
    erorr_content = dict()

    error_json["errors"] = [erorr_content]
    erorr_content["status"] = str(error_code)
    erorr_content["title"] = error_title
    erorr_content["detail"] = error_message

    print(json.dumps(error_json, indent=None, separators=(",", ":")))


def __check_hash(ppc_hash):
    hash_pattern = r"[A-Za-z0-9]{40}"

    if re.match(hash_pattern, ppc_hash):
        return True
    else:
        __raise_error(400, "Malformed PPC hash.")
        return False


def __check_id(uniprot_id):
    id_pattern = r"[A-Z0-9]{1,10}_[A-Z0-9]{1,5}"
    ac_pattern = r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"

    if re.match(ac_pattern, uniprot_id) or re.match(id_pattern, uniprot_id):
        return True
    else:
        __raise_error(400, "Malformed UniProt ID/AC.")
        return False


def __check_sequence(sequence):
    if sequence is None:
        __raise_error(500, "Could not retrieve protein sequence.")
        return False

    if re.match(r"[A-Za-z]+", sequence):
        return True
    else:
        __raise_error(400, "Malformed protein sequence.")
        return False


def __check_ppc_dir(ppc_root_dir):
    files = [
        "query.tmseg",
        "query.reprof",
        "query.mdisorder",
        "query.norsnet",
        "query.profbval",
        "query.consurf.grades",
        "query.prona",
        "query.disulfinder",
    ]

    if ppc_root_dir is not None:
        files_found = [Path(ppc_root_dir, file).exists() for file in files]

        if any(files_found):
            return True

    __raise_error(404, "No results available for this sequence.")
    return False


def main(ppc_hash, uniprot_id, fasta_file, methods, compact_json, temp_dir):
    if ppc_hash is not None:
        sequence = None

        if not __check_hash(ppc_hash):
            return 1
    elif uniprot_id is not None:
        uniprot_id = uniprot_id.upper()

        if not __check_id(uniprot_id):
            return 1

        sequence = __get_fasta(uniprot_id, temp_dir)

        if not __check_sequence(sequence):
            return 1
    elif fasta_file is not None:
        ff = Path(fasta_file)

        if not ff.exists():
            __raise_error(404, "FASTA file not found.")
            return 1

        sequence = __parse_fasta(ff)

        if not __check_sequence(sequence):
            return 1
    else:
        __raise_error(400, "No target information provided.")
        return 1

    ppc_root_dir = __get_ppc_root_dir(sequence, ppc_hash)
    # ppc_root_dir = "../"

    if not __check_ppc_dir(ppc_root_dir):
        return 1

    if sequence is None:
        sequence = __parse_fasta(Path(ppc_root_dir, "query.fasta"))

    parsers = [
        ("tmseg", __parse_tmseg),
        ("reprof", __parse_reprof),
        ("mdisorder", __parse_mdisorder),
        ("norsnet", __parse_norsnet),
        ("profbval", __parse_profbval),
        ("consurf", __parse_consurf),
        ("prona", __parse_prona),
        ("disulfind", __parse_disulfind),
    ]

    features = []

    if methods is None:
        for name, func in parsers:
            features.extend(func(ppc_root_dir))
    else:
        for name, func in parsers:
            if name in methods:
                features.extend(func(ppc_root_dir))

    pp_json = dict()

    pp_json["sequence"] = sequence
    pp_json["features"] = features

    __add_evidences(pp_json, ppc_root_dir.name)

    if compact_json:
        print(json.dumps(pp_json, indent=None, separators=(",", ":")))
    else:
        print(json.dumps(pp_json, indent=4, separators=(",", ": ")))

    return 0


if __name__ == "__main__":
    args_parser = argparse.ArgumentParser()

    args_parser.add_argument("-x", "--ppc-hash", default=None)
    args_parser.add_argument("-i", "--uniprot-id", default=None)
    args_parser.add_argument("-f", "--fasta-file", default=None)
    args_parser.add_argument("-m", "--methods", default=None)
    args_parser.add_argument("-t", "--temp-dir", default="./")
    args_parser.add_argument("-c", "--compact", action="store_true", default=False)

    args = args_parser.parse_args()

    ppc_hash = args.ppc_hash
    uniprot_id = args.uniprot_id
    fasta_file = args.fasta_file
    methods = args.methods
    temp_dir = args.temp_dir
    compact_json = args.compact

    if methods is not None:
        methods = set([m.strip() for m in methods.split(",")])

    main(ppc_hash, uniprot_id, fasta_file, methods, compact_json, temp_dir)
