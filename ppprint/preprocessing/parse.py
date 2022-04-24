"""
Parses some predict protein output.
This is mostly copypasta of `pp_to_json.py`, but for whole proteomes,
tailored to be used by ppprint.
"""

import json
import logging
from itertools import chain, groupby
from pathlib import Path
from typing import List, Callable

logger = logging.getLogger(__name__)


def retry_with_latin(f: Callable[..., List]):
    def inner(path, encoding="utf-8") -> List:
        try:
            return f(path)
        except UnicodeDecodeError:
            # If not utf-8, retry with latin-1
            logger.warning(f"File encoding is not utf-8, retrying with latin-1 for {path}")
            return f(path, encoding="latin-1")

    return inner


def group_segments(annotation):
    """Groups all consecutive elements of the same value in
    an annotation string and records the positions.
    """

    position = 1

    for type, g in groupby(annotation):
        # Amount of elements of the same value
        length = len(list(g))
        segment = dict()

        segment["type"] = type
        segment["begin"] = position
        segment["end"] = position + length - 1

        position = position + length

        yield segment


def filter_segments(segments, type_dict):
    """Sets type and description on filtered segments."""

    # Only allow the regions of a set of pre-defined types
    segments = [s for s in segments if s["type"] in type_dict]

    # According to type, set pre-defined description
    for segment in segments:
        segment["description"] = type_dict[segment.pop("type")]

    return segments


def get_sequence(path: Path) -> str:
    """Returns the sequence for the given .fasta file."""

    with path.open("r") as f:
        sequence = []

        count = 0
        for line in f:
            if not line.startswith(">"):
                sequence.append("".join(line.split()))
            else:
                count += 1

    # Warn if there is more than one sequence in the file
    if count > 1:
        logger.warning(f"FASTA file {path} contained {count} sequences. Expected 1.")
    return "".join(sequence)


@retry_with_latin
def parse_tmseg(path: Path, encoding="utf-8") -> List:
    """Parses a tmseg file."""

    # Set line categories as None to later discern lines
    header = None
    sequence = None
    annotation = None

    with path.open("r", encoding=encoding, errors="strict") as f:
        for line in f.readlines():
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

    if not annotation or not sequence or len(annotation) != len(sequence):
        # Annotation string is not present or does not fit sequence
        return []

    segments = group_segments(annotation)

    type_dict = {
        "S": "Signal Peptide",
        "H": "Transmembrane Helix",
        "1": "Cytoplasmic",
        "2": "Extracellular",
    }

    return filter_segments(segments, type_dict)


@retry_with_latin
def parse_prona(path: Path, encoding="utf-8") -> List:
    """Parses a prona file."""

    prona_pro = []
    prona_dna = []
    prona_rna = []

    with path.open("r", encoding=encoding, errors="strict") as input_file:
        for line in input_file:
            if line.lstrip().startswith("Res_"):
                data = line.split()

                if len(data) == 8:
                    if data[3] != "0":
                        # Residue predicted as protein binding
                        pro_ri = int(data[2])

                        # Category according to RI class
                        if pro_ri <= -1:
                            value = "X"
                        elif pro_ri <= 33:
                            value = "P0"
                        elif pro_ri <= 66:
                            value = "P1"
                        else:
                            value = "P2"
                        prona_pro.append(value)
                    else:
                        # Not protein binding
                        prona_pro.append(".")

                    if data[5] != "0":
                        # Residue predicted as DNA binding
                        dna_ri = int(data[4])

                        # Category according to RI class
                        if dna_ri <= -1:
                            value = "X"
                        elif dna_ri <= 33:
                            value = "D0"
                        elif dna_ri <= 66:
                            value = "D1"
                        else:
                            value = "D2"
                        prona_dna.append(value)
                    else:
                        # Not DNA binding
                        prona_dna.append(".")

                    if data[7] != "0":
                        # Residue predicted as RNA binding
                        rna_ri = int(data[6])

                        # Category according to RI class
                        if rna_ri <= -1:
                            value = "X"
                        elif rna_ri <= 33:
                            value = "R0"
                        elif rna_ri <= 66:
                            value = "R1"
                        else:
                            value = "R2"

                        prona_rna.append(value)
                    else:
                        # Not RNA binding
                        prona_rna.append(".")
                else:
                    break

    segments = chain(
        group_segments(prona_pro), group_segments(prona_dna), group_segments(prona_rna)
    )

    type_dict = {
        "P0": "Protein Binding (RI: 00-33)",
        "P1": "Protein Binding (RI: 34-66)",
        "P2": "Protein Binding (RI: 67-100)",
        "D0": "DNA Binding (RI: 00-33)",
        "D1": "DNA Binding (RI: 34-66)",
        "D2": "DNA Binding (RI: 67-100)",
        "R0": "RNA Binding (RI: 00-33)",
        "R1": "RNA Binding (RI: 34-66)",
        "R2": "RNA Binding (RI: 67-100)",
    }

    return filter_segments(segments, type_dict)


@retry_with_latin
def parse_mdisorder(path: Path, encoding="utf-8") -> List:
    """Parses a given mdisorder file."""

    num_col = -1
    mdisorder = []

    with path.open("r", encoding=encoding, errors="strict") as input_file:
        for line in input_file:
            if num_col > 0:
                data = line.split()

                if len(data) == num_col:
                    # Build annotation only from two-state prediction by MetaDisorder
                    mdisorder.append(data[10])
                else:
                    break
            elif line.lstrip().startswith("Number"):
                # Reached header
                data = line.split()
                num_col = len(data)

    segments = group_segments(mdisorder)

    type_dict = {"D": "Disordered Region"}

    return filter_segments(segments, type_dict)


@retry_with_latin
def parse_reprof(path: Path, encoding="utf-8") -> List:
    """
        Parses a given reprof file.
        Only the structure information is retained.
    """
    num_col = -1
    structure = []

    with path.open("r", encoding=encoding, errors="strict") as input_file:
        for line in input_file:
            if num_col > 0:
                data = line.split()

                if len(data) == num_col:
                    structure.append(data[2])
                else:
                    break
            elif line.lstrip().startswith("#"):
                continue
            elif line.lstrip().startswith("No"):
                data = line.split()
                num_col = len(data)

    segments = group_segments(structure)

    type_dict = {
        "H": "Helix",
        "E": "Strand",
        "L": "Other",
    }

    return filter_segments(segments, type_dict)


def parse_protein(base_path: Path, protein: str):
    """Parses information for a given protein."""

    def run(f, extension: str):
        """Runs a parser function with the specified file extension."""

        path = base_path / f"{protein}.{extension}"

        if path.exists():
            return f(path)
        else:
            return []

    # We can trust that sequence never gets a list
    sequence = run(get_sequence, "fasta")
    tmseg = run(parse_tmseg, "tmseg")
    prona = run(parse_prona, "prona")
    mdisorder = run(parse_mdisorder, "mdisorder")
    reprof = run(parse_reprof, "reprof")

    return (sequence, tmseg, prona, mdisorder, reprof)


def parse(base_path: Path):
    """Parses sequence, tmseg, prona and mdisorder and returns a generator."""

    for p in (p for p in base_path.iterdir() if p.is_dir()):
        # Identify the proteins based on the present .fasta files
        # (required existence of a .fasta for each protein)
        fasta_files = list(p.glob("*.fasta"))
        for protein in (p.stem for p in fasta_files):
            yield parse_protein(p, protein)


def write_json(base_path: Path, out_file: Path):
    """Writes a JSON file for a given proteome."""

    json.JSONEncoder(ensure_ascii=False, check_circular=False)
    with open(out_file, "w") as f:
        # We have to store everything in a list to make the json encoder happy :(
        # (Actually we don't have to, but it's required on loading anyways)
        data = [
            {
                "sequence": sequence,
                "tmseg": tmseg,
                "prona": prona,
                "mdisorder": mdisorder,
                "reprof": reprof,
            }
            for sequence, tmseg, prona, mdisorder, reprof in parse(base_path)
        ]

        for chunk in json.JSONEncoder().iterencode(data):
            f.write(chunk)


if __name__ == "__main__":
    write_json(Path("data"), Path("outfile.json"))
