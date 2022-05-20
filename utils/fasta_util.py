import gzip
import os
import subprocess
from dataclasses import dataclass
from typing import Iterator, List

from Bio import SeqIO

from utils import file_util

GZIP_SUFFIXES = [".gz", ".gzip"]
FASTA_SUFFIXES = [".fasta", ".fna", ".ffn", ".faa", ".frn"]


@dataclass
class GffRow:
    seqid: str
    source: str
    type_: str
    start: str
    end: str
    score: str
    strand: str
    phase: str
    attributes: str


# TODO: support cloud files automatically?
def open_fasta(fpath) -> Iterator[SeqIO.SeqRecord]:
    """Iterates of the records in a fasta file (even if it's gzipped)."""
    is_gzipped = False
    for gzip_suffix in GZIP_SUFFIXES:
        if fpath.lower().endswith(gzip_suffix):
            is_gzipped = True
            decompressed_fpath = fpath[: -len(gzip_suffix)]
            break
    if not is_gzipped:
        decompressed_fpath = fpath

    if not any([decompressed_fpath.lower().endswith(fasta_suffix) for fasta_suffix in FASTA_SUFFIXES]):
        raise ValueError(f"Unrecognized suffix on file {fpath}")

    if is_gzipped:
        with gzip.open(fpath, "rt") as fh:
            for record in SeqIO.parse(fh, "fasta"):
                yield record
    else:
        for record in SeqIO.parse(fpath, "fasta"):
            yield record


def parse_gff_file(gff_file) -> List[GffRow]:
    "Parse GFF file into fields from https://en.wikipedia.org/wiki/General_feature_format"
    annotations_fname = os.path.basename(gff_file)

    with file_util.tmp_copy_on_open(gff_file, annotations_fname) as local_file:
        if local_file.endswith(".bgz"):
            # replace suffix so `gunzip` utility works
            gzip_name = local_file[: -len(".bgz")] + ".gz"
            cmd = f"mv {local_file} {gzip_name} && gunzip {gzip_name}"
            _ = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
            local_file = gzip_name[: -len(".gz")]

        with open(local_file) as fh:
            gff_lines = [line.rstrip() for line in fh.readlines()]
        gff_lines = [line.rstrip() for line in gff_lines]
        assert gff_lines[0] == "##gff-version 3"
        gff_lines = gff_lines[1:]
        gff_rows = [GffRow(*line.split("\t")) for line in gff_lines]
    return gff_rows
