import gzip
from typing import Iterator

from Bio import SeqIO

GZIP_SUFFIXES = [".gz", ".gzip"]
FASTA_SUFFIXES = [".fasta", ".fna", ".ffn", ".faa", ".frn"]


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
