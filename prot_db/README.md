#Why BigTable?

This project involved analyzing billions of proteins from public data sources (namely MGnify).
Individual proteins can occur multiple times in our source data.   For example, proteins can occur many times
in one source file, could be conserved across species, or could occur in separate studies of the same microbe. 
We want to deduplicate these so that expensive analyses aren't repeated, and so that we have one canonical collection
of info per protein.

BigTable is a scalable and cost-effective means of dedpulicating many proteins and incrementally
adding per-protein analyses.  If you want to overlook a lot of complexity,
you can think of BigTable as a giant python `dict` keyed on the protein.
Each row contains the column-families described below.

# Design
##Row Key Schema
The row key expresses a 1:1 mapping with amino acid sequences. 

An initial idea was to store the protein sequence directly as the row key.
Protein sequences can vary from ~100 to ~30k amino acids.  4kb is the limit for BigTable rowkeys, so this approach
wasn't valid.

Another idea is to take the protein ID from NCBI, EMBL, or similar.  Our goal is to find novel proteins which,
by definition, should not yet have an ID in these databases.

We ultimately landed on a 64byte row key that's derived from the raw amino acid sequence.

The first 4 bytes are from murmurhash.  Given our scale of data, we might still expect some hash collisions. Therefore,
the next 27 bytes are the first 27 amino acids of the sequence (padded if necessary), 
the following 27 bytes from the end of the sequence (padded if necessary),
and the last 6 bytes are digits of the zero-padded sequence length (e.g. "001254")
"{murmurhash3}{seq[:27]}{seq[-27:]}{len(seq):06d}"

###Row Key Scans
One advantage of this row key design is that it enables table scans over a uniformly random sample of proteins by
leveraging the hash at the beginning of the rowkey.
One can specify a range covering a subset of the possible 4-byte integers / 32-bit numbers, and this will scan a
corresponding fraction of the overall table.

For example, scanning the rowkey range
`[00 00 00 00 ..., 00 D0 00 00 ...)`
is equivalent to scanning `[0, 13631488)` out of 4294967296 32-bit values.  Assuming murmurhash
uniformly shuffles the sequences, this is equivalent to uniformly sampling 0.317% of the bigtable sequences.

This assumption is used in bigtable_constants.get_row_key_boundary() to get uniform slices along binary slices of the
bigtable rowkey space  (e.g. `[0, 1, 2]` out of 2, for num_bits = 1, or `[0, 1, 2, 3, 4, 5, 6, 7, 8] for out of 8,
for num_bits = 3).

##Column Families

###CF_ID_SEQUENCES
This column family contains the raw amino acid sequence

###CF_ID_MGNIFY_STUDY
This column family contains one `MGnify_study_ID, MGnify_analysis_ID` for every MGnify study a protein is found in.
Ideally we would have a column for every analysis the protein appears.
BigTable recommends "billions of rows and thousands of columns".  MGnify has ~4k studies and ~300k analyses.
In order to bound this, we store at most one column per study in each row.


#Getting started with BigTable

Here are some useful resource.

* An emulator to try programs locally / cheaply:
    https://cloud.google.com/bigtable/docs/emulator

* Also dockerized:
    https://github.com/bigtruedata/docker-gcloud-bigtable-emulator

* https://cloud.google.com/bigtable/docs/schema-design

* https://cloud.google.com/bigtable/docs/samples-python-hello

* https://cloud.google.com/bigtable/docs/writing-data#python
