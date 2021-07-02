from genegraphdb import *
from genegraphdb import graphdb
from genegraphdb import _load
from Bio import SeqIO
import gzip
import re
import sys
import os
from os import remove
from os.path import abspath
import time
from collections import deque, ChainMap
from csv import reader

def load_CRISPR():
    pass
def merge_gff():
    # parse through minced.gff
    # cat 8156401/8156401.minced.gff | grep ID=CRISPR > temp.minced.gff
    # cat temp.prodigal.gff temp.minced.gff > temp2.merged.gff
    # sortBed -i temp.merged.gff > temp.merged.sorted.gff
    return merged_gff