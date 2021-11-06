from genegraphdb import *
from genegraphdb import graphdb
from genegraphdb import _load
from genegraphdb import testing
from Bio import SeqIO
import gzip
import re
import sys
import os
from os import remove
from os.path import abspath
import time
from collections import deque

def create_dir(dir_name):
    os.system("mkdir " + dir_name)

def download_dir(bucket_dir):
    bucket_path = "gs://jluo_bucket/" + bucket_dir
    command = "gsutil cp -m " + bucket_path + " ."
    os.system(command)