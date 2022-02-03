import os
import preprocess
import pfamannot
import interproscan
import parse

os.system("mkdir cct cct/OUTPUT hmmer hmmer/out_clean hmmer/out_raw interproscan INPUT INPUT/faa INPUT/fna INPUT/gff OUTPUT")

inactive_tnpBs = ['0648051c710ffc3f99',
  '0da5a5d756496ec587',
  '0f08a20cfe029876e3',
  '124625bd293e1e735d',
  '13d95e27e737ac62a9',
  '15ce911abbf167916e',
  '31428f0f03af2a9999',
  '476fa0561841de3b03',
  '58587edff93c9575ac',
  '6be956b3dd0dfba935',
  '70e9d1aa2ac3c57a62',
  '998f3b369603e3d75b',
  'a9a1793552f754a2c3',
  'b9cc7b8a78f75d3bc5',
  'c328360c354617b7a9',
  'ee5c498bf075a0ad33',
  'f6e74b8eea638f8073']

def update_gffs(ids_list):
    print("creating contig fnas for cctyper")
    create_inputs(ids_list)
    print("annotating pfam domains with cctyper and hmmsearch")
    annotate(ids_list)
    print("processing final gff output")
    final_gff_output(ids_list)
    
def create_inputs(ids_list):
    for input_id in ids_list:
        preprocess.create_gff(input_id)
        preprocess.single_gff_to_faa(input_id)
        preprocess.create_fna(input_id)

def annotate(ids_list):
    #os.system('conda run -n cctyper python cctyper.py')
    #pfamannot.hmmsearch_pool(ids_list)
    #interproscan.interproscan_pool(ids_list)
    for hitid in ids_list:
        pfamannot.format_hmmsearch_output(hitid)
        #interproscan.interproscan(hitid)
        pass
        
def final_gff_output(ids_list):
    for hitid in ids_list:
        parse.put_gff_together(hitid)

update_gffs(inactive_tnpBs)