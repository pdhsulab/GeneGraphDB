#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
def update_modules():
    for module in os.listdir():
        if os.path.isfile(module) and ".ipynb" in module:
            cmd = "jupyter nbconvert --to python " + module
            os.system(cmd)
    print("finished updating modules")
update_modules()


# In[ ]:


os.system("""mkdir cct cct/OUTPUT hmmer hmmer/out_clean hmmer/out_raw interproscan INPUT INPUT/faa INPUT/fna INPUT/gff OUTPUT minced""")


# In[2]:


import preprocess
import pfamannot
import parse
import minced
import os


# ### for each .examples.gff file (target p30 neighborhood)
# - get multifasta file from gff: input for hmmsearch and cctyper
# - parse hmmsearch output
# - parse cctyper outputs
# 
# ### repeat for each .examples.gff file in input directory
# - 
# 

# In[14]:


# load input files into input directory
# hits_near_inactive_tnpBs = ['7205dff302ff900300', '6c3b2ade31d3833745', '23422406293d40c201',
#        'c1050b21cc75640d51', 'aa6af7e9289c3558d3', 'd2246f26c16fb9eec9',
#        '6bf6c4c7da68779d7a']

# hits_near_inactive_tnpBs = ['0093c124a5b6ee038d',
#  '23422406293d40c201',
#  '27177e034d397f4d56',
#  '3a86b6b12ef4de1264',
#  '49c65ceeadb463ef7d',
#  '4a23bcab7c50591a6a',
#  '4e45823b6b83464d6e',
#  '6bf6c4c7da68779d7a',
#  '7205dff302ff900300',
#  '986ad2cde69643b021',
#  '9a8a6e441d6dea8617',
#  'aa6af7e9289c3558d3',
#  'ba436543bec6c30920',
#  'c1050b21cc75640d51',
#  'c15257c28ad9a20cdd',
#  'c7a1b8f01c5b7c48a2',
#  'd2246f26c16fb9eec9',
#  'd2248a31d0e133ef11',
#  'e8a14f064cce46f350']

# hits_near_inactive_tnpBs = ['f45a906808c45716e4',
#  'ee669bb85d7c4a9459', '031933be88d0f6dd71',
#  'cc3b1023fa38a3c1e6', '0204aff763d8b45bfe',
#  'e95a232b429ceb001b', '9038f7993863bcadca',
#  '6bf6c4c7da68779d7a', '0b921ace501a0fe980',
#  'd862301d73972ea41c', '2a89c5bbdcfcb0cd2a',
#  '4dcbca21fe7e9a5c29', 'e8e35986a4c05b4442',
#  '4e45823b6b83464d6e', '66e3c83c192b25b489',
#  'e8a14f064cce46f350', '2bce6e5989f64eecdf',
#  '23422406293d40c201', 'c1050b21cc75640d51',
#  '04bbfb58c48a0e8a29', '32de67ee40e11f3f2f',
#  'e36a64acf9381581d3', '88a7c58044ed1c0c01',
#  'c36b3b94eaf89826ba', '99f6125a673bcd700d',
#  '260273e052f76dffd2', '62fe559eac786c1400',
#  '3b91e97121e56aa047',
#  'c7de21b826534e1702', '7684dc2cf4ee624258',
#  'aa6af7e9289c3558d3', 'c7a1b8f01c5b7c48a2',
#  '0a3217e07a25d0cbe0', '1b0744e098d12904a7',
#  'f3a73b3acc03da9bf8', '56bf9d999b00b23351',
#  '214319fd8cc31efee9', '7c8d22120719f2db6d',
#  'a366ac26802335850e', '49c65ceeadb463ef7d',
#  '97127aa4b059e7df3b', 'a0d79afeb20bd32508',
#  '729ad68f76fcc9ad44', 'ad17062fc0b1953724',
#  '7e848117da228c93cd', 'd90f25484d732aa877']

hits_near_inactive_tnpBs = ['9038f7993863bcadca', '6bf6c4c7da68779d7a', 'd862301d73972ea41c', '2a89c5bbdcfcb0cd2a', 
                        '23422406293d40c201', 'c1050b21cc75640d51', '04bbfb58c48a0e8a29', '260273e052f76dffd2',
                        'c7de21b826534e1702', 'aa6af7e9289c3558d3', 'c7a1b8f01c5b7c48a2', '0a3217e07a25d0cbe0',
                        '214319fd8cc31efee9', '7c8d22120719f2db6d', 'a0d79afeb20bd32508', 'd90f25484d732aa877']
hits_not_fetched = []

input_files = [hitid + ".examples.gff" for hitid in hits_near_inactive_tnpBs]
input_filepaths = ["../ggdbfetch_output/" + hitpath for hitpath in input_files]
for i in range(len(input_filepaths)):
    file_path = input_filepaths[i]
    cmd = "cp " + file_path + " INPUT/gff/" + input_files[i]
    if os.path.isfile(file_path):
        os.system(cmd)
    else:
        hits_not_fetched.append(hits_near_inactive_tnpBs[i])
# to do - fetch these target gene neighbourhood close to semi-inactive tnpBs
for failed_hit in hits_not_fetched:
    hits_near_inactive_tnpBs.remove(failed_hit)


# In[16]:


infile_paths = ['INPUT/gff/' + file for file in os.listdir('INPUT/gff') if file != '.ipynb_checkpoints']


# In[23]:


def update_gffs(ids_list):
    print("creating contig fnas for cctyper")
    #create_inputs(infile_paths)
    print("annotating pfam domains with cctyper and hmmsearch")
    #annotate(ids_list)
    print("processing final gff output")
    final_gff_output(ids_list)
    
def create_inputs(infile_paths):
    for path in infile_paths:
        preprocess.single_gff_to_fna(path)
        preprocess.single_gff_to_faa(path)

def annotate(ids_list):
    os.system('conda run -n cctyper python cctyper.py')
    pfamannot.hmmsearch_pool(ids_list)
    minced.annot_crispr_pool(ids_list)
    for hitid in ids_list:
        pfamannot.format_hmmsearch_output(hitid)
        
def final_gff_output(ids_list):
    for hitid in ids_list:
        parse.put_gff_together(hitid) 
    #parse.gff_output_csv()


# In[24]:


update_gffs(hits_near_inactive_tnpBs)


# In[ ]:




