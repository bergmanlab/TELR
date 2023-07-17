#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import json
import glob


# In[2]:


# load data
eval_dir = "/scratch/sh60271/s2rplus/eval_telr_0830/out"
pattern = "/**/eval_liftover.json"
eval_json_files = glob.glob(eval_dir + pattern, recursive=True)

# load data
eval_data = []
for eval_json_file in eval_json_files:
    with open(eval_json_file) as f:
        data = json.load(f)
    eval_data.append(data)

# convert to pandas
df = pd.DataFrame(eval_data)

# write in csv format
output_file = "/scratch/sh60271/s2rplus/eval_telr_0830/out/telr_eval_coord.csv"
df.to_csv(
    output_file,
    sep=",",
    index=False,
    header=True,
)


# In[3]:


# load data
eval_dir = "/scratch/sh60271/s2rplus/eval_telr_0830/out"
pattern = "/**/*diploid*/telr_af_eval/telr_eval_af.json"
eval_json_files = glob.glob(eval_dir + pattern, recursive=True)

# load data
eval_data = []
for eval_json_file in eval_json_files:
    with open(eval_json_file) as f:
        data = json.load(f)
    eval_data.append(data)

# convert to pandas
df = pd.DataFrame(eval_data)

# write in csv format
output_file = "/scratch/sh60271/s2rplus/eval_telr_0830/out/telr_eval_diploid_af.csv"
df.to_csv(
    output_file,
    sep=",",
    index=False,
    header=True,
)


# In[4]:


# load data
eval_dir = "/scratch/sh60271/s2rplus/eval_telr_0830/out"
pattern = "/**/*tetra*/telr_af_eval/telr_eval_af.json"
eval_json_files = glob.glob(eval_dir + pattern, recursive=True)

# load data
eval_data = []
for eval_json_file in eval_json_files:
    with open(eval_json_file) as f:
        data = json.load(f)
    eval_data.append(data)

# convert to pandas
df = pd.DataFrame(eval_data)

# write in csv format
output_file = "/scratch/sh60271/s2rplus/eval_telr_0830/out/telr_eval_tetraploid_af.csv"
df.to_csv(
    output_file,
    sep=",",
    index=False,
    header=True,
)


# In[3]:


# load data
eval_dir = "/scratch/sh60271/s2rplus/eval_telr_0830/out"
pattern = "/**/*diploid*/telr_cn_eval/telr_eval_te_cn.json"
eval_json_files = glob.glob(eval_dir + pattern, recursive=True)

# load data
eval_data = []
for eval_json_file in eval_json_files:
    with open(eval_json_file) as f:
        data = json.load(f)
    eval_data.append(data)

# convert to pandas
df = pd.DataFrame(eval_data)

# write in csv format
output_file = "/scratch/sh60271/s2rplus/eval_telr_0830/out/telr_eval_diploid_cn.csv"
df.to_csv(
    output_file,
    sep=",",
    index=False,
    header=True,
)


# In[4]:


# load data
eval_dir = "/scratch/sh60271/s2rplus/eval_telr_0830/out"
pattern = "/**/*tetra*/telr_cn_eval/telr_eval_te_cn.json"
eval_json_files = glob.glob(eval_dir + pattern, recursive=True)

# load data
eval_data = []
for eval_json_file in eval_json_files:
    with open(eval_json_file) as f:
        data = json.load(f)
    eval_data.append(data)

# convert to pandas
df = pd.DataFrame(eval_data)

# write in csv format
output_file = "/scratch/sh60271/s2rplus/eval_telr_0830/out/telr_eval_tetraploid_cn.csv"
df.to_csv(
    output_file,
    sep=",",
    index=False,
    header=True,
)


# In[5]:


eval_dir = "/scratch/sh60271/s2rplus/eval_telr_0830/out"
pattern = "/**/telr_seq_eval/seq_eval.json"
eval_json_files = glob.glob(eval_dir + pattern, recursive=True)
# print(eval_json_files)

appended_eval_seq_data = []
for eval_json_file in eval_json_files:
    setting = eval_json_file.split("/")[6]
    with open(eval_json_file) as f:
        eval_seq_data = json.load(f)
        df = pd.DataFrame(eval_seq_data)
        df["setting"] = setting
        appended_eval_seq_data.append(df)
appended_eval_seq_data = pd.concat(appended_eval_seq_data)

# write in csv format
output_file = "/scratch/sh60271/s2rplus/eval_telr_0830/out/telr_eval_seqs.csv"
appended_eval_seq_data.to_csv(
    output_file,
    sep=",",
    index=False,
    header=True,
)


# In[ ]:




