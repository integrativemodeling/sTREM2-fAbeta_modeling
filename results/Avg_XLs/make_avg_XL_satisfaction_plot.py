#!/usr/bin/env python3

import pandas as pd
import glob
import sys
import os
import shutil
import random
import matplotlib.pyplot as plt
from itertools import islice
from itertools import chain
import re


xl_csv = glob.glob('../../analysis/XLs_dist_info_*.csv') # please provide correct path of your analysis directory
sorted_xl_csv =sorted(xl_csv)
a_dict = {}
common_col = 'MC_frame'
num_xl_type = 15
num_chain = 45
cwd = os.getcwd()


def extract_shortest_xl(xl_csv):

    df = pd.read_csv(xl_csv)
    col_list = list(df.columns)
    frame_num = df[common_col].to_list()
    a_dict[common_col] = frame_num
    # How many elements each list should have
    n = 45
    x = [col_list[i:i + n] for i in range(2, len(col_list), n)] 
    counter = 1
    final_col_name = []
    for f in x:
        f.insert(0, common_col)
        new_df = df[f]
        a_name =  '_'.join(f[1].strip().split('|')[2:7])
        final_col_name.append(a_name)
        a_dict[a_name] = new_df.min(axis=1).tolist()
    df_final = pd.DataFrame.from_dict(a_dict)

    return df_final, final_col_name

fig_name = '{}.pdf'
df_name = 'df_final_{}'

### combine all dataframes
df_list = []
for i in range(len(sorted_xl_csv)):
    df_dummy = df_name.format(str(i+1)) 
    print(df_dummy)
    df_dummy, final_col_name = extract_shortest_xl(sorted_xl_csv[i])
    df_list.append(df_dummy)

df_final = pd.concat(df_list)
df_col_names = df_final.columns
col_name_dic = {}
xls = []
for f in df_col_names:
    if f == 'MC_frame':
        pass
#        col_name_dic[f] = f
    else:
        f_list = f.split('_')[1:]
        print(f_list)
        res = " ".join(re.findall("[a-zA-Z]+", f_list[0]))
        f_new = ''.join(['A', f_list[1], 'T', f_list[3]])
        col_name_dic[f] = f_new
print(col_name_dic)

# finding duplicate values from a dictionary using set
rev_dict = {}
for key, value in col_name_dic.items():
    rev_dict.setdefault(value, set()).add(key)
  
  
result = set(chain.from_iterable(
         values for key, values in rev_dict.items()
         if len(values) > 1))
  
# randomly pop one from the duplicate

del_key = list(result).pop()
del col_name_dic[del_key]

df_for_avg = df_final.loc[:, df_final.columns != 'MC_frame']
df_ = df_for_avg.rename(col_name_dic, axis=1)
df_.pop(del_key)
df_.mean().plot.barh(fontsize=11)
plt.axvline(x=30, color='grey',ls='dashed',lw=3)
xlabel = 'Avg. XL distance (in ' + r'$\AA$' + ')'
plt.xlabel(xlabel, fontweight='bold')
plt.savefig('avg_dist.pdf')
plt.savefig('avg_dist.eps')


