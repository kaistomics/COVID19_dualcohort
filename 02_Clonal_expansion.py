#!/usr/bin/env python
# coding: utf-8



import pandas as pd


with open('subject_list.txt', 'r') as f:
    subjects = f.read().splitlines()


TRB_expansion_max = []

i = 0
for subject in subjects:
    for timepoint in [ 3,  5]:
        try:
            div = pd.read_csv('./TRB_expansion/'+
                                   'TRB_scBCR_{}_T{}1_VJCDR3clone_MAITremoved240328_PostNAunremoved.csv'.format(subject, timepoint))
            div['expansion'] = (div['0'].apply(lambda x: x if x <1e+8 else x/(1e+8)))
            div = div.loc[div['expansion'] != 0]
            
            TRB_expansion_max.append((subject, '{}/1'.format(timepoint),div['expansion'].max() ))
        except:
            pass
        
TRB_expansion_max_df = pd.DataFrame(TRB_expansion_max)
TRB_expansion_max_df.columns = ['subject', 'TD', 'values']
TRB_expansion_max_df = TRB_expansion_max_df.pivot(index = 'subject', columns='TD', values='values')
TRB_expansion_max_df.columns = ['T2', 'T4']
TRB_expansion_max_df.to_csv('./TRB_Expansion_max.csv')




IGH_expansion_max_df = []

i = 0
for subject in subjects:
    for timepoint in [3, 5]:
        try:
            div = pd.read_csv('./IGH_expansion/'+
                                   'IGH_scBCR_{}_T{}1_Lev02clone_240328_PostNAunremoved.csv'.format(subject, timepoint))
            div['expansion'] = (div['0'].apply(lambda x: x if x <1e+8 else x/(1e+8)))           
            IGH_expansion_max_df.append((subject, '{}/1'.format(timepoint),div['expansion'].max() ))
        except:
            pass
 
        
IGH_expansion_max_df = pd.DataFrame(IGH_expansion_max_df)
IGH_expansion_max_df.columns = ['subject', 'TD', 'values']
IGH_expansion_max_df = IGH_expansion_max_df.pivot(index = 'subject', columns='TD', values='values')
IGH_expansion_max_df.columns = ['T2', 'T4']
IGH_expansion_max_df.to_csv('./IGH_Expansion_max.csv')

