#!/usr/bin/env python
# coding: utf-8

import pandas as pd


i = 0
subject = 'Noninfected_741'
TRB_expansion_max = []
for timepoint in [2, 4]:

    div = pd.read_csv('./Example_data/input_data/'+
                           '{}_T{}0_TRB.csv'.format(subject, timepoint))
    div['expansion'] = (div['0'].apply(lambda x: x if x <1e+8 else x/(1e+8)))
    div = div.loc[div['expansion'] != 0]
    TRB_expansion_max.append((subject, 'T{}'.format(timepoint),div['expansion'].max() ))
        
TRB_expansion_max_df = pd.DataFrame(TRB_expansion_max)
TRB_expansion_max_df.columns = ['subject', 'TD', 'values']
TRB_expansion_max_df = TRB_expansion_max_df.pivot(index = 'subject', columns='TD', values='values')
TRB_expansion_max_df.to_csv('./Example_data/TRB_Expansion_max.csv')



i = 0
subject = 'Noninfected_741'
IGH_expansion_max_df = []
for timepoint in [2, 4]:
    div = pd.read_csv('./Example_data/input_data/'+
                           '{}_T{}0_IGH.csv'.format(subject, timepoint))
    div['expansion'] = (div['0'].apply(lambda x: x if x <1e+8 else x/(1e+8)))           
    IGH_expansion_max_df.append((subject, 'T{}'.format(timepoint),div['expansion'].max() ))

IGH_expansion_max_df = pd.DataFrame(IGH_expansion_max_df)
IGH_expansion_max_df.columns = ['subject', 'TD', 'values']
IGH_expansion_max_df = IGH_expansion_max_df.pivot(index = 'subject', columns='TD', values='values')
IGH_expansion_max_df.to_csv('./Example_data/IGH_Expansion_max.csv')

