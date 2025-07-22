#!/usr/bin/env python
import os 
import numpy as np
import pandas as pd

from statsmodels.stats.multitest import fdrcorrection
import subprocess
from subprocess import PIPE, Popen
from statsmodels.stats.multitest import multipletests

GWAS_summary = pd.read_csv('./Immune_response_metaGWAS.tsv', sep='\t', )

close_gene_per_SNP = pd.read_csv(
    './hg19_gene_TSS_gene_SNP_position.tsv', sep='\t',)


celltype_list = ['B_Memory',
 'B_Naive',
 'CD14+_Monocyte',
 'CD16+_Monocyte',
 'CD4+_T_Effector_Memory',
 'CD4+_T_Naive_CM',
 'CD56bright_NK',
 'CD56dim_NK',
 'CD8+_T_Effector_Memory',
 'CD8+_T_Naive_CM',
 'Classical_DC',
 'Plasmacytoid_DC',
 'Proliferating_Lymphocyte',
 'gd_T']


total_eQTL_result = []
for celltype in celltype_list:        
    tmp = pd.read_csv('./eQTL_results/'+celltype+'_eQTL.tsv', sep='\t', )
    tmp.columns = ['rs_id', 'gene_eQTL', 'beta_eQTL', 't-stat_eQTL', 'p-value_eQTL', 'FDR_eQTL']
    tmp['Cell_Type_eQTL'] = celltype
    
    mytmp = pd.merge(tmp, close_gene_per_SNP, on=['gene_eQTL', 'rs_id'], how='right')
    mytmp = mytmp.loc[~mytmp['beta_eQTL'].isna()]
    mytmp = mytmp[['rs_id', 'gene_eQTL', 'beta_eQTL', 't-stat_eQTL', 'p-value_eQTL', 'FDR_eQTL', 'Cell_Type_eQTL', 'gene_strand_gtf_check', ]].drop_duplicates()    
    mytmp = pd.merge(mytmp, GWAS_summary, on='rs_id', how='left')
    mytmp['p-value_eQTL'] = mytmp['p-value_eQTL'].astype(float)
    total_eQTL_result.append(mytmp)
total_eQTL_result_df = pd.concat(total_eQTL_result)    

#### filter ESF
expression_ESF = []
for celltype in celltype_list:
    tmp = pd.read_csv('./mean_expression/'+celltype+'_Total_mean_pc.csv', index_col=[0])
    #### Normal expression
    tmp = tmp.T.loc[(tmp.columns.str[5] == 'N')].T
    expression_ESF.append((tmp>0).mean(axis = 1))
expression_ESF = pd.concat(expression_ESF, axis = 1)
expression_ESF.columns = celltype_list
expression_ESF = pd.melt(expression_ESF.reset_index(), id_vars='index', var_name='gene_eQTL', value_name='ESF')
expression_ESF.columns = ['gene_eQTL', 'Cell_Type_eQTL', 'ESF']

ESF = 0.5
expression_ESF = expression_ESF.loc[expression_ESF['ESF']>=ESF].copy()

total_eQTL_result_df_orig = pd.merge(total_eQTL_result_df, expression_ESF, on=['gene_eQTL', 'Cell_Type_eQTL'], how='inner')

tmp = []
for celltype in celltype_list:
    cat_eqtl_GWASp_celltype_tmp = total_eQTL_result_df_orig.loc[total_eQTL_result_df_orig['Cell_Type_eQTL'] == celltype].copy()
    cat_eqtl_GWASp_celltype_tmp['FDR_eQTL'] = fdrcorrection(cat_eqtl_GWASp_celltype_tmp['p-value_eQTL'])[1]
    tmp.append(cat_eqtl_GWASp_celltype_tmp)
final_eQTL = pd.concat(tmp)
final_eQTL = final_eQTL.loc[final_eQTL['beta_eQTL']*final_eQTL['BETA_GWAS'] >0].copy()
final_eQTL.to_csv('./Immune_response_eQTL.tsv', sep='\t', index=False)

cat_eqtl_GWASp_tmp = total_eQTL_result_df_orig.loc[total_eQTL_result_df_orig['leadSNP'] == True]
celltype_list = list(set(total_eQTL_result_df_orig['Cell_Type_eQTL']))
tmp = []    
for celltype in celltype_list:
    cat_eqtl_GWASp_celltype_tmp = cat_eqtl_GWASp_tmp.loc[cat_eqtl_GWASp_tmp['Cell_Type_eQTL'] == celltype].copy()
    cat_eqtl_GWASp_celltype_tmp = cat_eqtl_GWASp_celltype_tmp.loc[cat_eqtl_GWASp_celltype_tmp['KOREAN_MAF']>=0.1].reset_index(drop=True)
    cat_eqtl_GWASp_celltype_tmp['FDR_eQTL'] = fdrcorrection(cat_eqtl_GWASp_celltype_tmp['p-value_eQTL'])[1]
    tmp.append(cat_eqtl_GWASp_celltype_tmp)
tmp = pd.concat(tmp)
tmp = tmp.loc[tmp['beta_eQTL']*tmp['BETA_GWAS'] >0].copy()

eGene_df = tmp[[ 'Cell_Type_eQTL', 'gene_eQTL',]].drop_duplicates()
eGene_df.columns = [ 'celltype', 'gene',]
eGene_df.to_csv('./Immune_response_eGene.tsv', sep='\t', index=False)