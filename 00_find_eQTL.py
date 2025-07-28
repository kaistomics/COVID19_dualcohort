#!/usr/bin/env python
 
import pandas as pd
from statsmodels.stats.multitest import multipletests

GWAS_summary = pd.read_csv('./Example_data/Immune_response_metaGWAS.tsv', sep='\t', )

close_gene_per_SNP = pd.read_csv(
    './Example_data/hg19_gene_TSS_gene_SNP_position.tsv', sep='\t',)

celltype_list = ['B_Memory', 'CD56bright_NK', 'CD8+_T_Effector_Memory']

total_eQTL_result = []
for celltype in celltype_list:        
    tmp = pd.read_csv('./Example_data/eQTL_results/'+celltype+'_eQTL.tsv', sep='\t', )
    tmp.columns = ['rs_id', 'gene_eQTL', 'beta_eQTL', 't-stat_eQTL', 'p-value_eQTL', 'FDR_eQTL']
    tmp['Cell_Type_eQTL'] = celltype
    
    mytmp = pd.merge(tmp, close_gene_per_SNP, on=['gene_eQTL', 'rs_id'], how='right')
    mytmp = mytmp.loc[~mytmp['beta_eQTL'].isna()]
    mytmp = mytmp[['rs_id', 'gene_eQTL', 'beta_eQTL', 't-stat_eQTL', 'p-value_eQTL', 'FDR_eQTL', 'Cell_Type_eQTL', 'gene_strand_gtf_check', ]].drop_duplicates()    
    mytmp = pd.merge(mytmp, GWAS_summary, on='rs_id', how='left')
    mytmp['p-value_eQTL'] = mytmp['p-value_eQTL'].astype(float)
    total_eQTL_result.append(mytmp)
total_eQTL_result_df = pd.concat(total_eQTL_result)    

expression_ESF = pd.read_csv('./Example_data/Gene_expressing_samples_fraction.tsv', sep='\t')

ESF = 0.5
expression_ESF = expression_ESF.loc[expression_ESF['ESF']>=ESF].copy()

total_eQTL_result_df_orig = pd.merge(total_eQTL_result_df, expression_ESF, on=['gene_eQTL', 'Cell_Type_eQTL'], how='inner')

tmp = []
for celltype in celltype_list:
    cat_eqtl_GWASp_celltype_tmp = total_eQTL_result_df_orig.loc[total_eQTL_result_df_orig['Cell_Type_eQTL'] == celltype].copy()
    cat_eqtl_GWASp_celltype_tmp['FDR_eQTL'] = multipletests(cat_eqtl_GWASp_celltype_tmp['p-value_eQTL'], method='fdr_bh')[1]
    tmp.append(cat_eqtl_GWASp_celltype_tmp)
final_eQTL = pd.concat(tmp)
final_eQTL = final_eQTL.loc[final_eQTL['beta_eQTL']*final_eQTL['BETA_GWAS'] >0].copy()
final_eQTL.to_csv('./Example_data/Immune_response_eQTL.tsv', sep='\t', index=False)

cat_eqtl_GWASp_tmp = total_eQTL_result_df_orig.loc[total_eQTL_result_df_orig['leadSNP'] == True]
celltype_list = list(set(total_eQTL_result_df_orig['Cell_Type_eQTL']))
tmp = []    
for celltype in celltype_list:
    cat_eqtl_GWASp_celltype_tmp = cat_eqtl_GWASp_tmp.loc[cat_eqtl_GWASp_tmp['Cell_Type_eQTL'] == celltype].copy()
    cat_eqtl_GWASp_celltype_tmp = cat_eqtl_GWASp_celltype_tmp.loc[cat_eqtl_GWASp_celltype_tmp['KOREAN_MAF']>=0.1].reset_index(drop=True)
    cat_eqtl_GWASp_celltype_tmp['FDR_eQTL'] = multipletests(cat_eqtl_GWASp_celltype_tmp['p-value_eQTL'], method='fdr_bh')[1]
    tmp.append(cat_eqtl_GWASp_celltype_tmp)
tmp = pd.concat(tmp)
tmp = tmp.loc[tmp['beta_eQTL']*tmp['BETA_GWAS'] >0].copy()

eGene_df = tmp[[ 'Cell_Type_eQTL', 'gene_eQTL',]].drop_duplicates()
eGene_df.columns = [ 'CellType', 'gene',]
eGene_df.to_csv('./Example_data/Immune_response_eGene.tsv', sep='\t', index=False)
