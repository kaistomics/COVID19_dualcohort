import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
import multiprocessing
import argparse


def _sparse_nanmean(X, axis):
    """
    np.nanmean equivalent for sparse matrices
    """
    if not issparse(X):
        raise TypeError("X must be a sparse matrix")

    ### count the number of nan elements per row/column (dep. on axis)
    Z = X.copy()
    Z.data = np.isnan(Z.data)
    Z.eliminate_zeros()
    n_elements = Z.shape[axis] - Z.sum(axis)

    ### set the nans to 0, so that a normal .sum() works
    Y = X.copy()
    Y.data[np.isnan(Y.data)] = 0
    Y.eliminate_zeros()

    ### the average
    s = Y.sum(axis, dtype='float64') 
    m = s / n_elements

    return m



def geneset_score_per_celltype(adata,
                        annotation_field_name, 
                        celltype,
                        gene_list,
                        gene_pool=None, *args, 
                        ctrl_size = 50, 
                        n_bins = 25,
                        score_name = 'score',
                        random_state = 0,
                        use_raw
                    ):
   
    use_raw = sc._utils._check_use_raw(adata, use_raw)
    if random_state is not None:
        np.random.seed(random_state)

    gene_list_in_var = []
    var_names = adata.raw.var_names if use_raw else adata.var_names
    genes_to_ignore = []
    
    for gene in gene_list:
        if gene in var_names:
            gene_list_in_var.append(gene)
        else:
            genes_to_ignore.append(gene)
    if len(genes_to_ignore) > 0:
        print(f'genes are not in var_names and ignored: {genes_to_ignore}')
    gene_list = set(gene_list_in_var[:])

    if len(gene_list) == 0:
        return None

    if gene_pool is None:                
        gene_pool = list(var_names)
    else:
        gene_pool = [x for x in gene_pool if x in var_names]
    if not gene_pool:
        raise ValueError("No valid genes were passed for reference set.")

    ### Trying here to match the Seurat approach in scoring cells.
    ### Basically we need to compare genes against random genes in a matched
    ### interval of expression.
    _adata = adata.raw if use_raw else adata


    celltype_indices = _adata.obs.loc[_adata.obs[annotation_field_name] == celltype].index

    _adata_subset = (
        _adata[celltype_indices, gene_pool] if len(gene_pool) < len(_adata.var_names) else _adata
    )
    
#    print('calculate average value of genes')
    if issparse(_adata_subset.X):
        obs_avg = pd.Series(                                          ### gene mean expression, all cells. 
            np.array(_sparse_nanmean(_adata_subset.X, axis=0)).flatten(),
            index=gene_pool,
        )  ### average expression of genes
    else:
        obs_avg = pd.Series(
            np.nanmean(_adata_subset.X, axis=0), index=gene_pool
        )  ### average expression of genes

    obs_avg = obs_avg[
        np.isfinite(obs_avg)
    ]  ### Sometimes (and I don't know how) missing data may be there, with nansfor
    
#    print('start control gene selection')
    n_items = int(np.round(len(obs_avg) / (n_bins - 1)))            ### number of gene rank groups. 
    obs_cut = obs_avg.rank(method='min') // n_items               ###  유전자들의 rank를 n_items로 나눈 몫. 즉, 유전자의 계급값. 
    control_genes = set()

    # now pick `ctrl_size` genes from every cut
    for cut in np.unique(obs_cut.loc[list(gene_list)]):       ### for 내가 관심있는 유전자들의 계급값들 
        r_genes = np.array(obs_cut[obs_cut == cut].index)       ### 그 계급값에 해당되는 모든 유전자
        np.random.shuffle(r_genes)
        # uses full r_genes if ctrl_size > len(r_genes)
        control_genes.update(set(r_genes[:ctrl_size]))          ### 그 계급값에 해당되는 모든 유전자 중에서 50개만 선정함. & for문 반복. 

    # To index, we need a list – indexing implies an order.
    control_genes = list(control_genes - gene_list)
    gene_list = list(gene_list)
#    print('finished selecting control genes')    
    
#    print('calculate subtraction')
    

    X_list = _adata[celltype_indices, gene_list].X
    if issparse(X_list):
        X_list = np.array(_sparse_nanmean(X_list, axis=1)).flatten()         #관심있는 유전자의 total cell expression. 각 세포 별로 계산.           
    else:
        X_list = np.nanmean(X_list, axis=1, dtype='float64')

    X_control = _adata[celltype_indices, control_genes].X                                #control 유전자의 total cell expression. 
    if issparse(X_control):
        X_control = np.array(_sparse_nanmean(X_control, axis=1)).flatten()
    else:
        X_control = np.nanmean(X_control, axis=1, dtype='float64')
    
    score = X_list - X_control              #두 값을 뺀 것. 
    tmp_series = pd.Series(score, index = celltype_indices)
    output_df = pd.DataFrame(tmp_series)
    output_df['celltype'] = celltype
    return output_df


def main():
    parser = argparse.ArgumentParser(description='Calculate geneset scores per cell type')
    parser.add_argument('--adata_file', type=str, required=True,
                       help='Path to the h5ad file containing the AnnData object, ex) ./Example_data/Adata_example.h5ad')
    parser.add_argument('--filename', type=str, required=True,
                       help='Path to the TSV file containing eGene data, ./Example_data/Immune_response_eGene.tsv')
    parser.add_argument('--n_process', type=int, required=True,
                       help='Number of processes to use for multiprocessing')
    parser.add_argument('--out_name', type=str, required=True,
                       help='Output filename for the results CSV')
    parser.add_argument('--annotation_name', type=str, required=True,
                       help='Name of the annotation field in adata.obs, ex) Final_Annotation')
    args = parser.parse_args()
    
    adata_file = args.adata_file
    filename = args.filename
    n_process = args.n_process
    out_name = args.out_name
    annotation_name = args.annotation_name

    ### Load data
    adata = sc.read_h5ad(adata_file)
    adata.raw.obs = adata.obs

    celltypes_in_adata = ['B Memory', 'CD56bright NK', 'CD8+ T Effector/Memory']
    celltypes_in_significant_df = ['B_Memory', 'CD56bright_NK', 'CD8+_T_Effector_Memory']

    def calculate_for_each_celltype(i):
        celltype_adata = celltypes_in_adata[i] # cell type name in adata
        celltype_sig = celltypes_in_significant_df[i] # cell type namd in eGene df
        celltype_posi_gene = list(eGene_df.loc[eGene_df['CellType'] == celltype_sig]['gene'])
        celltype_score = geneset_score_per_celltype(adata, annotation_name, celltype_adata, gene_list=celltype_posi_gene, 
                                        use_raw=True)
        return celltype_score

    eGene_df = pd.read_csv(filename, sep='\t')
    with multiprocessing.Pool(n_process) as pool:
        cellscores = pool.map(calculate_for_each_celltype, range(len(celltypes_in_adata)))

    cellscores_concat_df = pd.concat(cellscores, axis=0)
    cellscores_concat_df.columns = ['PTS', 'celltype']
    cellscores_concat_df.to_csv(out_name)
    
    print(f"Results saved to {out_name}")


if __name__ == "__main__":
    main()
    
    
