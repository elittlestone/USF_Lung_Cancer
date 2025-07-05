import argparse
import pandas as pd
import pyreadr
from scipy.io import mmread
import anndata as ad
import scanpy as sc
import os 
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

def main():

    # Directory where figures generated will be outputted 
    output_directory = "results/figures"
    os.makedirs(output_directory, exist_ok=True)
    sc.settings.figdir = output_directory

    parser = argparse.ArgumentParser(description = "Txrn1 expression comparison script")

    parser.add_argument("--expression_matrix_file", type = str, required = True, 
    help = "Sparse expression matrix file")

    parser.add_argument("--genes", type = str, required = True, 
    help = "Genes extracted from matrix file")

    parser.add_argument("--barcodes", type = str, required = True, 
    help = "Barcodes extracted from matrix file")

    parser.add_argument("--annotation_file", type = str, required = True, 
    help = "Annotated celltypes")

    args = parser.parse_args()
    matrix_file, genes, barcodes, annotation_file = args.expression_matrix_file, args.genes, args.barcodes, args.annotation_file

    read_exp_matrix(matrix_file, genes, barcodes, annotation_file)

def read_exp_matrix(matrix_file, genes, barcodes, annotation_file):
    # Load matrix
    matrix = mmread(matrix_file).tocsr().T

    # Load gene names and cell IDs
    genes = pd.read_table(genes, header = None)[0].astype(str).tolist()
    barcodes = pd.read_table(barcodes, header = None)[0].astype(str).tolist()


    # Create annotation data object 
    adata = sc.AnnData(matrix)
    adata.var_names = genes
    adata.obs_names = barcodes

    # Load cell-type annotations
    annotations_dict = pyreadr.read_r(annotation_file)
    annotations_df = annotations_dict[None]
    annotations_df = annotations_df.set_index("cell")

    # Rename the annotation column to "cell_type"
    annotations_df = annotations_df.rename(columns={"truth": "cell_type"})

    # Merge the annotations table with the anndata object 
    adata.obs = adata.obs.join(annotations_df, how="left")

    # Plot to compare expression across cell-types (cancer cells vs immune cells)
    sc.pl.violin(adata, keys="TXNRD1", groupby="cell_type", show = False, stripplot = True)
    plt.gcf().set_size_inches(12, 6)
    plt.xticks(rotation = 45, ha = 'right')  # rotate labels for better readability
    plt.savefig("results/figures/violin_cancer_vs_immune_cells.png")
    #plt.show()

    sc.pl.dotplot(adata, "TXNRD1", groupby="cell_type", save = "_cancer_vs_immune.png")


    cancer_expr = adata[adata.obs["cell_type"] == "Cancer cells", "TXNRD1"].X
    immune_expr = adata[adata.obs["cell_type"] != "", "TXNRD1"].X

    cancer_expr = cancer_expr.toarray().flatten() if hasattr(cancer_expr, "toarray") else cancer_expr.flatten()
    immune_expr = immune_expr.toarray().flatten() if hasattr(immune_expr, "toarray") else immune_expr.flatten()

    stat, pval = mannwhitneyu(cancer_expr, immune_expr, alternative="two-sided")
    print(f"Wilcoxon test p-value: {pval}")


    # Perform PCA, NN, UMAP
    #sc.tl.pca(adata)
    #sc.pp.neighbors(adata)
    #sc.tl.umap(adata)
    #sc.pl.umap(adata, color="TXNRD1", save="_TXNRD1_expression.png")



    

    #mean_exp = adata.obs.groupby("cell_type")["TXNRD1"].mean()
    #mean_exp.plot(kind='bar')
    #plt.ylabel("Mean TXNRD1 expression")
    #plt.savefig("mean_exp_barplot.png")
    #plt.show()



if __name__ == "__main__":
    main()