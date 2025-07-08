import argparse
import pandas as pd
import pyreadr
from scipy.io import mmread
import anndata as ad
import scanpy as sc
import os 
import matplotlib as mtplt
mtplt.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from tqdm import tqdm

def main():

    # Directory where figures generated will be outputted 
    output_directory = "results/figures"
    os.makedirs(output_directory, exist_ok=True)
    sc.settings.figdir = output_directory
    sc.settings.verbosity = 'debug'


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

    # Read in expression matrix and genes, barcodes
    adata = read_exp_matrix(matrix_file, genes, barcodes)
    adata = plots_and_mannwhitney(adata, annotation_file)
    post_processing_and_clustering(adata)


def read_exp_matrix(matrix_file, genes, barcodes):
    # Load matrix
    matrix = mmread(matrix_file).tocsr().T

    # Load gene names and cell IDs
    genes = pd.read_table(genes, header = None)[0].astype(str).tolist()
    barcodes = pd.read_table(barcodes, header = None)[0].astype(str).tolist()


    # Create annotation data object 
    adata = sc.AnnData(matrix)
    adata.var_names = genes
    adata.obs_names = barcodes

    # Normalize expression
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    # Scale 
    sc.pp.scale(adata)

    return adata

def plots_and_mannwhitney(adata, annotation_file):

    # Load cell-type annotations
    annotations_dict = pyreadr.read_r(annotation_file)
    annotations_df = annotations_dict[None]
    annotations_df = annotations_df.set_index("cell")

    # Rename the annotation column to "cell_type"
    annotations_df = annotations_df.rename(columns={"truth": "cell_type"})

    # Get counts for each cell type

    celltype_counts = annotations_df["cell_type"].value_counts()
    print(f"Cell type counts: {celltype_counts}")

    # Merge the annotations table with the anndata object 
    adata.obs = adata.obs.join(annotations_df, how="left")

    # Plot to compare expression across cell-types (cancer cells vs immune cells)
    genes_of_interest = ["TXNRD1", "TXN"]

    # Initalize results list to then output to CSV 
    results = []
    reference_group = "Cancer cells"

    #matches = [var for var in adata.var_names if "txn1" in var.lower()]
    for gene in tqdm(genes_of_interest):
        # First create summary plots, then perform statistical analyses

        # Violin plot
        sc.pl.violin(adata, keys=gene, groupby="cell_type", show = False, stripplot = True)
        plt.gcf().set_size_inches(12, 6)
        plt.xticks(rotation = 45, ha = 'right')  # rotate labels for better readability
        plt.savefig(f"results/figures/violin_{gene}_cancer_vs_immune_cells.png")
        plt.close("all")

        # Dotplot
        sc.pl.dotplot(adata, gene, groupby="cell_type", show = False)
        plt.savefig(f"results/figures/dotplot_{gene}_cancer_vs_immune_cells.png")
        plt.close("all")

        # Statistical tests 
        cancer_expr = adata[adata.obs["cell_type"] == "Cancer cells", gene].X
        cancer_expr = cancer_expr.toarray().flatten()
        cancer_mean = cancer_expr.mean()

        # Perform Mann-Whitney U for each cell type against cancer cells
        for cell_type in adata.obs["cell_type"].unique():
            if cell_type == reference_group:
                expression = cancer_expr
                mean_expression = cancer_mean
                stat, pval = None, None
            else:
                expression = adata[adata.obs["cell_type"] == f"{cell_type}", gene].X
                expression = expression.toarray().flatten()
                mean_expression = expression.mean()
                stat, pval = mannwhitneyu(cancer_expr, expression, alternative = "two-sided")
   
            
            # Store results for each cell-type and each gene in a dictionary
            results.append({
                "Gene": gene,
                "Cell Type": cell_type,
                "Mean Expression": mean_expression,
                "P-value": pval,
            })

    # Output results/statistics to CSV 
    results_df = pd.DataFrame(results)
    results_df.sort_values(by = ["Cell Type", "Gene"], inplace = True)
    results_df.to_csv("results/mannwhitney_celltype_comparison.csv", index = False)

    return adata


def post_processing_and_clustering(adata):

    # Perform PCA, kNN, UMAP
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    # Plot UMAP
    genes_of_interest = ["TXNRD1", "TXN"]
    for gene in tqdm(genes_of_interest):
        sc.pl.umap(adata, color="cell_type", show=False, save = f"_{gene}_celltypes.pdf")
        sc.pl.umap(adata, color = gene, show=False, cmap = "inferno", vmin = 0, size = 8, vmax = "p99.5", save = f"_{gene}_expression.pdf")

if __name__ == "__main__":
    main()