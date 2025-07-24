import argparse
import pandas as pd
import pyreadr
from scipy.io import mmread
import scanpy as sc
import os 
import matplotlib as mtplt
mtplt.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import seaborn as sns
from tqdm import tqdm

def main():

    # Directory where figures generated will be outputted 
    output_directory = "results/figures"
    os.makedirs(output_directory, exist_ok=True)
    sc.settings.figdir = output_directory
    sc.settings.verbosity = 'debug'

    # Command line arguments 
    parser = argparse.ArgumentParser(description = "Txrn1 expression comparison script")

    parser.add_argument("--expression_matrix_file", type = str, required = True, 
    help = "Sparse expression matrix file")

    parser.add_argument("--genes", type = str, required = True, 
    help = "Genes extracted from matrix file")

    parser.add_argument("--barcodes", type = str, required = True, 
    help = "Barcodes extracted from matrix file")

    parser.add_argument("--annotation_file", type = str, required = True, 
    help = "Annotated celltypes")

    parser.add_argument("--gene_list", type = str, required = True, 
    help = "List of genes for analysis")
   
    parser.add_argument("--csv_output", type = str, required = True,
                        help = "Output csv for statistical tests")
    args = parser.parse_args()

    # Read in genes to analyze from excel file 
    dna_repair_genes = pd.read_excel(args.gene_list)
    genes_of_interest = [gene for gene in dna_repair_genes["geneid"]]

    
    matrix_file, genes, barcodes, annotation_file, csv_output = (args.expression_matrix_file, args.genes, 
                                                                 args.barcodes, args.annotation_file, args.csv_output)

    # Read in expression matrix and genes, barcodes
    adata = read_exp_matrix(matrix_file, genes, barcodes)

    # Create summary plots
    adata = summary_plots(adata, annotation_file, genes_of_interest)

    # Perform DGEA
    statistical_tests(adata, genes_of_interest, csv_output)

    # Clustering and UMAP 
    post_processing_and_clustering(adata, genes_of_interest)


def read_exp_matrix(matrix_file, genes, barcodes):
    """
    Function for reading in the expression matrix from an R object (.rds files)
    @matrix_file : gene expression matrix
    @genes : gene names obtained from cell-annotation file 
    @barcodes: barcodes/cell_ids obtained from cell-annotation file
    """
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

    # Log transformation
    sc.pp.log1p(adata)

    # Store raw counts before scaling 
    adata.raw = adata.copy()

    # Scale 
    sc.pp.scale(adata)

    return adata

def summary_plots(adata, annotation_file, genes_of_interest):

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

    for gene in tqdm(genes_of_interest):

        # Violin plot
        sc.pl.violin(adata, keys=gene, groupby="cell_type", show = False, stripplot = True)
        plt.gcf().set_size_inches(12, 6)
        plt.xticks(rotation = 45, ha = 'right')  # rotate labels for better readability
        plt.savefig(f"results/figures/violin_{gene}_cancer_vs_immune_cells.png", bbox_inches='tight', dpi=300)
        plt.close("all")

        # Dotplot
        sc.pl.dotplot(adata, gene, groupby="cell_type", show = False)
        plt.savefig(f"results/figures/dotplot_{gene}_cancer_vs_immune_cells.png", bbox_inches='tight', dpi=300)
        plt.close("all")
    
    return adata
   
def statistical_tests(adata, genes_of_interest, csv_output):
    """
    Run Mannwhitneyu tests to test for difference in expression of given genes in cancer vs other immune cells
    @adata : Annotation data object (gene expression matrix + other metadata)
    @genes_of_interest : genes to be tested 
    """

    # Initalize results list to then output to CSV 
    results = []
    reference_group = "Cancer cells"

    for gene in genes_of_interest:

        # Statistical tests 
        cancer_expr = adata.raw[adata.obs["cell_type"] == "Cancer cells", gene].X
        cancer_expr = cancer_expr.toarray().flatten()
        cancer_mean = cancer_expr.mean()

        # Perform Mann-Whitney U for each cell type against cancer cells
        for cell_type in adata.obs["cell_type"].unique():
            if cell_type == reference_group:
                expression = cancer_expr
                mean_expression = cancer_mean
                stat, pval = None, None
            else:
                expression = adata.raw[adata.obs["cell_type"] == f"{cell_type}", gene].X
                expression = expression.toarray().flatten()
                mean_expression = expression.mean()
                stat, pval = mannwhitneyu(cancer_expr, expression, alternative = "greater")
   
            
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
    results_df.to_csv(csv_output, index = False)

    return adata


def post_processing_and_clustering(adata, genes_of_interest):

    # Perform PCA, kNN, UMAP
    sc.tl.pca(adata)
    sc.pl.pca_variance_ratio(adata, n_pcs = 50, log = True, save = "_.png")
    sc.pp.neighbors(adata, n_pcs = 20)
    sc.tl.leiden(adata)
    sc.tl.umap(adata)

    # Plot UMAPs
    
    # Create custom palette for easier differentiation of celltypes 
    n = len(adata.obs["cell_type"].unique())
    palette = sns.color_palette("tab20", n)

    # UMAP of celltypes 
    sc.pl.umap(adata, color = "cell_type", show = False, palette = palette, save = "_celltypes.png")
    
    for gene in tqdm(genes_of_interest):
        sc.pl.umap(adata, color = gene, show=False, cmap = "inferno", vmin = 0, size = 8, vmax = "p99.5", save = f"_{gene}_expression.png")
        sc.pl.heatmap(
                adata,
                var_names = genes_of_interest,
                groupby = "cell_type",
                use_raw = False,
                show = False,
                cmap = "viridis",
                swap_axes = False,
                save = f"_{gene}_heatmap.png" 
                )
if __name__ == "__main__":
    main()
