import argparse
import matplotlib.pyplot as plt
import numpy as np
import os 
import pandas as pd
import pyreadr
import scanpy as sc
from scipy.io import mmread
import scvi
import seaborn as sns
import torch
from tqdm import tqdm

def main():

    # Directory where figures generated will be outputted 
    output_directory = "results/E-MTAB-6149_figures"
    os.makedirs(output_directory, exist_ok=True)
    sc.settings.figdir = output_directory
    sc.settings.verbosity = 'debug'

    # Command line arguments 
    parser = argparse.ArgumentParser(description = "Gene expression comparison script")

    parser.add_argument("--expression_matrix_file", type = str, required = True, 
    help = "Sparse expression matrix file")

    parser.add_argument("--genes", type = str, required = True, 
    help = "Genes extracted from matrix file")

    parser.add_argument("--barcodes", type = str, required = True, 
    help = "Barcodes extracted from matrix file")

    parser.add_argument("--annotation_file", type = str, required = True, 
    help = "Annotated celltypes")

    parser.add_argument("--processed_h5ad_file", type = str, required = True,
                        help = "Processed h5ad file")

    parser.add_argument("--vae_file", type = str, required = True,
                        help = "Trained model")

    parser.add_argument("--gene_list", type = str, required = False, 
    help = "List of genes for analysis")
   
    parser.add_argument("--csv_output", type = str, required = False,
                        help = "Output csv for statistical tests")
    args = parser.parse_args()

    # Read in genes to analyze from excel file 
    #dna_repair_genes = pd.read_excel(args.gene_list)
    #genes_of_interest = [gene for gene in dna_repair_genes["geneid"]]

    
    matrix_file, genes, barcodes, annotation_file, processed_h5ad_file, vae_file, gene_list, csv_output = (
            args.expression_matrix_file, args.genes, args.barcodes, args.annotation_file, 
                           args.processed_h5ad_file, args.vae_file, args.gene_list, args.csv_output)

    # Read in expression matrix and genes, barcodes
    adata = read_exp_matrix(matrix_file, genes, barcodes)

    # Merge pre-identified celltype annotations with expression matrix 
    adata = merge_annotations_with_expression_matrix(adata, annotation_file)

    # Perform QC steps
    adata = quality_control(adata)

    # Preprocess/setup AnnData for scvi traiing
    adata = pre_processing(adata)

    # Train scvi model 
    scvi_analysis_and_clustering(adata, processed_h5ad_file, vae_file)
    
    # Perform DGEA
    #statistical_tests(adata, genes_of_interest, csv_output)



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


    return adata


def merge_annotations_with_expression_matrix(adata, annotation_file):

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

    # Add batch labels to anndata object 
    adata.obs["batch"] = adata.obs_names.str.split("_").str[-1].astype(int)
    
    return adata 


def quality_control(adata):
    """
    Calculate QC metrics on data
    @adata: AnnData object containing concatenated gene expression matrices for all samples 
    """
    print("Calculating QC metrics...")

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.startswith("^HB[^(P)]")

    # Calculate QC metrics 
    sc.pp.calculate_qc_metrics(adata,
                               qc_vars = ["mt", "ribo", "hb"],
                               percent_top = None,
                               log1p = False,
                               inplace = True)


    #  QC plots
    
    # Violin plot
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, save = "_qc.png")

    # Scatter plot: total counts vs. mito %
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save = "_total_counts_vs_per_mito_counts.png")

    # Scatter: total counts vs gene counts
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save = "_total_counts_vs_gene_counts.png")


    print("Before filtering:")
    print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")
    print(adata.obs["pct_counts_mt"].describe())
    # Filter cells with too few genes or too many genes 
    sc.pp.filter_cells(adata, min_genes = 200)
    sc.pp.filter_genes(adata, min_cells = 3)

    print("After removing cells without a minimum of 200 genes expressed, and genes not present in a minimum of 3 cells:")
    print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")


    # Filter cells with too high mito gene percentage
    adata = adata[adata.obs['pct_counts_mt'] < 20, :]

    print(f"After filtering cells with mito counts percentage > 20\nCells: {adata.n_obs}, Genes: {adata.n_vars}")
  
    return adata


def pre_processing(adata):
    """
    Setup AnnData object for scvi-tools analysis
    """

    adata = adata.copy()

    # Setup AnnData object for scVI
    print("Setting up AnnData object for scvi-tools...")
    scvi.model.SCVI.setup_anndata(
            adata,
            layer = None,
            batch_key = "batch",
            labels_key = None,
            continuous_covariate_keys = ["total_counts", "pct_counts_mt"]
            )

    return adata

def scvi_analysis_and_clustering(adata, processed_h5_file, vae_file):
    """
    Perform scvi-tools analysis: train mode, get latent representation, perform clustering and UMAP
    """
    print("Starting scvi-tools analysis...")

    # Check if CUDA is available

    print(f"CUDA available: {torch.cuda.is_available()}") 
    if torch.cuda.is_available():
        print(f"Using GPU: {torch.cuda.get_device_name()}")

    # Train scvi model
    print("Training scVI model...")
    vae = scvi.model.SCVI(adata, n_layers = 2, n_latent = 30, gene_likelihood = "nb")

    # Train the model
    vae.train(max_epochs = 400, plan_kwargs = {"lr": 1e-3}, check_val_every_n_epoch = 10, accelerator = "gpu")
    vae.save(vae_file)

    print("Getting latent representation...")
    adata.obsm["X_scVI"] = vae.get_latent_representation()

    # Get normalized expression
    adata.layers["scvi_normalized"] = vae.get_normalized_expression(library_size = 1e4)

    # Perform clustering on scVi latent space
    print("Performing Clustering...")
    sc.pp.neighbors(adata, use_rep = "X_scVI", n_neighbors = 15, n_pcs = 30)
    sc.tl.leiden(adata, resolution = 0.5, key_added = "leiden_scvi")

    # Compute UMAP on scVI latent space 
    print("Computing UMAP...")
    sc.tl.umap(adata, min_dist = 0.3)

    # Plot results
    print("Creating visualizations...")

    # UMAP colored by samples 
    sc.pl.umap(adata, color = "batch", show = False, save = "_samples_scvi.png")

    # UMAP colored by clusters 
    sc.pl.umap(adata, color = "leiden_scvi", show = False, save = "_leiden_scvi_clusters.png")

    # UMAP colored by QC metrics 
    sc.pl.umap(adata, color = ["total_counts", "n_genes_by_counts", "pct_counts_mt"],
               show = False, save = "_qc_metrics_scvi.png")
   
    os.makedirs(os.path.dirname(processed_h5_file), exist_ok = True)
    adata.write(processed_h5_file)
    

    return adata, vae



def summary_plots(adata, genes_of_interest):
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
    pass


    return adata




if __name__ == "__main__":
    main()
