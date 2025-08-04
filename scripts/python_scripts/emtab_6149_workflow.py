import argparse
import matplotlib.pyplot as plt
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

    # Command line arguments 
    parser = argparse.ArgumentParser(description = "Gene expression comparison script")

    parser.add_argument("--expression_matrix_file", type = str, required = False, 
    help = "Sparse expression matrix file")

    parser.add_argument("--genes", type = str, required = False, 
    help = "Genes extracted from matrix file")

    parser.add_argument("--barcodes", type = str, required = False, 
    help = "Barcodes extracted from matrix file")

    parser.add_argument("--annotation_file", type = str, required = False, 
    help = "Annotated celltypes")

    parser.add_argument("--processed_h5ad_file", type = str, required = True,
                        help = "Processed h5ad file")

    parser.add_argument("--vae_file", type = str, required = False,
                        help = "Trained model")

    parser.add_argument("--gene_list", type = str, required = False, 
    help = "List of genes for analysis")
   
    parser.add_argument("--figures_output", type = str, required = False,
                        help = "Output figures directory")
    
    parser.add_argument("--csv_output_dir", type = str, required = False,
                        help = "Output csv for statistical tests")

    args = parser.parse_args()


    # Directory where figures generated will be outputted 
    figures_output_directory = args.figures_output
    os.makedirs(figures_output_directory, exist_ok=True)
    sc.settings.figdir = figures_output_directory
    sc.settings.verbosity = 'debug'
    # Define CLI arguments 
    
    matrix_file, genes, barcodes, annotation_file, processed_h5ad_file, vae_file, gene_list, csv_output_dir = (
            args.expression_matrix_file, args.genes, args.barcodes, args.annotation_file, args.processed_h5ad_file, 
            args.vae_file, args.gene_list, args.csv_output_dir)
    
    # Read in genes to analyze from excel file 
    dna_repair_genes = pd.read_excel(gene_list)
    genes_of_interest = [dna_repair_gene for dna_repair_gene in dna_repair_genes["geneid"]]

    # Read in expression matrix and genes, barcodes
    #adata = read_exp_matrix(matrix_file, genes, barcodes)

    # Merge pre-identified celltype annotations with expression matrix 
    #adata = merge_annotations_with_expression_matrix(adata, annotation_file)

    # Perform QC steps
    #adata = quality_control(adata)

    # Preprocess/setup AnnData for scvi traiing
    #adata = pre_processing(adata)

    # Train scvi model 
    #adata = scvi_analysis_and_clustering(adata, processed_h5ad_file, vae_file)
    
    # Caluclate Plots for DNA Repair Genes
    adata = gene_specific_plots(processed_h5ad_file, genes_of_interest, figures_output_directory)

    # Perform DGEA
    os.makedirs(csv_output_dir, exist_ok = True)
    #scvi_differential_expression(adata, vae_file, genes_of_interest, csv_output_dir)
    statistical_tests(adata, genes_of_interest, csv_output_dir)



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
    """
    Read in expression matrix and join with annotated celltypes file
    """
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
    adata.obs["batch"] = adata.obs_names.str.split("_").str[-1].astype(str)
    
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
    vae = scvi.model.SCVI(adata, n_layers = 2, n_latent = 30, gene_likelihood = "zinb")

    # Train the model
    vae.train(accelerator = "mps")
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


def gene_specific_plots(processed_h5ad_file, genes_of_interest, figures_output_dir):
    """
    Create UMAPs, Violin, and Dotplots for genes of interest
    """

    adata = sc.read_h5ad(processed_h5ad_file)
    n = len(adata.obs["cell_type"].unique())
    palette = sns.color_palette("tab20", n)
    
    # UMAP of celltypes 
    sc.pl.umap(adata, color = "cell_type", show = False, palette = palette, save = "_celltypes.png")
   
    # All genes in one Dotplot 
    sc.pl.dotplot(adata, var_names = genes_of_interest, groupby = "cell_type", standard_scale = "var", show = False, save = "all_genes.png")
            
    # Generate UMAPs and Dotplots of Genes of Interest 
    for gene in tqdm(genes_of_interest):
        sc.pl.umap(adata, color = gene, show = False, cmap = "inferno", vmin = 0, size = 8, vmax = "p99.5", save = f"_{gene}.png" )
    

        # Violin Plot 
        sc.pl.violin(adata, keys=gene, groupby="cell_type", show = False, stripplot = True)
        plt.gcf().set_size_inches(12, 6)
        plt.xticks(rotation = 45, ha = 'right')  # rotate labels for better readability
        plt.savefig(f"{figures_output_dir}/violin_{gene}_cancer_vs_immune_cells.png", bbox_inches='tight', dpi=300)
        plt.close("all")
    
    return adata 

def statistical_tests(adata, genes_of_interest, csv_output):
    """
    Run Mannwhitneyu tests to test for difference in expression of given genes in cancer vs other immune cells
    @adata : Annotation data object (gene expression matrix + other metadata)
    @genes_of_interest : genes to be tested 
    """
    # Define reference group vs other celltypes to test against 
    reference_group = "Cancer cells"
    all_cell_types = adata.obs["cell_type"].unique().tolist()
    non_cancer_cells = [celltype for celltype in all_cell_types if celltype != reference_group]
   
    adata.layers["counts"] = adata.X.copy()

    # Normalize expression values 
    sc.pp.normalize_total(adata, target_sum = 1e4)

    # Log transform 
    sc.pp.log1p(adata)

    # Run wilcoxon rank sum test between cancer cells and other celltypes 
    
    sc.tl.rank_genes_groups(adata, groupby = "cell_type", groups = non_cancer_cells, 
                             reference = reference_group, method = "wilcoxon")
    
    # Loop through each gene to get test results per celltype 
    all_results = []
    for group in non_cancer_cells:
        df = sc.get.rank_genes_groups_df(adata, group = group)
        df = df[df['names'].isin(genes_of_interest)]
        df['tested_group'] = group
        all_results.append(df)

    # Output all statistical test results to csv 
    combined_df = pd.concat(all_results, ignore_index = True)
    combined_df.to_csv(os.path.join(csv_output, "wilcoxon_degs.csv"), index = False)

def scvi_differential_expression(adata, vae_file, genes_of_interest, csv_output_dir):
    """
    Differential expression testing with Bayesian approach 
    """
    # Load trained model 
    vae = scvi.model.SCVI.load(vae_file, adata)

    # Define reference group vs other celltypes 
    reference_group = "Cancer cells"
    all_cell_types = adata.obs["cell_type"].unique().tolist()
    non_cancer_cells = [celltype for celltype in all_cell_types if celltype != reference_group]

    # Initialize list to store all results of DGE 
    all_results = []

    for celltype in non_cancer_cells:
        de_results = vae.differential_expression(
            groupby = "cell_type",
            group1 = reference_group,
            group2 = celltype,
            mode = "change",
            batch_correction = True
            )
        print(de_results.head())
        print(de_results.columns)

        de_results_filtered = de_results.loc[de_results.index.isin(genes_of_interest)].copy()
        de_results_filtered['comparison'] = f"{reference_group}_vs_{celltype}"
        all_results.append(de_results_filtered)

    combined_results = pd.concat(all_results)
    combined_results.to_csv(os.path.join(csv_output_dir, "scvi_bayesian_degs.csv"))


if __name__ == "__main__":
    main()
