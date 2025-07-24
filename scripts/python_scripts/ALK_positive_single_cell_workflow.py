import argparse
import celltypist
import glob
import matplotlib as mtplt
mtplt.use('Agg')
import os
import pandas as pd
import scanpy as sc
import scvi
import torch
from tqdm import tqdm

def main():

    # Directory where figures generated will be outputted 
    output_directory = "results/geo_figures"
    os.makedirs(output_directory, exist_ok=True)
    sc.settings.figdir = output_directory
    sc.settings.verbosity = 'debug'

    # Command line arguments 
    parser = argparse.ArgumentParser(description = "Gene expression comparison script")

    parser.add_argument("--h5_directory", type = str, required = False, 
    help = "Expression matrix file (required if not skipping scvi training steps)")

    parser.add_argument("--pre_processed_h5_file", type = str, required = True, 
                        help = "Pre processed and concatenated h5 file")

    parser.add_argument("--vae_file", type = str, required = True, 
                        help = "Trained model file")

    parser.add_argument("--skip_scvi_training", action = "store_true",
                        help = "Skip scVI training and load preprocessed AnnData from file")

    parser.add_argument("--csv_output_dir", type = str, required = False,
                        help = "Output csv for statistical tests")

    args = parser.parse_args()

    
    # Read in expression matrix files
    h5_directory = args.h5_directory
    processed_and_concatenated_h5 = args.pre_processed_h5_file
    vae_file = args.vae_file
    csv_output_dir = args.csv_output_dir


    if args.skip_scvi_training:
        print("Skipping scVI training; loading preprocessed AnnData")
        adata_combined = sc.read_h5ad(processed_and_concatenated_h5)
        celltype_annotation_with_celltypist(adata_combined)
        differential_gene_expression(adata_combined, vae_file, csv_output_dir)
    else:
        adata_combined = read_and_concat_h5_files(h5_directory)
        adata_combined = quality_control(adata_combined)
        adata_combined = pre_processing(adata_combined)
        adata_combined, vae = scvi_analysis_and_clustering(adata_combined, processed_and_concatenated_h5, vae_file)
        #celltype_annotation_with_celltypist(adata_combined)
        #differential_gene_expression(adata_combined, vae, csv_output_dir)

def read_and_concat_h5_files(h5_directory):
    """
    Function for reading multiple h5 files and concatenating them into one AnnData object
    @h5_directory : directory containing h5 gene expression matrix files 
    """
    # Get all files in directory 
    h5_files = glob.glob(os.path.join(h5_directory, "*.h5"))

    adata_list = []

    for h5_file in tqdm(h5_files, desc = "Reading in H5 files..."):
        # Load matrix into anndata object
        adata = sc.read_10x_h5(filename = h5_file)
        adata.var_names_make_unique()

        # Get sample info from filename
        filename = os.path.basename(h5_file)
        filename_parts = filename.replace(".h5", "").split("_")
        gsm_id = filename_parts[0]
        sample_name = filename_parts[1]

        # Add sample metadata to observations
        adata.obs["sample_id"] = sample_name
        adata.obs["gsm_id"] = gsm_id
        adata.obs["batch"] = f"{gsm_id}_{sample_name}"

        # Add sample prefix to barcodes so they are unique when concatenated into single anndata object
        adata.obs_names = [f"{sample_name}_{barcode}" for barcode in adata.obs_names]
        print(f"Loaded {adata.n_obs} cells from {sample_name}")
        adata_list.append(adata)

    # Concatenate all AnnData objects 
    print("Concatenating all samples...")
    adata_combined = sc.concat(adata_list,
                               axis = 0,
                               join = "outer",
                               index_unique = None,
                               fill_value = 0,
                               pairwise = False)
    print(f"Combined dataset: {adata_combined.n_obs} cells x {adata_combined.n_vars} genes")
    print(f"Samples in dataset: {adata_combined.obs['sample_id'].unique()}")

    return adata_combined

def quality_control(adata_combined):
    """
    Calculate QC metrics on data
    @adata_combined: AnnData object containing concatenated gene expression matrices for all samples 
    """
    print("Calculating QC metrics...")

    adata_combined.var["mt"] = adata_combined.var_names.str.startswith("MT-")
    adata_combined.var["ribo"] = adata_combined.var_names.str.startswith(("RPS", "RPL"))
    adata_combined.var["hb"] = adata_combined.var_names.str.startswith("^HB[^(P)]")

    # Calculate QC metrics 
    sc.pp.calculate_qc_metrics(adata_combined,
                               qc_vars = ["mt", "ribo", "hb"],
                               percent_top = None,
                               log1p = False,
                               inplace = True)


    #  QC plots
    
    # Violin plot
    sc.pl.violin(adata_combined, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, save = "_qc.png")

    # Scatter plot: total counts vs. mito %
    sc.pl.scatter(adata_combined, x='total_counts', y='pct_counts_mt', save = "_total_counts_vs_per_mito_counts.png")

    # Scatter: total counts vs gene counts
    sc.pl.scatter(adata_combined, x='total_counts', y='n_genes_by_counts', save = "_total_counts_vs_gene_counts.png")


    print("Before filtering:")
    print(f"Cells: {adata_combined.n_obs}, Genes: {adata_combined.n_vars}")
    print(adata_combined.obs["pct_counts_mt"].describe())
    # Filter cells with too few genes or too many genes 
    sc.pp.filter_cells(adata_combined, min_genes = 200)
    sc.pp.filter_genes(adata_combined, min_cells = 3)

    print("After removing cells without a minimum of 200 genes expressed, and genes not present in a minimum of 3 cells:")
    print(f"Cells: {adata_combined.n_obs}, Genes: {adata_combined.n_vars}")


    # Filter cells with too high mito gene percentage
    adata_combined = adata_combined[adata_combined.obs['pct_counts_mt'] < 20, :]

    print(f"After filtering cells with mito counts percentage > 20\nCells: {adata_combined.n_obs}, Genes: {adata_combined.n_vars}")
  
    return adata_combined 


def pre_processing(adata_combined):
    """
    Setup AnnData object for scvi-tools analysis 
    @adata_combined: AnnData object with concatenated gene expression matrices from each sample
    @pre_processed_h5_file: AnnData written to h5ad file 
    """

    # Setup AnnData object for SCVI
    print("Setting up AnnData object for scvi-tools...")

    adata_combined.raw = adata_combined.copy()

    # Only use highly variable genes 
    sc.pp.highly_variable_genes(
            adata_combined,
            flavor = "seurat_v3",
            n_top_genes = 2000,
            batch_key = "sample_id",
            subset = True
            )

    scvi.model.SCVI.setup_anndata(
            adata_combined,
            layer = None,
            batch_key = "sample_id",
            labels_key = None,
            continuous_covariate_keys = ["total_counts", "pct_counts_mt"]
            )


    return adata_combined 



def scvi_analysis_and_clustering(adata_combined, processed_h5_file, vae_file):
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
    vae = scvi.model.SCVI(adata_combined, n_layers = 2, n_latent = 30, gene_likelihood = "nb")

    # Train the model
    vae.train(max_epochs = 400, plan_kwargs = {"lr": 1e-3}, check_val_every_n_epoch = 10, accelerator = "mps")
    vae.save(vae_file)

    print("Getting latent representation...")
    adata_combined.obsm["X_scVI"] = vae.get_latent_representation()

    # Get normalized expression
    adata_combined.layers["scvi_normalized"] = vae.get_normalized_expression(library_size = 1e4)

    # Perform clustering on scVi latent space
    print("Performing Clustering...")
    sc.pp.neighbors(adata_combined, use_rep = "X_scVI", n_neighbors = 15, n_pcs = 30)
    sc.tl.leiden(adata_combined, resolution = 0.5, key_added = "leiden_scvi")

    # Compute UMAP on scVI latent space 
    print("Computing UMAP...")
    sc.tl.umap(adata_combined, min_dist = 0.3)

    # Plot results
    print("Creating visualizations...")

    # UMAP colored by samples 
    sc.pl.umap(adata_combined, color = "sample_id", show = False, save = "_samples_scvi.png")

    # UMAP colored by clusters 
    sc.pl.umap(adata_combined, color = "leiden_scvi", show = False, save = "_leiden_scvi_clusters.png")

    # UMAP colored by QC metrics 
    sc.pl.umap(adata_combined, color = ["total_counts", "n_genes_by_counts", "pct_counts_mt"],
               show = False, save = "_qc_metrics_scvi.png")
   
    os.makedirs(os.path.dirname(processed_h5_file), exist_ok = True)
    adata_combined.write(processed_h5_file)
    
    return adata_combined, vae

def celltype_annotation_with_celltypist(adata_combined, model_name = "Immune_All_Low.pkl", 
                                        leiden_key = "leiden_scvi", results_key = "celltypist"):
    """
    Annotate cell types using CellTypist and add results to AnnData object 
    """

    # Load CellTypist pre-trained model of choice 
    print("Loading CellTypist model...")

    model = celltypist.models.download_models(model = model_name)
    
    # Perform cell-type annotation on our data using pretrained model

    # Log transform normalized expression matrix
    adata_for_celltypist = adata_combined.copy()
    adata_for_celltypist.X = adata_combined.layers["scvi_normalized"].copy()
    sc.pp.log1p(adata_for_celltypist)



    print("Running CellTypist prediction...")
    cluster_labels = adata_for_celltypist.obs[leiden_key].to_list()
    predictions = celltypist.annotate(
            adata_for_celltypist,
            model = model,
            majority_voting = True, 
            over_clustering = cluster_labels)
    
    print(type(predictions.predicted_labels))
    print(predictions.predicted_labels.shape)
    print(predictions.predicted_labels.head())

    # Add predicted labels to AnnData 
    adata_for_celltypist.obs[f"{results_key}_predicted_labels"] = predictions.predicted_labels["majority_voting"].values
    
    # Plot UMAP colored by predicted labels 
    print("Plotting UMAP with CellTypist annotations...")
    sc.pl.umap(
            adata_for_celltypist, 
            color = f"{results_key}_predicted_labels",
            legend_loc = "on data",
            title = "CellTypist Cell Type Annotations",
            save = f"_{results_key}.png"
            )

    return adata_for_celltypist



def differential_gene_expression(adata_combined, vae, csv_output_dir):
    celltype_key = "celltypist_predicted_labels"
    output_dir = csv_output_dir
    vae = scvi.model.SCVI.load(vae, adata = adata_combined)
    os.makedirs(output_dir, exist_ok = True)
    celltypes = adata_combined.obs[celltype_key].unique()
    
    all_diff_expr_results = []

    # Iterate through all celltypes 
    for celltype in celltypes:
        print(f"Running scVI DGEA for cell type: {celltype}")

        # Boolean masks for group and reference
        group_mask = adata_combined.obs[celltype_key] == celltype
        ref_mask = adata_combined.obs[celltype_key] != celltype

        # Run differential expression : group vs rest of cells
        differential_expr_results = vae.differential_expression(
                group1 = group_mask,
                group2 = ref_mask,
                mode = "change",
                delta = 0.5,
                output_file = os.path.join(csv_output_dir, f"DGEA_{celltype}.csv")
                )
    
        differential_expr_results["celltype"] = celltype
        all_diff_expr_results.append(differential_expr_results)

    all_differential_expr_results_dataframe = pd.concat(all_diff_expr_results)
    
    # Sort by celltype and then gene name 
    all_differential_expr_results_dataframe = all_differential_expr_results_dataframe.sort_values(
            by = ["celltype", all_differential_expr_results_dataframe.index.name or "gene"], ascending = [True, True])

    # Save to CSV file 
    all_differential_expr_results_dataframe.to_csv(os.path.join(output_dir, "all_celltypes_DGEA.csv"))
if __name__ == "__main__":
    main()
