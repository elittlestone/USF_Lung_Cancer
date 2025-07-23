import argparse
import pandas as pd
import pyreadr
from scipy.io import mmread
import anndata as ad
import scvi
import scanpy as sc
import os 
import matplotlib as mtplt
mtplt.use('Agg')
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
import torch

def main():

    # Directory where figures generated will be outputted 
    output_directory = "results/figures"
    os.makedirs(output_directory, exist_ok=True)
    sc.settings.figdir = output_directory
    sc.settings.verbosity = 'debug'

    # Set scvi-tools settings
    scvi.settings.seed = 42
    scvi.settings.dl_num_workers = 0  # Use 0 for compatibility
    scvi.settings.progress_bar_style = "tqdm"

    parser = argparse.ArgumentParser(description = "scvi-tools enhanced expression comparison script")

    parser.add_argument("--expression_matrix_file", type = str, required = True, 
    help = "Sparse expression matrix file")

    parser.add_argument("--genes", type = str, required = True, 
    help = "Genes extracted from matrix file")

    parser.add_argument("--barcodes", type = str, required = True, 
    help = "Barcodes extracted from matrix file")

    parser.add_argument("--annotation_file", type = str, required = True, 
    help = "Annotated celltypes")

    parser.add_argument("--gene_list", type = str, required = True, 
    help = "Genes of interest for analysis")

    parser.add_argument("--n_epochs", type = int, default = 100,
    help = "Number of training epochs for scVI model")

    parser.add_argument("--n_latent", type = int, default = 30,
    help = "Dimensionality of latent space")

    parser.add_argument("--n_layers", type = int, default = 1,
    help = "Number of hidden layers in encoder/decoder")

    args = parser.parse_args()
    
    # Read in expression matrix and genes, barcodes (I/O)
    adata = read_exp_matrix(args.expression_matrix_file, args.genes, args.barcodes)

    # Add annotations before modeling or plotting
    annotations_dict = pyreadr.read_r(args.annotation_file)
    annotations_df = annotations_dict[None].set_index("cell")
    annotations_df = annotations_df.rename(columns={"truth": "cell_type"})
    adata.obs = adata.obs.join(annotations_df, how="left")

    # Quality control and preprocessing
    adata = quality_control(adata)
    
    # Probabilistic modeling with scvi-tools
    adata, model = probabilistic_modeling(adata, args.n_latent, args.n_layers, args.n_epochs)
    
    # Generate plots using scVI representations
    plots(adata, args.gene_list)
    
    # Dimensionality reduction and clustering using scVI
    post_processing_and_clustering(adata, args.gene_list)
    
    # Differential expression analysis using scVI
    run_de_analysis(adata, model)
    
    # Additional scvi-tools analyses
    run_additional_scvi_analyses(adata, model)


def read_exp_matrix(matrix_file, genes, barcodes):
    """Load expression matrix and create AnnData object"""
    # Load matrix
    matrix = mmread(matrix_file).tocsr().T

    # Load gene names and cell IDs
    genes = pd.read_table(genes, header = None)[0].astype(str).tolist()
    barcodes = pd.read_table(barcodes, header = None)[0].astype(str).tolist()

    # Create annotation data object 
    adata = sc.AnnData(matrix)
    adata.var_names = genes
    adata.obs_names = barcodes
    
    # Make variable names unique (required for scvi-tools)
    adata.var_names_make_unique()

    return adata


def quality_control(adata):
    """Enhanced quality control using scanpy and scvi-tools best practices"""
    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
    
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    
    # Add mitochondrial and ribosomal gene percentages
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
    
    print(f"Number of cells before filtering: {adata.n_obs}")
    print(f"Number of genes before filtering: {adata.n_vars}")
    
    # Basic filtering - adjust thresholds as needed
    sc.pp.filter_cells(adata, min_genes=200)  # Filter cells with too few genes
    sc.pp.filter_genes(adata, min_cells=3)    # Filter genes expressed in too few cells
    
    print(f"Number of cells after filtering: {adata.n_obs}")
    print(f"Number of genes after filtering: {adata.n_vars}")
    
    return adata


def probabilistic_modeling(adata, n_latent=30, n_layers=1, n_epochs=100):
    """Enhanced probabilistic modeling with scvi-tools"""
    
    # Setup AnnData for scVI - specify categorical covariates
    scvi.model.SCVI.setup_anndata(
        adata, 
        batch_key=None,  # Add batch key if you have batch effects
        categorical_covariate_keys=["cell_type"],  # Use cell type as covariate
        continuous_covariate_keys=None
    )
    
    # Create scVI model with enhanced parameters
    model = scvi.model.SCVI(
        adata, 
        n_latent=n_latent,
        n_layers=n_layers,
        dropout_rate=0.1,
        dispersion="gene",
        gene_likelihood="zinb"  # Zero-inflated negative binomial
    )
    
    # Train the model
    print("Training scVI model...")
    model.train(
        max_epochs=n_epochs,
        accelerator="auto",  # Automatically detect available accelerator
        enable_progress_bar=True,
        check_val_every_n_epoch=10
    )
    
    # Get model outputs
    print("Extracting latent representations and normalized expression...")
    adata.obsm["X_scVI"] = model.get_latent_representation()
    adata.layers["scvi_normalized"] = model.get_normalized_expression(library_size=10e4)
    
    # Get imputed values (denoised gene expression)
    adata.layers["scvi_imputed"] = model.get_normalized_expression(
        transform_batch=None,
        library_size="latent"
    )
    
    # Compute neighbors and UMAP using scVI representation
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15)
    sc.tl.umap(adata, min_dist=0.5)
    
    return adata, model


def plots(adata, gene_list):
    """Generate plots using both raw and scVI-processed data"""
    
    # Load DNA repair genes from Excel file
    dna_repair_genes = pd.read_excel(gene_list)
    genes_of_interest = [gene for gene in dna_repair_genes["geneid"] if gene in adata.var_names]
    
    print(f"Found {len(genes_of_interest)} genes of interest in the dataset")
    
    # Create a copy with log-transformed data for visualization
    adata_log = adata.copy()
    sc.pp.normalize_total(adata_log, target_sum=1e4)
    sc.pp.log1p(adata_log)
    
    for gene in tqdm(genes_of_interest, desc="Creating plots"):
        # Violin plot using log-transformed data
        sc.pl.violin(adata_log, keys=gene, groupby="cell_type", show=False, stripplot=True)
        plt.gcf().set_size_inches(14, 8)
        plt.xticks(rotation=45, ha='right')
        plt.title(f'{gene} Expression Distribution by Cell Type')
        plt.savefig(f"results/figures/violin_{gene}_scvi_analysis.png", bbox_inches='tight', dpi=300)
        plt.close("all")

        # Dotplot using log-transformed data
        sc.pl.dotplot(adata_log, [gene], groupby="cell_type", show=False)
        plt.title(f'{gene} Expression Summary by Cell Type')
        plt.savefig(f"results/figures/dotplot_{gene}_scvi_analysis.png", bbox_inches='tight', dpi=300)
        plt.close("all")
        
        # Violin plot using scVI normalized data
        adata_scvi_norm = adata.copy()
        adata_scvi_norm.X = adata.layers["scvi_normalized"]
        sc.pl.violin(adata_scvi_norm, keys=gene, groupby="cell_type", show=False, stripplot=True)
        plt.gcf().set_size_inches(14, 8)
        plt.xticks(rotation=45, ha='right')
        plt.title(f'{gene} Expression (scVI Normalized) by Cell Type')
        plt.savefig(f"results/figures/violin_{gene}_scvi_normalized.png", bbox_inches='tight', dpi=300)
        plt.close("all")


def run_de_analysis(adata, model):
    """Enhanced differential expression analysis using scvi-tools"""
    
    # Load genes of interest
    dna_repair_genes = pd.read_excel("gene_list.xlsx")
    genes_of_interest = [gene for gene in dna_repair_genes["geneid"] if gene in adata.var_names]
    
    results = []
    reference_group = "Cancer cells"
    
    print("Running differential expression analysis...")
    
    # Get all unique cell types
    cell_types = adata.obs["cell_type"].unique()
    
    for cell_type in tqdm(cell_types, desc="Processing cell types"):
        if cell_type == reference_group:
            continue
            
        # Run differential expression for all genes of interest at once
        try:
            de_results = model.differential_expression(
                groupby="cell_type",
                group1=reference_group,
                group2=cell_type,
                delta=0.25,  # Minimum fold change threshold
                batch_size=None,
                all_stats=True,
                batch_correction=False,
                batchid1=None,
                batchid2=None,
                fdr_target=0.05,
                silent=True
            )
            
            # Filter for genes of interest
            de_filtered = de_results[de_results.index.isin(genes_of_interest)]
            
            # Add results for each gene
            for gene in genes_of_interest:
                if gene in de_filtered.index:
                    gene_stats = de_filtered.loc[gene]
                    bayes_factor = gene_stats['bayes_factor']
                    lfc_mean = gene_stats['lfc_mean']
                    lfc_std = gene_stats['lfc_std']
                    raw_mean1 = gene_stats['raw_mean1']  # Reference group mean
                    raw_mean2 = gene_stats['raw_mean2']  # Comparison group mean
                else:
                    bayes_factor = None
                    lfc_mean = None
                    lfc_std = None
                    raw_mean1 = None
                    raw_mean2 = None
                
                results.append({
                    "Gene": gene,
                    "Reference_Group": reference_group,
                    "Comparison_Group": cell_type,
                    "Reference_Mean": raw_mean1,
                    "Comparison_Mean": raw_mean2,
                    "Log_Fold_Change_Mean": lfc_mean,
                    "Log_Fold_Change_Std": lfc_std,
                    "Bayes_Factor": bayes_factor
                })
                
        except Exception as e:
            print(f"Error processing {cell_type}: {str(e)}")
            continue
    
    # Also add reference group data
    for gene in genes_of_interest:
        if gene in adata.var_names:
            ref_mean = adata[adata.obs["cell_type"] == reference_group, gene].X.mean()
        else:
            ref_mean = None
            
        results.append({
            "Gene": gene,
            "Reference_Group": reference_group,
            "Comparison_Group": reference_group,
            "Reference_Mean": ref_mean,
            "Comparison_Mean": ref_mean,
            "Log_Fold_Change_Mean": 0.0,
            "Log_Fold_Change_Std": 0.0,
            "Bayes_Factor": None
        })

    # Save results
    results_df = pd.DataFrame(results)
    results_df.sort_values(by=["Gene", "Comparison_Group"], inplace=True)
    results_df.to_csv("results/scvi_enhanced_de_results.csv", index=False)
    
    print(f"Saved differential expression results for {len(genes_of_interest)} genes")


def run_additional_scvi_analyses(adata, model):
    """Additional analyses leveraging scvi-tools capabilities"""
    
    # 1. Generate posterior samples for uncertainty quantification
    print("Generating posterior samples for uncertainty analysis...")
    posterior_samples = model.get_latent_representation(give_mean=False, n_samples=10)
    adata.obsm["X_scVI_samples"] = posterior_samples
    
    # 2. Compute gene expression programs/factors
    print("Computing gene loadings and expression programs...")
    try:
        gene_loadings = model.get_loadings()
        np.save("results/scvi_gene_loadings.npy", gene_loadings)
        
        # Save top genes for each latent dimension
        loadings_df = pd.DataFrame(gene_loadings, index=adata.var_names)
        top_genes_per_factor = {}
        
        for factor in range(gene_loadings.shape[1]):
            top_genes = loadings_df.iloc[:, factor].abs().nlargest(20)
            top_genes_per_factor[f"Factor_{factor}"] = top_genes.index.tolist()
        
        pd.DataFrame(dict([(k, pd.Series(v)) for k, v in top_genes_per_factor.items()])).to_csv(
            "results/scvi_top_genes_per_factor.csv", index=False
        )
        
    except Exception as e:
        print(f"Could not compute gene loadings: {str(e)}")
    
    # 3. Model evaluation metrics
    print("Computing model evaluation metrics...")
    try:
        # Reconstruction error
        reconstruction_error = model.get_reconstruction_error()
        adata.obs['reconstruction_error'] = reconstruction_error['reconstruction_loss']
        
        # Marginal log likelihood
        marginal_ll = model.get_marginal_ll(n_mc_samples=1000)
        print(f"Marginal log likelihood: {marginal_ll}")
        
        # Save model evaluation metrics
        eval_metrics = {
            'marginal_log_likelihood': marginal_ll,
            'mean_reconstruction_error': reconstruction_error['reconstruction_loss'].mean(),
            'std_reconstruction_error': reconstruction_error['reconstruction_loss'].std()
        }
        
        pd.DataFrame([eval_metrics]).to_csv("results/scvi_model_evaluation.csv", index=False)
        
    except Exception as e:
        print(f"Could not compute all evaluation metrics: {str(e)}")


def post_processing_and_clustering(adata, gene_list):
    """Post-processing and visualization using scVI representations"""

    # Load genes of interest
    dna_repair_genes = pd.read_excel(gene_list)
    genes_of_interest = [gene for gene in dna_repair_genes["geneid"] if gene in adata.var_names]
    
    # Perform Leiden clustering on scVI representation
    sc.tl.leiden(adata, resolution=0.5, key_added="leiden_scvi")
    
    for gene in tqdm(genes_of_interest, desc="Creating UMAP plots"):
        # UMAP colored by cell types
        sc.pl.umap(adata, color="cell_type", show=False, 
                  title=f"Cell Types - {gene} Analysis", 
                  save=f"_{gene}_celltypes_scvi.pdf")
        
        # UMAP colored by Leiden clusters
        sc.pl.umap(adata, color="leiden_scvi", show=False, 
                  title=f"Leiden Clusters - {gene} Analysis",
                  save=f"_{gene}_leiden_clusters_scvi.pdf")
        
        # UMAP colored by gene expression (raw)
        sc.pl.umap(adata, color=gene, show=False, cmap="inferno", 
                  vmin=0, size=8, vmax="p99.5",
                  title=f"{gene} Raw Expression",
                  save=f"_{gene}_expression_raw_scvi.pdf")
        
        # UMAP colored by scVI-normalized expression
        adata_temp = adata.copy()
        adata_temp.X = adata.layers["scvi_normalized"]
        sc.pl.umap(adata_temp, color=gene, show=False, cmap="inferno", 
                  vmin=0, size=8, vmax="p99.5",
                  title=f"{gene} scVI Normalized Expression",
                  save=f"_{gene}_expression_scvi_normalized.pdf")


if __name__ == "__main__":
    main()
