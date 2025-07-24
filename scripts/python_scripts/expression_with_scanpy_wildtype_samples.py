import argparse
import pandas as pd
import scanpy as sc
import os 
import glob
import matplotlib as mtplt
mtplt.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

def main():

    # Directory where figures generated will be outputted 
    output_directory = "results/geo_figures"
    os.makedirs(output_directory, exist_ok=True)
    sc.settings.figdir = output_directory
    sc.settings.verbosity = 'debug'

    # Command line arguments 
    parser = argparse.ArgumentParser(description = "Gene expression comparison script")

    parser.add_argument("--h5_directory", type = str, required = True, 
    help = "Sparse expression matrix file")

    parser.add_argument("--pre_processed_h5_file", type = str, required = True, 
                        help = "Pre processed and concatenated h5 file")

    parser.add_argument("--gene_list", type = str, required = False, 
    help = "List of genes for analysis")
   
    parser.add_argument("--csv_output", type = str, required = False,
                        help = "Output csv for statistical tests")
    args = parser.parse_args()

    # Read in genes to analyze from excel file 
    #dna_repair_genes = pd.read_excel(args.gene_list)
    #genes_of_interest = [gene for gene in dna_repair_genes["geneid"]]

    
    # Read in expression matrix files
    h5_directory = args.h5_directory
    pre_processed_and_concatenated_h5 = args.pre_processed_h5_file
    adata = read_and_concat_h5_files(h5_directory)
    pre_processing(adata, pre_processed_and_concatenated_h5)

    # Create summary plots
    #adata = summary_plots(adata, genes_of_interest)

    # Perform DGEA
    #statistical_tests(adata, genes_of_interest, csv_output)

    # Clustering and UMAP 
    #post_processing_and_clustering(adata, genes_of_interest)


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

def pre_processing(adata_combined, pre_processed_h5_file):

    # Normalize expression
    #sc.pp.normalize_total(adata_combined)

    # Log transformation
    #sc.pp.log1p(adata_combined)

    # Store raw counts before scaling 
    adata_combined.raw = adata_combined.copy()

    # Scale 
    #sc.pp.scale(adata_combined)
    
    os.makedirs(os.path.dirname(pre_processed_h5_file), exist_ok = True)
    adata_combined.write(pre_processed_h5_file)

    return adata_combined 

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



def statistical_tests(adata, genes_of_interest, csv_output):
    pass
if __name__ == "__main__":
    main()
