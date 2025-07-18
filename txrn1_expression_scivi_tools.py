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

    # Read in expression matrix and genes, barcodes (I/O)
    adata = read_exp_matrix(matrix_file, genes, barcodes)

    # Add annotations before modeling or plotting
    annotations_dict = pyreadr.read_r(annotation_file)
    annotations_df = annotations_dict[None].set_index("cell")
    annotations_df = annotations_df.rename(columns={"truth": "cell_type"})
    adata.obs = adata.obs.join(annotations_df, how="left")


    adata, model = probabilistic_modeling(adata)
    plots(adata)
    post_processing_and_clustering(adata)
    run_de_analysis(adata, model)


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

    return adata



def probabilistic_modeling(adata):
    scvi.settings.dl_num_workers = 1
    scvi.model.SCVI.setup_anndata(adata, batch_key=None)
    model = scvi.model.SCVI(adata)
    model.train(accelerator = "mps")

    adata.obsm["X_scVI"] = model.get_latent_representation()
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.umap(adata)
    return adata, model



def plots(adata):
    genes_of_interest = ["TXNRD1", "TXN"]
    adata_log = adata.copy()

    # Log transform to make plots more interpretable 
    sc.pp.log1p(adata_log)

    for gene in tqdm(genes_of_interest):
        sc.pl.violin(adata_log, keys=gene, groupby="cell_type", show=False, stripplot=True)
        plt.gcf().set_size_inches(12, 6)
        plt.xticks(rotation=45, ha='right')
        plt.savefig(f"results/figures/violin_{gene}_cancer_vs_immune_cells.png", bbox_inches='tight', dpi=300)
        plt.close("all")

        sc.pl.dotplot(adata_log, gene, groupby="cell_type", show=False)
        plt.savefig(f"results/figures/dotplot_{gene}_cancer_vs_immune_cells.png", bbox_inches='tight', dpi=300)
        plt.close("all")

def run_de_analysis(adata, model):
    results = []
    reference_group = "Cancer cells"
    genes_of_interest = ["TXNRD1", "TXN"]

    for gene in tqdm(genes_of_interest):
        for cell_type in adata.obs["cell_type"].unique():
            if cell_type == reference_group:
                mean_expression = adata[adata.obs["cell_type"] == cell_type, gene].X.mean()
                bayes_factor = None
            else:
                de = model.differential_expression(
                    groupby="cell_type",
                    group1=reference_group,
                    group2=cell_type,
                    gene_list=[gene],
                )
                bayes_factor = de["bayes_factor"].values[0]
                mean_expression = adata[adata.obs["cell_type"] == cell_type, gene].X.mean()

            results.append({
                "Gene": gene,
                "Cell Type": cell_type,
                "Mean Expression": mean_expression,
                "Bayes Factor": bayes_factor
            })

    results_df = pd.DataFrame(results)
    results_df.sort_values(by=["Cell Type", "Gene"], inplace=True)
    results_df.to_csv("results/scvi_de_results.csv", index=False)


def post_processing_and_clustering(adata):

    # Plot UMAP
    genes_of_interest = ["TXNRD1", "TXN"]
    for gene in tqdm(genes_of_interest):
        sc.pl.umap(adata, color="cell_type", show=False, save = f"_{gene}_celltypes.pdf")
        sc.pl.umap(adata, color = gene, show=False, cmap = "inferno", vmin = 0, size = 8, vmax = "p99.5", save = f"_{gene}_expression.pdf")

if __name__ == "__main__":
    main()
