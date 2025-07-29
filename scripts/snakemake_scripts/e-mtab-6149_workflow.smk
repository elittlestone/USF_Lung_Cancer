
rule run_scvi_processing:
  input:
    expression_matrix = "data/raw/emtab_6149/matrix.mtx",
    genes = "data/raw/emtab_6149/genes.csv",
    barcodes = "data/raw/emtab_6149/barcodes.csv",
    annotation_file = "data/raw/emtab_6149/E-MTAB-6149_cell_types.rds"
  output:
    processed_h5ad = "data/processed/emtab_6149/processed_anndata.h5ad",
    vae_file = directory("data/processed/emtab_6149/vae_file"),
    violin_qc = "results/emtab_figures/violin_qc.png",
    scatter_mito = "results/emtab_figures/scatter_total_counts_vs_per_mito_counts.png",
    scatter_genes_by_counts = "results/emtab_figures/scatter_total_counts_vs_gene_counts.png",
    samples_umap = "results/emtab_figures/umap_samples_scvi.png",
    leiden_umap = "results/emtab_figures/umap_leiden_scvi_clusters.png",
    qc_umap = "results/emtab_figures/umap_qc_metrics_scvi.png"
  shell:
    """
    python scripts/python_scripts/emtab_6149_workflow.py \
    --expression_matrix_file {input.expression_matrix} --genes {input.genes} \
    --barcodes {input.barcodes} --annotation_file {input.annotation_file} \
    --processed_h5ad_file {output.processed_h5ad} \
    --vae_file {output.vae_file}
    """

