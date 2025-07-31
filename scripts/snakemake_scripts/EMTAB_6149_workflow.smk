
rule run_scvi_processing:
  input:
    expression_matrix = "data/raw/EMTAB_6149/matrix.mtx",
    genes = "data/raw/EMTAB_6149/genes.csv",
    barcodes = "data/raw/EMTAB_6149/barcodes.csv",
    annotation_file = "data/raw/EMTAB_6149/E-MTAB-6149_cell_types.rds"
  output:
    processed_h5ad = protected("data/processed/EMTAB_6149/processed_anndata.h5ad"),
    vae_file = directory("data/processed/EMTAB_6149/vae_file"),
    violin_qc = "results/EMTAB_6149/figures/violin_qc.png",
    scatter_mito = "results/EMTAB_6149/figures/scatter_total_counts_vs_per_mito_counts.png",
    scatter_genes_by_counts = "results/EMTAB_6149/figures/scatter_total_counts_vs_gene_counts.png",
    samples_umap = "results/EMTAB_6149/figures/umap_samples_scvi.png",
    leiden_umap = "results/EMTAB_6149/figures/umap_leiden_scvi_clusters.png",
    qc_umap = "results/EMTAB_6149/figures/umap_qc_metrics_scvi.png"
  shell:
    """
    python scripts/python_scripts/EMTAB_6149_workflow.py \
    --expression_matrix_file {input.expression_matrix} --genes {input.genes} \
    --barcodes {input.barcodes} --annotation_file {input.annotation_file} \
    --processed_h5ad_file {output.processed_h5ad} \
    --vae_file {output.vae_file}
    """

rule gene_specific_plots:
  input:
    processed_h5ad = "data/processed/EMTAB_6149/processed_anndata.h5ad",
    gene_list = "data/raw/EMTAB_6149/dna_repair_genes.xlsx"
  params:
    figures_dir = directory("results/EMTAB_6149/figures"),
    csv_dir = directory("results/EMTAB_6149/csvs"),
    vae = directory("data/processed/EMTAB_6149/vae_file")
  shell:
    """
    python scripts/python_scripts/EMTAB_6149_workflow.py \
    --processed_h5ad_file {input.processed_h5ad} --gene_list {input.gene_list} \
    --figures_output {params.figures_dir} --vae_file {params.vae} \
    --csv_output {params.csv_dir}
    """
