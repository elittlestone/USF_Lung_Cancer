WILDTYPE_GSM_IDS = ["GSM8462132","GSM8462130","GSM8462127","GSM8462126","GSM8462125"]
WILDTYPE_SAMPLES = ["AL11", "AL08", "AL04", "AL03", "AL02"]

# Generate URLs directly
FILTERED_H5_URLS = [f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:7]}nnn/{gsm}/suppl/{gsm}_{sample}_GEX_filtered_feature_bc_matrix.h5" 
                    for gsm, sample in zip(WILDTYPE_GSM_IDS, WILDTYPE_SAMPLES)]

rule download_h5_files:
  params:
        filtered_h5_file = FILTERED_H5_URLS
  shell:
    """
    for url in {params.filtered_h5_file}; do
    wget "$url" -P "data/wildtype_h5_files/"
    done
    """

rule run_scanpy_script:
  input:
    h5_dir = "data/wildtype_h5_files"
  output:
    conc_anndata_h5 = "data/processed/concatenated_anndata.h5ad",
    vae_file = "data/processed/vae_file",
    violin_qc = "results/geo_figures/violin_qc.png",
    scatter_mito = "results/geo_figures/scatter_total_counts_vs_per_mito_counts.png",
    scatter_genes_by_counts = "results/geo_figures/scatter_total_counts_vs_gene_counts.png",
    samples_umap = "results/geo_figures/umap_samples_scvi.png",
    leiden_umap = "results/geo_figures/umap_leiden_scvi_clusters.png",
    qc_umap = "results/geo_figures/umap_qc_metrics_scvi.png"
  shell:
    """
    python scripts/python_scripts/ALK_positive_single_cell_workflow.py \
    --h5_directory {input.h5_dir} --pre_processed_h5_file {output.conc_anndata_h5} \
    --vae_file {output.vae_file}
    """

rule run_celltype_annotation:
  input:
    processed_h5ad_file = "data/processed/concatenated_anndata.h5ad",
    vae_file = directory("data/processed/vae_file")
  output:
    celltypist_map = "results/geo_figures/umap_celltypist.png",
    dgea_csv = "data/processed/all_celltypes_DGEA.csv"
  params:
    csv_output_dir = "data/processed"
  shell:
    """
    python scripts/python_scripts/ALK_positive_single_cell_workflow.py \
    --pre_processed_h5_file {input.processed_h5ad_file} --skip_scvi_training \
    --vae_file {input.vae_file} --csv_output_dir {params.csv_output_dir}
    """
