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
    conc_anndata_h5 = "data/processed/concatenated_anndata.h5ad"
  shell:
    """
    python scripts/python_scripts/expression_with_scanpy_wildtype_samples.py \
    --h5_directory {input.h5_dir} --pre_processed_h5_file {output.conc_anndata_h5}
    """
