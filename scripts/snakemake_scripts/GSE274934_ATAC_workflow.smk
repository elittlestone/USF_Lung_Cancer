gsm_ids = config["single_cell_ATAC"]["gsm_ids"]
sample_ids = config["single_cell_ATAC"]["sample_ids"]


ATAC_FILES = ["filtered_peak_bc_matrix.h5"]


def atac_url(gsm, sample, filename):
  return f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:7]}nnn/{gsm}/suppl/{gsm}_{sample}_ATAC_{filename}"


def output_path(gsm, sample, filename):
    return f"data/raw/GSE274934/ATAC/{gsm}_{sample}_ATAC_{filename}"



all_files = [
    f"data/raw/GSE274934/ATAC/{gsm}_{sample}_ATAC_{filename}"
    for gsm, sample in zip(gsm_ids, sample_ids)
    for filename in ATAC_FILES
]

rule all:
  input:
    all_files

rule download_atac_files:
  output:
    path = "data/raw/GSE274934/ATAC/{gsm}_{sample}_ATAC_{filename}"
  params:
    url = lambda wildcards: atac_url(wildcards.gsm, wildcards.sample, wildcards.filename)
  shell:
    """
    wget -nc {params.url} -O {output.path}
    """
