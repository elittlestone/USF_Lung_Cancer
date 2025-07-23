

rule sc_rna_seq_analysis:
  input:
    matrix_file = "data/matrix.mtx",
    genes_file = "data/genes.csv",
    barcodes_file = "data/barcodes.csv",
    annotation_file = "data/E-MTAB-6149_cell_types.rds",
    gene_list = "data/dna_repair_genes.xlsx"
  output:
    csv = "results/mannwhitney_celltype_comparison.csv",
    figures_dir = directory("results/figures")
  shell:
    "python scripts/python_scripts/expression_with_scanpy.py \
    --expression_matrix_file {input.matrix_file} \
    --genes {input.genes_file} \
    --barcodes {input.barcodes_file} \
    --annotation_file {input.annotation_file} \
    --gene_list {input.gene_list}"


