from pptx import Presentation
from pptx.util import Inches, Pt
from collections import defaultdict
import os
import re

# Folder for plots to be stored in
plots_dir = "results/figures"
plot_files = [file for file in os.listdir(plots_dir) if not file.startswith("umap_celltypes")]
cell_type_umap_path = os.path.join(plots_dir, "umap_celltypes.png")

# Group plots by gene
gene_to_plots = defaultdict(dict)

for file in plot_files:
    match = re.match(r"(umap|violin|dotplot)_(.+?)_", file)
    if match:
        plot_type, gene = match.groups()
        gene_to_plots[gene][plot_type] = os.path.join(plots_dir, file)

pres = Presentation()
blank_slide_layout = pres.slide_layouts[6]

for gene, plots in gene_to_plots.items():
    slide = pres.slides.add_slide(blank_slide_layout)
    
    # Title of slide
    left = top = Inches(0.5)
    width = Inches(9)
    height = Inches(0.5)
    title_shape = slide.shapes.add_textbox(left, top, width, height)
    tf = title_shape.text_frame
    tf.text = f"{gene} Expression Overview"
    tf.paragraphs[0].font.size = Pt(28)
    
    # Add cell type UMAP 
    slide.shapes.add_picture(cell_type_umap_path, Inches(0.5), Inches(1.0), width = Inches(4.5))
    if "umap" in plots:
        slide.shapes.add_picture(plots["umap"], Inches(5.2), Inches(1.0), width = Inches(4.5))
    if "violin" in plots:
        slide.shapes.add_picture(plots["violin"], Inches(0.5), Inches(4.5), width = Inches(4))
    if "dotplot" in plots:
        slide.shapes.add_picture(plots["dotplot"], Inches(5.2), Inches(4.5), width = Inches(4.5))

pres.save("results/scanpy_gene_expression_summary.pptx")

