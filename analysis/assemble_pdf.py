"""Assemble All_figures.pdf (17 pages) from the figure folder, in presentation order.

Pages keep their own size (7.2 in wide); the four pages of fig_published.py are placed on 7.2 x 7.2 in pages.
Usage:  python assemble_pdf.py <figure folder> <output pdf>
"""
import sys
from pypdf import PdfReader, PdfWriter, PageObject, Transformation

FIG, OUT = sys.argv[1], sys.argv[2]
ORDER = ["Fig0_cartoon_DL_draft", "A_tracking_depth_and_n", "Final_Figure_Mouse", "Final_Figure_Human",
         "F1_published_thresholds_population_page", "F2_published_thresholds_inheritance_page",
         "F3_lysosome_coinheritance_page", "F4_lysohigh_lysolow_daughters_page",
         "B_pctdiff_distributions_by_birth", "C_changes_over_generations", "D_lineage_pattern_counts",
         "E_sisters_cousins_by_founder", "Part1_Fig3_mouse_summary_percent", "Part2_Fig3_human_summary_percent",
         "Part1_Fig2_mouse_by_measure_and_marker", "Part2_Fig2_human_by_measure_and_marker", "G5_check_generation4"]
SIDE, MARGIN = 7.2 * 72, 14
w = PdfWriter()
for name in ORDER:
    page = PdfReader(f"{FIG}/{name}.pdf").pages[0]
    if name.endswith("_page"):  # scale onto a square page, top-aligned
        pw, ph = float(page.mediabox.width), float(page.mediabox.height)
        s = min((SIDE - 2 * MARGIN) / pw, (SIDE - 2 * MARGIN) / ph)
        tx = (SIDE - pw * s) / 2 - float(page.mediabox.left) * s
        ty = SIDE - MARGIN - ph * s - float(page.mediabox.bottom) * s
        new = PageObject.create_blank_page(width=SIDE, height=SIDE)
        new.merge_transformed_page(page, Transformation().scale(s, s).translate(tx, ty))
        page = new
    w.add_page(page)
w.write(OUT)
print(f"{OUT}: {len(w.pages)} pages")
