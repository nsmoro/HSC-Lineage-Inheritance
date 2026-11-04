# Inheritance of asymmetric cell division along HSC lineages

Do haematopoietic stem cells (HSCs) pass on how they were born? When a cell divides, a marker (e.g. lysosomes,
active mitochondria, surface markers) can be shared equally (symmetric) or unequally (asymmetric) between the two
daughters. This repository asks whether a cell born in an asymmetric division also divides asymmetrically itself,
using curated single-cell lineage trees from long-term time-lapse imaging of mouse and human HSCs.

The code reproduces `All_figures.pdf` (17 pages) from the exported lineage tables.

## Data

The human data come from the study

> Loeffler D, Schneiter F, et al.
> **Asymmetric organelle inheritance predicts human blood stem cell fate.**
> *Blood* 2022; 139(13): 2011. https://doi.org/10.1182/blood.2020009778

- **Human:** CD49f+ HSCs (experiments 191125FS13, 191207FS14, 191216FS13, 191216FS14; HSC conditions only).
  Markers: LysoBrite (lysosomes), TMRM (active mitochondria), CD71, CD201, CD33, CD49c (integrin α3), CD34.
- **Mouse:** HSC time-lapse experiments 181117DL14 and 181219DL13 from the same group, with two antibody panels
  (Sca1 / CD41 / FcγR and CD105 / CD201 / CD71), all with LysoBrite. The mouse LysoBrite thresholds follow
  Loeffler D et al., *Nature* 2019, https://doi.org/10.1038/s41586-019-1531-6.

Only **curated** lineage trees (tracked and corrected in Matilda) are used: mouse cells of generations 0–3 and human
cells of generations 0–2 (a few generation-3 cells in one experiment are included). Raw images and tracking files
are not part of this repository.

## Repository structure

```
config/markers.yaml           channel -> marker assignment per experiment and condition
export/                       step 0: curated Matilda trees -> one cell table per experiment (needs raw data)
  export_lineage_cells.py
  io_readers.py, lineage.py   readers for Matilda .mat / TATexp.xml files and tree reconstruction
analysis/
  core.py                     divisions and mother/daughter pairs, four asymmetry measures
  dataset.py                  restriction to curated generations
  01_build_dataset.py         -> dataset.pkl
  02_birth_vs_division.py     do daughters pass on what they received? lineage enrichment vs random
  03_lineage_analyses.py      n numbers, threshold-free distributions, generation trends,
                              lineage pattern counts, sisters and cousins
  04_published_conventions.py published thresholds, LysoHigh/LysoLow daughters, per-experiment statistics
  fig_*.py                    figures (one script per figure group)
  assemble_pdf.py             -> All_figures.pdf
run_all.sh                    runs steps 1-4, the generation-4 check and all figures
```

## Usage

```bash
pip install -r requirements.txt

# Step 0 (cluster, raw data needed): export curated trees to <tables>/<experiment>_cells.csv
python export/export_lineage_cells.py --output tables/

# Steps 1-4 and all figures (about 5 minutes)
bash run_all.sh tables/ results/
```

Results: CSV tables in `results/`, figures (PDF and 300-dpi PNG) in `results/figures/`, and
`results/All_figures.pdf`. The generation-4 robustness run is written to `results/check_gen4/`.

## Methods in brief

**Asymmetry between sisters.** For every division and marker, the first value of each daughter after division
(within 4 frames in mouse, 6 frames in human) is compared with four measures: ratio (log2), % difference, absolute
difference, and absolute difference adjusted for the amount of marker in both sisters.

**Symmetric / asymmetric.**
- Main inheritance analyses: thirds of each measure within experiment × condition × marker × generation
  (most even third = symmetric, most unequal third = asymmetric; chance = 33%). This works the same for every
  measure and removes level differences between conditions and generations.
- Population plots and published-convention analyses: published thresholds where available (human LysoBrite
  <1.5× / >3×, human TMRM <1.1× / >1.5×, mouse LysoBrite <1.2× / >1.5×); other markers ≤1.1× / ≥1.8×.
- Lineage pattern counts: halves (A = more unequal than the median division, S = more equal), written A→A→A.

**Statistics.** 95% confidence intervals by resampling whole colonies; random expectations by permuting division
classes within experiment × condition × marker × generation; Mann–Whitney tests on pooled cells and paired t-tests
across independent experiments (human); mean ± SEM of experiments.

## Main results

- **Mouse:** cells born in an asymmetric division divide asymmetrically more often, for all 7 markers
  (6–20 percentage points with the ratio; LysoBrite with published thresholds 44% vs 38%).
- **Human:** no consistent effect with essentially one curated generation step; CD71 significant, small shifts
  for TMRM and CD71 without thresholds.
- **Both:** markers are co-inherited with lysosomes and active mitochondria segregate to the other daughter, as
  reported by Loeffler et al. 2022; LysoLow daughters partition CD71 more unequally in their own division.

## Author

Prepared by: Silva Morosini, Natalia  
Date: November, 2026

## Requirements

Python ≥ 3.10 with the packages in `requirements.txt`. Figures use the font Liberation Sans (falls back to the
matplotlib default if not installed).
