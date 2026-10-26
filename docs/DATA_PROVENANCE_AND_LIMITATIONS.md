# Data provenance and interpretation limits

## Source

Values in `data/` were transcribed from the supplied 38-page PDF, `enzymology_1793005839159.pdf`. The source PDF is not copied into this package. Page numbers in the CSV files point to the corresponding source tables.

The report contains one listed absorbance per condition. No replicate-level readings, standard deviations, raw spectrophotometer export, or plate/tube-level metadata were supplied.

## Calibration check

The nine maltose standards on page 8 were refit by ordinary least squares. The data give:

`A540 = 0.81833 × C − 0.01578`, `R² = 0.99070`

where `C` is maltose concentration in mg/mL. This agrees with the negative-intercept equation printed in the graph and later calculation, but the source also prints a positive-intercept version. The discrepancy is recorded here rather than silently choosing the inconsistent text.

The source's reported concentration tables are retained in `assay_results.csv`. The generated `results/recalculated_product_yield.csv` uses the fitted standard curve and raw absorbances. Several ethanol-table values differ from recalculation by more than rounding; use absorbance as the primary reported measurement and treat all concentration estimates as maltose-equivalent estimates.

## Experimental design limits

- There are no reported replicates, error estimates, or inferential statistical tests. Plots therefore show points and connecting lines without error bars or significance claims.
- DNS detects reducing sugars, not maltose specifically. Concentrations calculated from a maltose calibration should be described as maltose equivalents.
- The endpoint is absorbance after a fixed incubation, not a normalized enzyme rate such as U/mL or U/mg protein.
- The source compares 1% commercial fungal amylase with diluted saliva, without a protein- or activity-normalization step. Their absolute signals are not directly comparable as intrinsic catalytic efficiencies.
- pH 3 and 5 use sodium acetate buffer, while pH 7, 9, and 12 use potassium phosphate. Buffer chemistry therefore changes with pH.
- Calcium and ethanol concentration labels may refer to prepared solutions rather than the final reaction mixture; the source does not clearly state final concentrations after all additions.
- Starch signal rises through the highest tested concentration (0.5%). A saturation plateau is not demonstrated, so Michaelis–Menten constants cannot be estimated from this range.
- The PDF gives no raw instrument files or independent calibration replicates, so transcription and assay quality cannot be independently verified.

## Reading the figures

Figures 2–6 plot the source-reported absorbance at 540 nm, the direct tabulated endpoint. Figure 1 uses the tabulated standard points and a regression calculated by the included script. Recalculated maltose-equivalent estimates are provided separately and are not substituted for the source measurements.
