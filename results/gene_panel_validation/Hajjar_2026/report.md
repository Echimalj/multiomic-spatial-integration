# Hajjar 2026 biomarker panel

_Generated: 2026-08-03 14:54 EDT_

## Panel overview

- Panel identifier: `Hajjar_2026`
- Genes analyzed: 40
- Genes missing from the reference signatures: 1
- Headline genes: 2

## Evidence summary

- 27 of 40 genes had at least one FDR-significant spatial perturbation.
- No lineage interaction passed global FDR; 15 genes had nominal exploratory lineage evidence.
- 34 of 40 genes had at least one globally FDR-significant fine-cell-type association.

## Ranking framework

The evidence score is a transparent prioritization tool:

- Spatial FDR evidence: 3 points
- Spatial nominal-only evidence: 1 point
- Fine-cell-type global FDR evidence: 3 points
- Lineage global FDR evidence: 3 points
- Exploratory nominal lineage evidence: 1 point
- Panel-designated headline target: 1 point

The score is intended for prioritization and does not represent an independent statistical test.

## Top 15 ranked genes

Rank | Gene | Tier | Score | Reference_lineage | Spatial_region | Spatial_contrast | Spatial_effect | Spatial_FDR | Fine_celltype_AmyloidFree | Fine_celltype_Amyloid
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
1 | MAPT | Tier A: convergent evidence | 7 | ExcitatoryNeuron | Parenchyma | MaxPathology_effect | 1.59 | 0.000000433 | ExNeuron3 | ExNeuron3
2 | SOD1 | Tier A: convergent evidence | 7 | ExcitatoryNeuron | Capillaries | Amyloid_effect | -1.34 | 0.0000378 | OPC3 | ExNeuron10
3 | NTRK2 | Tier A: convergent evidence | 7 | InhibitoryNeuron | Capillaries | Amyloid_effect | -1.34 | 0.000892 | Astrocytes1 | Astrocytes2
4 | NGF | Tier A: convergent evidence | 7 | InhibitoryNeuron | Capillaries | Amyloid_effect | 0.940 | 0.00135 | Oligodendrocytes6 | OPC1
5 | LGALS8 | Tier A: convergent evidence | 7 | ExcitatoryNeuron | Capillaries | Amyloid_effect | 0.867 | 0.00142 | Oligodendrocytes3 | InhNeuron1
6 | TNFRSF12A | Tier A: convergent evidence | 7 | Vascular | Capillaries | Amyloid_effect | 0.810 | 0.00440 | Fibroblast | OPC3
7 | NRP2 | Tier A: convergent evidence | 7 | Vascular | Capillaries | Amyloid_effect | 0.746 | 0.00763 | OPC3 | ExNeuron8
8 | SPARCL1 | Tier A: convergent evidence | 7 | ExcitatoryNeuron | Capillaries | Amyloid_effect | -1.48 | 0.0102 | Astrocytes2 | Astrocytes2
9 | NCAM1 | Tier A: convergent evidence | 7 | ExcitatoryNeuron | Capillaries | Amyloid_effect | -0.891 | 0.0135 | Fibroblast | Fibroblast
10 | PTPRS | Tier A: convergent evidence | 7 | InhibitoryNeuron | Capillaries | Disease_effect | 0.911 | 0.0168 | ExNeuron10 | Microglia5
11 | VEGFA | Tier A: convergent evidence | 7 | Astrocyte | Capillaries | Amyloid_effect | 0.568 | 0.0260 | InhNeuron4 | SMC
12 | JAM2 | Tier A: convergent evidence | 7 | InhibitoryNeuron | Parenchyma | MaxPathology_effect | 0.674 | 0.0411 | Oligodendrocytes3 | InhNeuron3
13 | PLAU | Tier A: convergent evidence | 7 | Vascular | Capillaries | Amyloid_effect | 0.588 | 0.0413 | Astrocytes1 | SMC
14 | ICAM1 | Tier A: convergent evidence | 7 | Vascular | Capillaries | Amyloid_effect | 1.01 | 0.000351 | InhNeuron2 | ExNeuron2
15 | VCAM1 | Tier A: convergent evidence | 7 | Vascular | Capillaries | Amyloid_effect | 0.839 | 0.00763 | ExNeuron2 | Astrocytes4

## Figures

### Reference lineage attribution

![Reference lineage attribution](../../figures/gene_panel_validation/Hajjar_2026/gene_panel_reference_lineage_heatmap.png)

### Spatial perturbation

![Spatial perturbation](../../figures/gene_panel_validation/Hajjar_2026/gene_panel_spatial_perturbation_heatmap.png)

### Lineage-associated perturbation — exploratory nominal display

![Lineage-associated perturbation — exploratory nominal display](../../figures/gene_panel_validation/Hajjar_2026/gene_panel_lineage_attribution_heatmap_exploratory_nominal.png)

### Fine-cell-type discovery correlations

![Fine-cell-type discovery correlations](../../figures/gene_panel_validation/Hajjar_2026/gene_panel_fine_celltype_correlation_heatmap.png)

## Output files

- Integrated summary: `results/gene_panel_validation/Hajjar_2026/integrated_gene_panel_summary.csv`
- Ranked priorities: `results/gene_panel_validation/Hajjar_2026/ranked_gene_priorities.csv`

## Interpretation note

Reference attribution indicates the lineage in which a gene is preferentially expressed. Spatial perturbation tests region-specific expression changes. Lineage interactions assess whether the expression–lineage association changes with local amyloid pathology. Fine-cell-type correlations are scan-aware associations and should not be interpreted as direct cellular causation.
