# AGORA high-nomination targets

_Generated: 2026-08-03 14:54 EDT_

## Panel overview

- Panel identifier: `AGORA_High_Nomination`
- Genes analyzed: 38
- Genes missing from the reference signatures: 0
- Headline genes: 14

## Evidence summary

- 13 of 38 genes had at least one FDR-significant spatial perturbation.
- No lineage interaction passed global FDR; 14 genes had nominal exploratory lineage evidence.
- 37 of 38 genes had at least one globally FDR-significant fine-cell-type association.

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
1 | CLU | Tier A: convergent evidence | 8 | Astrocyte | Parenchyma | MaxPathology_effect | 2.26 | 0.000000674 | Astrocytes2 | Astrocytes2
2 | PLEC | Tier A: convergent evidence | 8 | ExcitatoryNeuron | Capillaries | Disease_effect | 1.31 | 0.0000620 | ExNeuron1 | OPC3
3 | PRDX1 | Tier A: convergent evidence | 8 | Astrocyte | Capillaries | Amyloid_effect | -1.16 | 0.00374 | Oligodendrocytes5 | Fibroblast
4 | HTRA1 | Tier A: convergent evidence | 8 | OPC | Capillaries | Amyloid_effect | -0.954 | 0.0386 | Astrocytes2 | Astrocytes2
5 | NDUFA10 | Tier A: convergent evidence | 7 | ExcitatoryNeuron | Capillaries | Amyloid_effect | -0.968 | 0.00433 | Microglia5 | OPC1
6 | DBI | Tier A: convergent evidence | 7 | Astrocyte | Capillaries | Amyloid_effect | -0.846 | 0.00594 | Astrocytes5 | ExNeuron9
7 | SEPTIN5 | Tier A: convergent evidence | 7 | ExcitatoryNeuron | Parenchyma | MaxPathology_effect | 1.19 | 0.0221 | ExNeuron3 | ExNeuron3
8 | INPP5D | Tier A: convergent evidence | 7 | Microglia | Capillaries | Amyloid_effect | 0.842 | 0.00184 | Microglia1 | VLMC1
9 | BIN1 | Tier A: convergent evidence | 7 | Microglia | Arteries | Amyloid_effect | -0.862 | 0.00305 | SMC | Oligodendrocytes7
10 | APOE | Tier A: convergent evidence | 7 | Astrocyte | Capillaries | Amyloid_effect | -1.09 | 0.00521 | Oligodendrocytes4 | InhNeuron5
11 | PDHB | Tier B: strong evidence | 6 | ExcitatoryNeuron | Parenchyma | Disease_effect | -0.499 | 0.340 | ExNeuron3 | ExNeuron9
12 | RUFY3 | Tier B: strong evidence | 6 | ExcitatoryNeuron | Capillaries | Amyloid_effect | -1.13 | 0.00184 | Microglia5 | Fibroblast
13 | MDK | Tier B: strong evidence | 6 | Vascular | Arteries | MaxPathology_effect | 0.783 | 0.00925 | ExNeuron12 | ExNeuron1
14 | SYK | Tier B: strong evidence | 6 | Microglia | Capillaries | Amyloid_effect | 0.664 | 0.0386 | VLMC1 | VLMC1
15 | PTN | Tier B: strong evidence | 5 | Astrocyte | Capillaries | Amyloid_effect | -0.652 | 0.0694 | Astrocytes4 | Astrocytes1

## Figures

### Reference lineage attribution

![Reference lineage attribution](../../figures/gene_panel_validation/AGORA_High_Nomination/gene_panel_reference_lineage_heatmap.png)

### Spatial perturbation

![Spatial perturbation](../../figures/gene_panel_validation/AGORA_High_Nomination/gene_panel_spatial_perturbation_heatmap.png)

### Lineage-associated perturbation — exploratory nominal display

![Lineage-associated perturbation — exploratory nominal display](../../figures/gene_panel_validation/AGORA_High_Nomination/gene_panel_lineage_attribution_heatmap_exploratory_nominal.png)

### Fine-cell-type discovery correlations

![Fine-cell-type discovery correlations](../../figures/gene_panel_validation/AGORA_High_Nomination/gene_panel_fine_celltype_correlation_heatmap.png)

## Output files

- Integrated summary: `results/gene_panel_validation/AGORA_High_Nomination/integrated_gene_panel_summary.csv`
- Ranked priorities: `results/gene_panel_validation/AGORA_High_Nomination/ranked_gene_priorities.csv`

## Interpretation note

Reference attribution indicates the lineage in which a gene is preferentially expressed. Spatial perturbation tests region-specific expression changes. Lineage interactions assess whether the expression–lineage association changes with local amyloid pathology. Fine-cell-type correlations are scan-aware associations and should not be interpreted as direct cellular causation.
