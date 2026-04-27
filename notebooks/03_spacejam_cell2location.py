# 03 — Bayesian Spatial Deconvolution (SpaceJam / Cell2Location)

#This notebook performs spatial deconvolution of GeoMx WTA data using:

#- snRNA-seq-derived cell-type signatures
#- a modified Pyro-based Cell2Location/SpaceJam model

#The output is a cell-type abundance estimate for each ROI.

import torch
import scanpy as sc
import pandas as pd

from models.LocationModelWTAMultiExperimentHierarchicalGeneLevel_Modified import LocationModelPyro
from python.spacejam_pyro_utils import (
    prepare_spacejam_inputs,
    train_spacejam_model,
    extract_posterior_abundance
)

## Load AnnData (GeoMx WTA)
adata_wta = sc.read_h5ad("data/CAA-AD_AnnData.h5ad")
adata_wta

## Load regression-derived cell-type signatures
inf_aver = pd.read_csv(
    "data/Regression-model/AD+CAA_inferred_signatures.csv",
    index_col=0
)

inf_aver.head()

## Prepare inputs for SpaceJam model

#This step aligns:

#- WTA counts
#- negative probes
#- cell-type signatures
#- ROI metadata

inputs = prepare_spacejam_inputs(
    adata=adata_wta,
    signatures=inf_aver,
    neg_probes_key="negProbes",
    counts_layer="counts",
    nuclei_col="Nuclei"
)

## Train SpaceJam / Pyro model
model = train_spacejam_model(
    inputs,
    model_class=LocationModelPyro,
    n_steps=20000,
    learning_rate=0.005,
    device="cuda"
)

## Extract posterior cell-type abundances
cell_abundance = extract_posterior_abundance(
    model,
    adata=adata_wta,
    normalize=True
)

cell_abundance.head()

## Save results
cell_abundance.to_csv("results/roi_celltype_abundance.csv")
