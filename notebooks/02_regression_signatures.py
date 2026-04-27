import sys
import torch
import scanpy as sc

from python.regression_utils import (
    load_seurat_mtx_as_anndata,
    standardize_regression_obs,
    round_counts_layer,
    minimal_filter_anndata,
    disable_lightning_mpi_detection,
    train_regression_models_by_condition,
    save_regression_models,
    load_and_export_regression_posterior,
    export_inferred_signatures,
)

print(sys.executable)
print(torch.cuda.is_available())

data_path = "/N/u/echimal/Quartz/Desktop/CLR_MRI/Human_GeoMx_Sep2025/"
rresults = f"{data_path}/Regression-model/"

adata = load_seurat_mtx_as_anndata(
    counts_mtx=f"{data_path}/AD_CAA_counts.mtx",
    metadata_csv=f"{data_path}/AD_CAA_meta.csv",
    genes_csv=f"{data_path}/AD_CAA_genes.csv",
    gene_col="gene_id",
    transpose_counts=True,
    counts_layer="counts",
)

adata = standardize_regression_obs(
    adata,
    celltype_col="New_Idents",
    new_celltype_col="celltype",
    experiment_col="Experiment",
    experiment_value="Batch1",
)

adata = round_counts_layer(
    adata,
    counts_layer="counts",
    replace_x=True
)

adata = minimal_filter_anndata(
    adata,
    min_genes=1,
    min_cells=1,
    counts_layer="counts"
)

adata

disable_lightning_mpi_detection()

trained_models = train_regression_models_by_condition(
    adata,
    condition_col="FDX",
    conditions=["AD+CAA", "Control"],
    labels_key="celltype",
    batch_key="Experiment",
    layer="counts",
    max_epochs=100,
    batch_size=1024,
    lr=0.01,
    accelerator="gpu"
)

save_regression_models(
    trained_models,
    output_dir=rresults
)

adata_AD = adata[adata.obs["FDX"] == "AD+CAA"].copy()
adata_CTRL = adata[adata.obs["FDX"] == "Control"].copy()

adata_AD, mod_AD = load_and_export_regression_posterior(
    adata_AD,
    model_dir=f"{rresults}/AD+CAA_regression_model/",
    num_samples=1000,
    batch_size=2500
)

adata_CTRL, mod_CTRL = load_and_export_regression_posterior(
    adata_CTRL,
    model_dir=f"{rresults}/Control_regression_model/",
    num_samples=1000,
    batch_size=2500
)

inf_AD, avg_AD = export_inferred_signatures(
    adata_AD,
    covariate_col="celltype",
    model_name="AD+CAA",
    output_folder=rresults
)

inf_CTRL, avg_CTRL = export_inferred_signatures(
    adata_CTRL,
    covariate_col="celltype",
    model_name="Control",
    output_folder=rresults
)
