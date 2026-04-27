"""
Utilities for running SpaceJam / Cell2Location-style Pyro models.
"""

import torch
import numpy as np
import pandas as pd


def prepare_spacejam_inputs(
    adata,
    signatures,
    neg_probes_key="negProbes",
    counts_layer="counts",
    nuclei_col="Nuclei"
):
    """
    Align AnnData and signature matrix for spatial modeling.
    """

    counts = adata.layers[counts_layer]
    if not isinstance(counts, np.ndarray):
        counts = counts.toarray()

    neg = adata.obsm[neg_probes_key]

    genes_shared = adata.var_names.intersection(signatures.index)

    counts = counts[:, adata.var_names.isin(genes_shared)]
    signatures = signatures.loc[genes_shared]

    return {
        "counts": torch.tensor(counts, dtype=torch.float32),
        "signatures": torch.tensor(signatures.values, dtype=torch.float32),
        "neg_probes": torch.tensor(neg, dtype=torch.float32),
        "nuclei": torch.tensor(adata.obs[nuclei_col].values, dtype=torch.float32)
    }


def train_spacejam_model(
    inputs,
    model_class,
    n_steps=20000,
    learning_rate=0.005,
    device="cpu"
):
    """
    Train Pyro spatial model.
    """

    model = model_class(**inputs)

    model.train_model(
        n_steps=n_steps,
        lr=learning_rate,
        device=device
    )

    return model


def extract_posterior_abundance(
    model,
    adata,
    normalize=True
):
    """
    Extract cell-type abundance from trained model.
    """

    abundance = model.get_cell_abundance()

    df = pd.DataFrame(
        abundance,
        index=adata.obs_names,
        columns=model.cell_types
    )

    if normalize:
        df = df.div(df.sum(axis=1), axis=0)

    return df
