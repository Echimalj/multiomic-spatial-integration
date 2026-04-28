"""
Exploratory plotting utilities for spatial cell-type abundance tables.
"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu


def load_abundance_tables(
    wide_file,
    long_file,
):
    """
    Load wide and long cell-type abundance tables.

    Parameters
    ----------
    wide_file : str
        Path to wide abundance table.
    long_file : str
        Path to long abundance table.

    Returns
    -------
    tuple
        wide_df, long_df
    """
    wide_df = pd.read_csv(wide_file)
    long_df = pd.read_csv(long_file)

    return wide_df, long_df


def summarize_roi_counts(
    df_long,
    roi_col="ROI",
    group_col="disease_status",
):
    """
    Count unique ROIs per group.
    """
    return (
        df_long
        .groupby(group_col)[roi_col]
        .nunique()
        .reset_index(name="n_rois")
    )


def plot_celltype_by_group(
    df_long,
    celltype,
    group="disease_status",
    abundance_col="rel_abundance",
    hue=None,
    title=None,
    figsize=(6, 4),
    save_file=None,
):
    """
    Plot one cell type by group using boxplot + stripplot.
    """
    tmp = df_long[df_long["celltype"] == celltype].copy()

    plt.figure(figsize=figsize)

    sns.boxplot(
        data=tmp,
        x=group,
        y=abundance_col,
        hue=hue
    )

    sns.stripplot(
        data=tmp,
        x=group,
        y=abundance_col,
        hue=hue,
        dodge=True if hue is not None else False,
        alpha=0.4,
        color="black"
    )

    plt.title(title or f"{celltype} by {group}")
    plt.tight_layout()

    if save_file is not None:
        Path(save_file).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_file, dpi=300, bbox_inches="tight")

    plt.show()


def plot_celltype_by_pathology_and_disease(
    df_long,
    celltype,
    abundance_col="rel_abundance",
    save_file=None,
):
    """
    Plot one cell type by pathology, faceted by disease status.
    """
    tmp = df_long[df_long["celltype"] == celltype].copy()

    g = sns.catplot(
        data=tmp,
        x="pathology",
        y=abundance_col,
        col="disease_status",
        kind="box",
        height=4,
        aspect=1
    )

    g.fig.suptitle(f"{celltype} by pathology", y=1.05)

    if save_file is not None:
        Path(save_file).parent.mkdir(parents=True, exist_ok=True)
        g.savefig(save_file, dpi=300, bbox_inches="tight")

    plt.show()


def plot_celltype_by_region_and_disease(
    df_long,
    celltype,
    abundance_col="rel_abundance",
    save_file=None,
):
    """
    Plot one cell type by region and disease status.
    """
    tmp = df_long[df_long["celltype"] == celltype].copy()

    plt.figure(figsize=(8, 4))

    sns.boxplot(
        data=tmp,
        x="region",
        y=abundance_col,
        hue="disease_status"
    )

    plt.xticks(rotation=30, ha="right")
    plt.title(f"{celltype} by region and disease status")
    plt.tight_layout()

    if save_file is not None:
        Path(save_file).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_file, dpi=300, bbox_inches="tight")

    plt.show()


def plot_mean_composition_stacked(
    df_long,
    group_col="disease_status",
    abundance_col="rel_abundance",
    save_file=None,
):
    """
    Plot mean cell-type composition as stacked barplot.
    """
    mean_mix = (
        df_long
        .groupby([group_col, "celltype"])[abundance_col]
        .mean()
        .reset_index()
    )

    mix_piv = (
        mean_mix
        .pivot(
            index=group_col,
            columns="celltype",
            values=abundance_col
        )
        .fillna(0)
    )

    ax = mix_piv.plot(
        kind="bar",
        stacked=True,
        figsize=(10, 5)
    )

    ax.set_ylabel("Mean relative abundance")
    ax.set_title(f"Mean cell-type composition by {group_col}")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")

    plt.tight_layout()

    if save_file is not None:
        Path(save_file).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_file, dpi=300, bbox_inches="tight")

    plt.show()

    return mix_piv


def run_mannwhitney_screen(
    df_long,
    group_col="disease_status",
    group_1="AD-CAA",
    group_2="Control",
    abundance_col="rel_abundance",
    min_n=5,
):
    """
    Quick non-parametric screening test by cell type.

    This is exploratory only. Formal inference should use mixed-effects models.
    """
    stats = []

    for ct in df_long["celltype"].unique():
        tmp = df_long[df_long["celltype"] == ct]

        g1 = tmp[tmp[group_col] == group_1][abundance_col]
        g2 = tmp[tmp[group_col] == group_2][abundance_col]

        if len(g1) > min_n and len(g2) > min_n:
            u, p = mannwhitneyu(g1, g2, alternative="two-sided")
            stats.append([ct, p, g1.mean(), g2.mean(), len(g1), len(g2)])

    stats_df = pd.DataFrame(
        stats,
        columns=[
            "celltype",
            "p_value",
            f"mean_{group_1}",
            f"mean_{group_2}",
            f"n_{group_1}",
            f"n_{group_2}"
        ]
    )

    stats_df["q_BH"] = stats_df["p_value"].rank(method="first")
    stats_df["q_BH"] = stats_df["p_value"] * len(stats_df) / stats_df["q_BH"]
    stats_df["q_BH"] = stats_df["q_BH"].clip(upper=1)

    return stats_df.sort_values("p_value")
