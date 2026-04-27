# Requires: pyro-ppl, torch, numpy, pandas
# pip install pyro-ppl torch

import torch
import torch.nn as nn
import pyro
import pyro.distributions as dist
from pyro.infer import SVI, Trace_ELBO, Predictive
from pyro.infer.autoguide import AutoNormal
from pyro.optim import Adam
import numpy as np

pyro.set_rng_seed(0)


class LocationModelPyro:
    """
    Pyro reimplementation of LocationModelWTAMultiExperimentHierarchicalGeneLevel.
    Inputs expected as numpy arrays (will be converted to torch tensors).
    """

    def __init__(
        self,
        cell_state_mat,   # genes x n_factors (numpy)
        X_data,           # n_rois x n_genes (numpy)  -> gene probes
        Y_data,           # n_rois x n_npro  (numpy)  -> negative probes
        spot2sample_mat,  # n_rois x n_exper (numpy)
        device=None,
        dtype=torch.float32,
        gene_level_prior={"mean": 1 / 2, "sd": 1 / 4, "sample_alpha": 20},
        cell_number_prior={"cells_per_spot": 8, "factors_per_spot": 7, "combs_per_spot": 2.5},
        n_comb=50,
        phi_hyp_prior={"mean": 3, "sd": 1},
        spot_fact_mean_var_ratio=0.5,
    ):
        # device
        self.device = device or ("cuda" if torch.cuda.is_available() else "cpu")
        self.dtype = dtype

        # convert arrays to torch tensors on device
        self.cell_state = torch.tensor(cell_state_mat, dtype=self.dtype, device=self.device)  # genes x n_fact
        self.X_data = torch.tensor(X_data, dtype=self.dtype, device=self.device)            # n_rois x n_genes
        self.Y_data = torch.tensor(Y_data, dtype=self.dtype, device=self.device)            # n_rois x n_npro
        self.spot2sample = torch.tensor(spot2sample_mat, dtype=self.dtype, device=self.device)  # n_rois x n_exper

        # dims
        self.n_rois = self.X_data.shape[0]
        self.n_genes = self.X_data.shape[1]
        self.n_npro = self.Y_data.shape[1]
        self.n_exper = self.spot2sample.shape[1]
        self.n_fact = self.cell_state.shape[1]
        self.n_comb = n_comb

        # total number of gene counts per ROI divided by 1e5 (as in original)
        l_r_np = np.array([np.sum(X_data[i, :]) for i in range(self.n_rois)]).reshape(self.n_rois, 1) * 1e-5
        self.l_r = torch.tensor(l_r_np, dtype=self.dtype, device=self.device)

        # save priors
        self.gene_level_prior = gene_level_prior.copy()
        self.cell_number_prior = cell_number_prior.copy()
        self.phi_hyp_prior = phi_hyp_prior.copy()
        self.spot_fact_mean_var_ratio = spot_fact_mean_var_ratio

        # model + guide placeholders
        self.guide = None
        self.svi = None
        self.optim = None

        # set large alpha for negative probes (approx Poisson)
        self.big_alpha = 1e10

    def model(self):
        """Pyro model using the exact hierarchical structure where possible."""
        # ------------------ Negative probe binding --------------------
        # b_n_hyper: Gamma(3,1) and Gamma(1,1) shape=2 in original -> treat as two hyperparams
        b_n_hyper = pyro.sample("b_n_hyper", dist.Gamma(torch.tensor(3.0, device=self.device), torch.tensor(1.0, device=self.device)).expand([2]).to_event(1))
        # b_n: shape (n_exper, n_npro) Gamma with mean b_n_hyper[0] and sd b_n_hyper[1]
        # Pyro does not have Gamma(mu, sigma) direct; we use alpha,beta: mean = alpha/beta.
        # Convert mu,sigma -> alpha,beta
        # For simplicity, assume alpha = (mu^2)/(sigma^2), beta = mu/(sigma^2)
        # When hyperparams are scalars, we broadcast.
        mu_b = b_n_hyper[0]
        sd_b = b_n_hyper[1]
        alpha_b = (mu_b ** 2) / (sd_b ** 2 + 1e-9)
        beta_b = mu_b / (sd_b ** 2 + 1e-9)
        b_n = pyro.sample("b_n", dist.Gamma(alpha_b.expand([self.n_exper, self.n_npro]), beta_b.expand([self.n_exper, self.n_npro])).to_event(2))
        # y_rn = spot2sample @ b_n * l_r  -> shape n_rois x n_npro
        y_rn = torch.matmul(self.spot2sample, b_n) * self.l_r  # broadcasting l_r

        # ------------------ gene_add (non-specific binding gene probes) ------------------
        # use same hyperparams: draw gene_add per (n_exper, n_genes)
        gene_add = pyro.sample("gene_add", dist.Gamma(alpha_b.expand([self.n_exper, self.n_genes]), beta_b.expand([self.n_exper, self.n_genes])).to_event(2))

        # ------------------ Gene-level scaling priors ------------------
        # compute shape/rate from mean/sd as in original
        mean_gl = torch.tensor(self.gene_level_prior["mean"], dtype=self.dtype, device=self.device)
        sd_gl = torch.tensor(self.gene_level_prior["sd"], dtype=self.dtype, device=self.device)
        shape = (mean_gl ** 2) / (sd_gl ** 2 + 1e-9)
        rate = mean_gl / (sd_gl ** 2 + 1e-9)
        # gene_level_alpha_hyp / gene_level_beta_hyp as Gamma's (small shape arrays)
        gene_level_alpha_hyp = pyro.sample(
            "gene_level_alpha_hyp",
            dist.Gamma(
                shape,
                torch.sqrt(shape / (self.gene_level_prior.get("mean_var_ratio", 1.0)) + 1e-9)
            )
        )

        gene_level_beta_hyp = pyro.sample(
            "gene_level_beta_hyp",
            dist.Gamma(
                rate,
                torch.sqrt(rate / (self.gene_level_prior.get("mean_var_ratio", 1.0)) + 1e-9)
            )
        )
        # gene_level: Gamma(shape=gene_level_alpha_hyp, rate=gene_level_beta_hyp) shape (1, n_genes)
        # convert gene_level_alpha_hyp/gene_level_beta_hyp to alpha/beta for Gamma
        # They were used as mu and sigma in PyMC; here approximate by drawing gamma with alpha/beta such that mean ~ mu
        # For simplicity, draw gene_level around mean_gl with some dispersion:
        gene_level = pyro.sample("gene_level", dist.Gamma((mean_gl ** 2) / (sd_gl ** 2 + 1e-9), mean_gl / (sd_gl ** 2 + 1e-9)).expand([1, self.n_genes]).to_event(2))

        # gene_level_independent: narrow Gamma(100,100): mean=1
        gene_level_independent = pyro.sample("gene_level_independent", dist.Gamma(torch.tensor(100.0, device=self.device), torch.tensor(100.0, device=self.device)).expand([self.n_exper, self.n_genes]).to_event(2))

        # gene_level_e: per experiment capture efficiency
        sample_alpha = float(self.gene_level_prior.get("sample_alpha", 20.0))

        gene_level_e = pyro.sample(
            "gene_level_e",
            dist.Gamma(
                torch.tensor(sample_alpha, dtype=self.dtype, device=self.device),
                torch.tensor(sample_alpha, dtype=self.dtype, device=self.device),
            )
            .expand([self.n_exper, 1])
            .to_event(2)
        )


        # gene_factors deterministic: cell_state (genes x n_fact) transposed to n_fact x genes
        gene_factors = self.cell_state.T  # n_fact x n_genes

        # ------------------ Spot factors priors ------------------
        cells_per_spot_mean = float(self.cell_number_prior["cells_per_spot"])
        cells_mvr = float(self.cell_number_prior.get("cells_mean_var_ratio", 1.0))

        cells_per_spot = pyro.sample(
            "cells_per_spot",
            dist.Gamma(
                torch.tensor(cells_per_spot_mean, dtype=self.dtype, device=self.device),
                torch.sqrt(
                    torch.tensor(cells_per_spot_mean / cells_mvr, dtype=self.dtype, device=self.device)
                ),
            )
            .expand([self.n_rois, 1])
            .to_event(2)
        )


        combs_per_spot_mean = float(self.cell_number_prior["combs_per_spot"])
        combs_mvr = float(self.cell_number_prior.get("combs_mean_var_ratio", 1.0))

        combs_per_spot = pyro.sample(
            "combs_per_spot",
            dist.Gamma(
                torch.tensor(combs_per_spot_mean, dtype=self.dtype, device=self.device),
                torch.sqrt(
                    torch.tensor(combs_per_spot_mean / combs_mvr, dtype=self.dtype, device=self.device)
                ),
            )
            .expand([self.n_rois, 1])
            .to_event(2)
        )

        # combs_factors ~ Gamma(shape = combs_per_spot / n_comb, rate = (1/cells_per_spot)*combs_per_spot)
        shape_cf = combs_per_spot / float(self.n_comb)
        rate_cf = (1.0 / (cells_per_spot + 1e-9)) * combs_per_spot
        combs_factors = pyro.sample("combs_factors", dist.Gamma(shape_cf, rate_cf).expand([self.n_rois, self.n_comb]).to_event(2))

        # factors_per_combs
        fpc_mean = float(self.cell_number_prior["factors_per_combs"])
        fpc_mvr = float(self.cell_number_prior.get("factors_mean_var_ratio", 1.0))

        factors_per_combs = pyro.sample(
            "factors_per_combs",
            dist.Gamma(
                torch.tensor(fpc_mean, dtype=self.dtype, device=self.device),
                torch.sqrt(
                    torch.tensor(fpc_mean / fpc_mvr, dtype=self.dtype, device=self.device)
                ),
            )
            .expand([self.n_comb, 1])
            .to_event(2)
        )


        c2f_shape = factors_per_combs / float(self.n_fact)
        comb2fact = pyro.sample("comb2fact", dist.Gamma(c2f_shape, factors_per_combs.expand([self.n_comb, self.n_fact])).to_event(2))

        # spot_factors: Gamma with mean = combs_factors @ comb2fact, sigma = sqrt(mean / spot_fact_mean_var_ratio)
        eps = 1e-6  # numerical safety

        mean_spot_factors = torch.matmul(combs_factors, comb2fact)
        mean_spot_factors = torch.clamp(mean_spot_factors, min=eps)

        sigma_spot_factors = torch.sqrt(
            mean_spot_factors / (self.spot_fact_mean_var_ratio + eps)
        )
        sigma_spot_factors = torch.clamp(sigma_spot_factors, min=eps)

        concentration = mean_spot_factors ** 2 / (sigma_spot_factors ** 2 + eps)
        rate = mean_spot_factors / (sigma_spot_factors ** 2 + eps)

        concentration = torch.clamp(concentration, min=eps)
        rate = torch.clamp(rate, min=eps)

        spot_factors = pyro.sample(
            "spot_factors",
            dist.Gamma(concentration, rate).to_event(2)
        )


        # spot_additive
        spot_add_hyp = pyro.sample(
            "spot_add_hyp",
            dist.Gamma(
                torch.tensor(1.0, dtype=self.dtype, device=self.device),
                torch.tensor(1.0, dtype=self.dtype, device=self.device),
            )
            .expand([2])
            .to_event(1)
        )

        spot_add = pyro.sample("spot_add", dist.Gamma(spot_add_hyp[0], spot_add_hyp[1]).expand([self.n_rois, 1]).to_event(2))

        # gene overdispersion phi_hyp and gene_E ~ Exponential(phi_hyp)
        phi_mean = float(self.phi_hyp_prior["mean"])
        phi_sd = float(self.phi_hyp_prior["sd"])

        phi_hyp = pyro.sample(
            "phi_hyp",
            dist.Gamma(
                torch.tensor(phi_mean, dtype=self.dtype, device=self.device),
                torch.tensor(1.0 / (phi_sd + 1e-9), dtype=self.dtype, device=self.device),
            )
            .expand([1, 1])
            .to_event(1)
        )

        # gene_E: Exponential(phi_hyp) per experiment x gene
        gene_E = pyro.sample("gene_E", dist.Exponential(phi_hyp.expand([self.n_exper, self.n_genes])).to_event(2))

        # ------------------ Expected expression -----------------------
        # mu_biol = spot_factors @ gene_factors.T * gene_level * (spot2sample @ gene_level_e) * (spot2sample @ gene_level_independent)
        #          + (spot2sample @ gene_add) * l_r + spot_add
        # shapes: spot_factors (n_rois x n_fact), gene_factors (n_fact x n_genes)
        mu_biol = torch.matmul(spot_factors, gene_factors) \
                  * gene_level.expand([self.n_rois, self.n_genes]) \
                  * torch.matmul(self.spot2sample, gene_level_e).expand([self.n_rois, self.n_genes]) \
                  * torch.matmul(self.spot2sample, gene_level_independent) \
                  + torch.matmul(self.spot2sample, gene_add) * self.l_r \
                  + spot_add

        # concatenate negative probes expectation y_rn with mu_biol along genes axis
        mu_concat = torch.cat([y_rn, mu_biol], dim=1)  # n_rois x (n_npro + n_genes)

        # compute alpha (overdispersion) concat: for negatives large alpha, for genes spot2sample @ (1/(gene_E^2))
        alpha_genes = torch.matmul(self.spot2sample, 1.0 / (gene_E * gene_E + 1e-9))  # n_rois x n_genes
        alpha_concat = torch.cat([torch.full_like(y_rn, fill_value=self.big_alpha), alpha_genes], dim=1)  # n_rois x total_features

        # ------------------ Likelihood (Negative Binomial via logits) ------------------
        #change of NB to logits ECJ 06Jan2026
        eps = 1e-6

        mu = mu_concat.clamp(min=eps)
        r = alpha_concat.clamp(min=eps)

        logits = torch.log(mu_concat + 1e-8) - torch.log(r + 1e-8)

        obs_concat = torch.cat([self.Y_data, self.X_data], dim=1)

        pyro.sample(
            "data_target",
            dist.NegativeBinomial(
                total_count=r,
                logits=logits
            ).to_event(2),
            obs=obs_concat
        )

        # store some deterministic values for later retrieval via Predictive if desired
        pyro.deterministic("mu_biol", mu_biol)
        pyro.deterministic("mu_concat", mu_concat)

        factor_gene_load = (gene_factors * gene_level).sum(dim=1)
        nUMI_factors = (spot_factors * factor_gene_load.unsqueeze(0)).sum(dim=1)

        pyro.deterministic("nUMI_factors", nUMI_factors)


    def fit(self, n_steps=5000, lr=1e-3, warmup=0):
        """
        Fit the model with SVI using an AutoNormal guide.
        Returns the trained svi object and history.
        """
        pyro.clear_param_store()
        # Auto guide
        self.guide = AutoNormal(self.model)
        self.optim = Adam({"lr": lr})
        self.svi = SVI(self.model, self.guide, self.optim, loss=Trace_ELBO())

        loss_hist = []
        for step in range(n_steps):
            loss = self.svi.step()
            loss_hist.append(loss)
            if step % max(1, n_steps // 10) == 0:
                print(f"Step {step}/{n_steps} Loss = {loss:.1f}")
        return loss_hist

    def posterior_predictive(self, num_samples=100):
        """
        Sample posterior predictive using the guide as the approximate posterior.
        Returns a dict with sampled deterministic traces and samples from the model.
        """
        # use Predictive with the guide to sample latent variables then conditionally sample observation nodes
        predictive = Predictive(self.model, guide=self.guide, num_samples=num_samples,
                                return_sites=["mu_biol", "spot_factors", "nUMI_factors", "mu_concat"])
        samples = predictive()
        return samples

    def get_param_states(self):
        """Return pyro.param_store as a dict for inspection."""
        return {k: v.detach().cpu().numpy() for k, v in pyro.get_param_store().items()}

