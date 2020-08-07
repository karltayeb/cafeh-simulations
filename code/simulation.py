import numpy as np
import pandas as pd

def compute_sigma2(prediction, pve):
    """
    compute residual variance for fixed pve
    """
    var = np.var(prediction)
    sigma2_t = var/pve - var
    if sigma2_t == 0:
        # if variance is 0, there were no causal variants-- dont care what the variance is
        sigma2_t = 1.0
    return sigma2_t


def sim_expression_single_study(X, causal, pve, effect_distribution='unit_normal'):
    """
    X: m x n
    return simulated_expression, true_effects, residual_variance
    """
    true_effects = np.zeros(X.shape[1])
    if effect_distribution is 'unit_normal':
        true_effects[causal] = np.random.normal(size=np.atleast_1d(causal).size)

    prediction =  X @ true_effects
    residual_variance = compute_sigma2(prediction, pve)
    expression = prediction + np.random.normal(size=prediction.size) * np.sqrt(residual_variance)
    return expression, true_effects, residual_variance


def sim_n_causal_per_study(X, n_study, n_causal, n_causal_per_study, pve):
    """
    n_study: number of studies to simulate
    n_causal_per_tissue: number of causal snps in each tissue
    n_causal: total pool of causal snps
    pve: percent variance explained by genotype
    """
    n_variants = X.shape[1]
    causal_snps = np.random.choice(n_variants, n_causal)

    results = []
    for t in range(n_study):
        causal_in_study = np.random.choice(causal_snps, n_causal_per_study, replace=False)
        results.append(sim_expression_single_study(X, causal_in_study, pve))

    expression = np.atleast_2d(np.array([x[0] for x in results]))
    true_effects = np.atleast_2d(np.array([x[1] for x in results]))
    residual_variance = np.array([x[2] for x in results])
    return {
        'expression': expression,
        'true_effects': true_effects,
        'residual_variance': residual_variance,
        'causal_snps': causal_snps
    }
