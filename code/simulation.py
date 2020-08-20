import numpy as np
import pandas as pd
from scipy.special import gamma

def get_n_causal(n_causal_per_study, prop_colocalizing):
    """
    find the number of causal snps to get close to specified prop_colocalizing
    """
    nCk = lambda n, k: gamma(n+1)/gamma(k+1)/gamma(n-k+1)

    def p_coloc(n):
        m = n_causal_per_study
        if m > n:
            return np.inf
        else:
            return 1 - nCk(n-m, m)/nCk(n, m)

    n_causal = np.argmin(np.abs(np.array([p_coloc(n) for n in range(100)]) - prop_colocalizing))

    if np.isclose(prop_colocalizing, 1):
        return n_causal_per_study
    elif np.isclose(n_causal, n_causal_per_study):
        return n_causal + 1
    else:
        return n_causal


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


def sim_expression_single_study(X, causal, pve, effect_distribution='normal'):
    """
    X: m x n
    return simulated_expression, true_effects, residual_variance
    """
    true_effects = np.zeros(X.shape[1])
    if effect_distribution is 'normal':
        true_effects[causal] = np.random.normal(size=np.atleast_1d(causal).size)

    if effect_distribution is 'constant':
        true_effects[causal] = 1.0

    prediction = X @ true_effects
    residual_variance = compute_sigma2(prediction, pve)
    expression = prediction + np.random.normal(size=prediction.size) * np.sqrt(residual_variance)
    return expression, true_effects, residual_variance


def sim_n_causal_per_study(X, n_study, prop_colocalizing, n_causal_per_study, pve, effect_distribution):
    """
    n_study: number of studies to simulate
    n_causal_per_tissue: number of causal snps in each tissue
    n_causal: total pool of causal snps
    pve: percent variance explained by genotype
    """
    n_variants = X.shape[1]

    if prop_colocalizing > 0:
        causal_snps = np.random.choice(n_variants, get_n_causal(n_causal_per_study, prop_colocalizing))
        results = []
        for t in range(n_study):
            causal_in_study = np.random.choice(causal_snps, n_causal_per_study, replace=False)
            results.append(sim_expression_single_study(X, causal_in_study, pve, effect_distribution))
    else:
        # special case with no colocalization: give each study its own set of causal snps
        results = []
        valid_snps = np.arange(n_variants)
        for t in range(n_study):
            causal_in_study = np.random.choice(valid_snps, n_study, replace=False)
            valid_snps = np.delete(valid_snps, causal_in_study)
            results.append(sim_expression_single_study(X, causal_in_study, pve, effect_distribution))

    expression = np.atleast_2d(np.array([x[0] for x in results]))
    true_effects = np.atleast_2d(np.array([x[1] for x in results]))
    residual_variance = np.array([x[2] for x in results])

    # trim down to the causal snps we actually used
    causal_snps = np.arange(n_variants)[np.any(true_effects != 0, 0)]

    tril = np.tril_indices(n_study, k=-1)
    true_coloc = (true_effects @ true_effects.T != 0)[tril]
    return {
        'expression': expression,
        'true_effects': true_effects,
        'true_coloc': true_coloc,
        'residual_variance': residual_variance,
        'causal_snps': causal_snps,
        'n_causal': causal_snps.size,
        'K': int(np.ceil(causal_snps.size/10) * 10)
    }
