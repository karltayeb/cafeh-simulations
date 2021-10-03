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

    n_causal = np.argmin(np.abs(np.array(
        [p_coloc(n) for n in range(100)]) - prop_colocalizing))

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
        # if variance is 0, there were no causal variants--
        # dont care what the variance is
        sigma2_t = 1.0
    return sigma2_t


def sim_expression_single_study(X, afreq, causal, pve, effect_distribution='normal'):
    """
    X: m x n
    return simulated_expression, true_effects, residual_variance

    SPECIAL: IF PVE = 0 just set residual variance to 1
    this is useful for looking at the mixture case
    where we need effect size to be interpretable
    """
    true_effects = np.zeros(X.shape[1])
    if effect_distribution == 'normal':
        true_effects[causal] = np.random.normal(
            size=np.atleast_1d(causal).size)

    elif effect_distribution == 'normal-mixture':
        true_effects[causal] = np.random.normal(
            size=np.atleast_1d(causal).size)
        # rescale effects to be drawn from a mixture of normals
        # x ~ 1/3 N(0, 0.1) + 1/3 N(0, 1) + 1/3 N(0, 10)
        true_effects = true_effects * np.random.choice(np.sqrt([0.01, 0.1, 1.0, 10]),
            true_effects.size).reshape(true_effects.shape)

    elif effect_distribution == 'point-mixture':
        true_effects[causal] = 1.0
        # rescale effects to be drawn from a mixture of normals
        # x ~ 1/3 N(0, 0.1) + 1/3 N(0, 1) + 1/3 N(0, 10)
        true_effects = true_effects * np.random.choice(np.sqrt([0.01, 0.1, 1.0, 10]),
            true_effects.size).reshape(true_effects.shape)

    elif effect_distribution == 'constant':
        true_effects[causal] = 1.0
    else:
        assert(False)

    # sample effect size var \propto 1/p(1-p)
    true_effects = true_effects  / np.sqrt(afreq * (1 - afreq))
    prediction = X @ true_effects

    if pve == 0:
        residual_variance = 1.0
    else:
        residual_variance = compute_sigma2(prediction, pve)

    # normalize so that residual variance is one
    # std = np.sqrt(residual_variance)
    # prediction = prediction / std
    # true_effects = true_effects / std
    # residual_variance = 1.0
    
    expression = prediction + np.random.normal(
        size=prediction.size) * np.sqrt(residual_variance)
    return expression, true_effects, residual_variance


def select_causal_snps(R2, n_causal, min_r2, max_r2, active=None):
    """
    pick `n_causal` snps with pairwise R2 between `min_r2` and `max_r2`
    `active`: boolean array True if variant can be causal
    """
    causal_snps = []
    
    # initialize p
    n = R2.shape[0]
    if active is None:
        p = np.ones(n) / n
    else:
        p = active / active.sum()

    # select variants
    restarts = 0
    while(len(causal_snps) < n_causal):
        next_snp = np.random.choice(n, p=p)
        p[next_snp] = 0
        p = p * (R2[next_snp] < max_r2) * (R2[next_snp] > min_r2)
        if p.sum() == 0:
            print('restart')
            causal_snps = []
            p = np.ones(n) / n
            restarts += 1
            assert(restarts < 1e3)
            continue
        p = p / p.sum()
        causal_snps.append(next_snp)
    causal_snps = np.array(causal_snps)
    return causal_snps


def sim_block_study(X, afreq, ldscore, n_study, n_blocks, n_causal_per_block, block_p, pve, effect_distribution, min_r2, max_r2, min_ldscore, max_ldscore):
    """
    each block has a causal snps, each study assigned to a main block
    tissues within a block share the causal snp
    tissues out of block have the causal snp with probability block_p

    min_r2, max_r2: pairwise r2 values allowed for causal snps
    min_lscore, max_ldscore: ldscore range allowed for causal variants
    """
    n_variants = X.shape[1]
    n_samples = X.shape[0]
    n_causal = n_blocks * n_causal_per_block

    R2 = np.corrcoef(X.T) ** 2

    active = (ldscore > min_ldscore) & (ldscore < max_ldscore)
    causal_snps = select_causal_snps(R2, n_causal, min_r2, max_r2, active).reshape(
        n_causal_per_block, -1)

    # draw block ids and causal snps
    # ensure each block gets at least one tissue
    block_id = np.sort(np.concatenate(
        [np.arange(n_blocks), np.random.choice(n_blocks, n_study - n_blocks)]))

    # make block probability matrix
    causal_p = np.eye(n_blocks)
    causal_p[causal_p == 0] = block_p

    results = []
    for t in range(n_study):
        # sample causal snps for tissue
        causal_idx = np.random.binomial(1, causal_p[block_id[t]]) == 1
        causal_in_study = np.concatenate([cs[causal_idx] for cs in causal_snps])
        results.append(sim_expression_single_study(
            X, afreq, causal_in_study, pve, effect_distribution))

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
        'K': int(np.ceil(causal_snps.size/10) * 10),
        'ldscore': ldscore
    }
