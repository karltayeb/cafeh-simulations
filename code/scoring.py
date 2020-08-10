import numpy as np

p_coloc = lambda active, t1, t2: 1 - \
    np.exp(np.sum(np.log(1e-10 + 1 - active[t1] * active[t2])))

def score_coloc_cafeh(active, true_effects, thresh=0.9):
    """
    compute true/false positive/negative frequency
    from q(s) for all components
    """
    n_studies = active.shape[0]
    tril = np.tril_indices(n_studies, k=-1)
    true_coloc = (true_effects @ true_effects.T != 0)[tril]

    model_coloc = np.concatenate(
        [[p_coloc(active, t1, t2) for t1 in range(t2)]
         for t2 in range(n_studies)])
    model_coloc = model_coloc > thresh
    return {
        'true_positive': (true_coloc & model_coloc).sum(),
        'false_positive': (~true_coloc & model_coloc).sum(),
        'true_negative': (~true_coloc & ~model_coloc).sum(),
        'false_negative': (true_coloc & ~model_coloc).sum()
    }

def score_finemapping_cafeh(credible_sets, purity, true_effects):
    """
    compute finemapping scores
    """

    # get the causal snps
    causal_snps = np.where(np.any(true_effects != 0, 0))[0]

    # how many of these components are pure?
    pure = np.array([purity[k] > 0.5 for k in range(len(purity))])

    # how many component contain a causal snp
    cs_contains_causal = np.array([np.any(np.isin(credible_sets[k], causal_snps))
                                   for k in range(len(credible_sets))])
    try:
        all_credible_snps = np.unique(np.concatenate(
            [credible_sets[k] for k in range(len(credible_sets)) if pure[k]]))
    except ValueError:
        all_credible_snps = np.array([])

    # how many causal snps are captured in a pure component
    causal_in_credible_set = np.isin(causal_snps, all_credible_snps)

    # how many causal snps are the first snp in a component
    top_credible_snps = np.unique([credible_sets[k][0] for k in range(len(credible_sets)) if pure[k]])
    top_causal = np.isin(causal_snps, top_credible_snps)

    return {
        'n_components': pure.sum(),
        'n_components_with_causal': (cs_contains_causal & pure).sum(),
        'n_causal_in_cs': causal_in_credible_set.sum(),
        'n_top_causal': top_causal.sum(),
        'n_causal': causal_snps.size
    }
