from types import SimpleNamespace
import numpy as np
import pandas as pd

p_coloc = lambda active, t1, t2: 1 - \
    np.exp(np.sum(np.log(1e-10 + 1 - active[t1] * active[t2])))

def score_coloc_cafeh(active, purity, true_coloc, thresh=0.9, filter_purity=False):
    """
    compute true/false positive/negative frequency
    from q(s) for all components
    """
    if filter_purity:
        pure = np.array([purity[k] > 0.5 for k in purity])
        active = active * pure[None]

    p_coloc = lambda t1, t2: 1 - \
        np.exp(np.sum(np.log(1e-10 + 1 - active[t1] * active[t2])))

    n_studies = active.shape[0]
    tril = np.tril_indices(n_studies, k=-1)
    model_coloc = np.concatenate(
        [[p_coloc(t1, t2) for t1 in range(t2)] for t2 in range(n_study)])
    return {
        'p_coloc': model_coloc,
        'true_positive': (true_coloc & (model_coloc > thresh)).sum(),
        'false_positive': (~true_coloc & (model_coloc > thresh)).sum(),
        'true_negative': (~true_coloc & ~(model_coloc > thresh)).sum(),
        'false_negative': (true_coloc & ~(model_coloc > thresh)).sum()
    }


def score_coloc_coloc(coloc_out, true_coloc, thresh=0.9):
    n_study = np.unique(
        np.array([np.array(k) for k in coloc_out])).size
    tril = np.tril_indices(n_study, k=-1)
    p_coloc = np.concatenate(
        [[coloc_out[(t1, t2)].pph4 for t1 in range(t2)] for t2 in range(n_study)])
    return {
        'p_coloc': p_coloc,
        'true_positive': (true_coloc & (p_coloc > thresh)).sum(),
        'false_positive': (~true_coloc & (p_coloc > thresh)).sum(),
        'true_negative': (~true_coloc & ~(p_coloc > thresh)).sum(),
        'false_negative': (true_coloc & ~(p_coloc > thresh)).sum()
    }


def score_coloc_ecaviar(ecaviar_out, true_coloc, thresh=0.9):
    n_study = np.unique(
        np.array([np.array(k) for k in ecaviar_out])).size
    tril = np.tril_indices(n_study, k=-1)
    p_coloc = np.concatenate(
        [[ecaviar_out[(t1, t2).max()] for t1 in range(t2)] for t2 in range(n_study)])
    return {
        'p_coloc': p_coloc,
        'true_positive': (true_coloc & (p_coloc > thresh)).sum(),
        'false_positive': (~true_coloc & (p_coloc > thresh)).sum(),
        'true_negative': (~true_coloc & ~(p_coloc > thresh)).sum(),
        'false_negative': (true_coloc & ~(p_coloc > thresh)).sum()
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

def score_finemapping_caviar(caviar_out, true_effects):
    pip = pd.concat([c.posterior.pip for c in caviar_out], axis=1).T.values
    credible_sets = {i: caviar_out[i].credible_set for i in range(len(caviar_out))}

    n_variants = true_effects.shape[1]
    causal_snps = {i: np.arange(n_variants)[te != 0] for i, te in enumerate(true_effects)}
    n_causal_in_cs = np.array([np.isin(credible_sets[i], causal_snps[i]) for i in causal_snps])
    n_causal = np.array([causal_snps[i].size for i in causal_snps])
    return SimpleNamespace(
        study_pip=pip, credible_sets=credible_sets, n_causal_in_cs=n_causal, n_causal=n_causal)
