from cafeh.cafeh_ss import CAFEH
from cafeh.cafeh_ss_simple import CAFEHSimple
from cafeh.independent_model_ss import CAFEHG
from cafeh.fitting import forward_fit_procedure, weight_ard_active_fit_procedure, weight_active_fit_procedure
from coloc import coloc
import numpy as np
import pandas as pd
import subprocess
import os
from itertools import combinations
from types import SimpleNamespace
import yaml

try:
    config = yaml.load(open('config.yml', 'r'))
except Exception as e:
    print(e)

def _fit(model, update_ard, update_active, update_variance):
    fit_args = {
        'update_weights': True,
        'update_pi': True,
        'update_variance': update_variance,
        'ARD_weights': False,
        'update_active': False,
        'max_iter': 50
    }
    
    model.fit(**fit_args)

    if update_ard:
        fit_args['ARD_weights'] = True
        model.fit(**fit_args)

    if update_active:
        fit_args['update_active'] = True
        model.fit(**fit_args)


def _fit_soft_init(model, update_ard, update_active, update_variance):
    fit_args = {
        'update_weights': True,
        'update_pi': True,
        'update_variance': update_variance,
        'ARD_weights': False,
        'update_active': False,
        'max_iter': 50
    }
    
    w_prior_variance = model.weight_precision_b

    wpv = 1e-6
    while wpv < np.max(w_prior_variance):
        model.weight_precision_b = np.ones_like(w_prior_variance) * wpv
        model.fit(**fit_args)
        wpv = wpv * 10

    model.weight_precision_b = w_prior_variance
    model.fit(**fit_args)

    if update_ard:
        fit_args['ARD_weights'] = True
        model.fit(**fit_args)

    if update_active:
        fit_args['update_active'] = True
        model.fit(**fit_args)

def get_param_dict(model, compress=True):
    param_dict = {}
    if compress:
        model._compress_model()
    for key in model.__dict__:
        if key not in ['LD', 'B', 'X', 'Y', 'precompute', 'records']:
            param_dict[key] = model.__dict__[key]
    if compress:
        model._decompress_model()
    return param_dict

def fit_cafeh_genotype(X, Y, K, p0k, w_prior_variance, standardize, update_ard, update_active,
    update_variance, **kwargs):
    if standardize:
        X = (X - X.mean(1)[:, None]) / X.std(1)[:, None]

    model = CAFEHG(X=X, Y=Y, K=K)
    model.prior_activity = np.ones(K) * p0k
    model.weight_precision_b = np.ones_like(model.weight_precision_b) * w_prior_variance

    _fit(model, update_ard, update_active, update_variance)
    model.clear_precompute()
    return model

def fit_cafeh_genotype_soft_init(X, Y, K, p0k, w_prior_variance, standardize, update_ard, update_active,
    update_variance, **kwargs):
    if standardize:
        X = (X - X.mean(1)[:, None]) / X.std(1)[:, None]

    model = CAFEHG(X=X, Y=Y, K=K)
    model.prior_activity = np.ones(K) * p0k
    model.weight_precision_b = np.ones_like(model.weight_precision_b) * w_prior_variance

    _fit_soft_init(model, update_ard, update_active, update_variance)
    model.clear_precompute()
    return model


def fit_cafeh_genotype_pairwise(X, Y, K, p0k, w_prior_variance, standardize, update_ard, update_active,
    update_variance, **kwargs):
    n_study = Y.shape[0]

    if standardize:
        X = (X.T / X.T.std(0)).T

    _p_coloc = lambda t1, t2, m: 1 - \
        np.exp(np.sum(np.log(1e-10 + 1 - m.active[t1] * m.active[t2])))
    p_coloc = {}
    for i, j in combinations(np.arange(n_study), 2):
        model = CAFEHG(X=X, Y=Y[[i, j]], K=K)
        model.prior_activity = np.ones(K) * p0k
        model.weight_precision_b = np.ones_like(model.weight_precision_b) *w_prior_variance
        _fit(model, update_ard, update_active, update_variance)
        p_coloc[(i, j)] = _p_coloc(0, 1, model)

    p_coloc = np.concatenate(
            [[p_coloc[(t1, t2)] for t1 in range(t2)] for t2 in range(n_study)])
    return {'p_coloc': p_coloc}


def fit_susie_genotype(X, Y, K, p0k, w_prior_variance, standardize, update_ard, update_active,
    update_variance, **kwargs):
    expected_effects = []
    study_pip = []
    credible_sets = []
    purity = []
    params = []

    if standardize:
        X = (X.T / X.T.std(0)).T

    for y in Y:
        model = CAFEHG(X=X, Y=y[None], K=K)
        model.prior_activity = np.ones(K) * p0k
        model.weight_precision_b = np.ones_like(model.weight_precision_b) *w_prior_variance
        _fit(model, update_ard, update_active, update_variance)
        model.clear_precompute()

        expected_effects.append(model.expected_effects)
        study_pip.append(model.get_study_pip().values.flatten())
        credible_sets.append(model.credible_sets)
        purity.append(model.purity)
        params.append(get_param_dict(model))

    expected_effects = np.array(expected_effects)
    study_pip = np.array(study_pip)
    return SimpleNamespace(
        expected_effects=expected_effects, study_pip=study_pip, credible_sets=credible_sets,
        purity=purity, params=params)


def fit_cafeh_summary(LD, B, se, S, K, p0k, w_prior_variance, standardize, update_ard, update_active, **kwargs):
    if standardize:
        B = B / S
        S = np.ones_like(S)

    model = CAFEH(LD=LD, B=B, S=S, K=K)
    model.prior_activity = np.ones(K) * p0k
    model.weight_precision_b = np.ones_like(model.weight_precision_b) *w_prior_variance
    _fit(model, update_ard, update_active, False)
    model.clear_precompute()
    return model


def fit_susie_summary(LD, B, se, S, K, p0k, w_prior_variance, standardize, update_ard, update_active, **kwargs):
    expected_effects = []
    study_pip = []
    credible_sets = []
    purity = []
    params = []

    if standardize:
        B = B / S
        S = np.ones_like(S)

    for i in range(B.shape[0]):
        model = CAFEH(LD=LD, B=B[[i]], S=S[[i]], K=K)
        model.prior_activity = np.ones(K) * p0k
        model.weight_precision_b = np.ones_like(model.weight_precision_b) * w_prior_variance
        _fit(model, update_ard, update_active, False)
        model.clear_precompute()

        expected_effects.append(model.expected_effects)
        study_pip.append(model.get_study_pip().values.flatten())
        credible_sets.append(model.credible_sets)
        purity.append(model.purity)
        params.append(get_param_dict(model))

    expected_effects = np.array(expected_effects)
    study_pip = np.array(study_pip)
    return SimpleNamespace(
        expected_effects=expected_effects, study_pip=study_pip, credible_sets=credible_sets,
        purity=purity, params=params)


def fit_cafeh_summary_simple(LD, B, se, S, K, p0k, w_prior_variance, standardize, fit, **kwargs):
    if standardize:
        B = B / S
        S = np.ones_like(S)

    model = CAFEHSimple(LD=LD, B=B, S=S, K=K)
    model.prior_activity = np.ones(K) * p0k
    model.weight_precision_b = np.ones_like(model.weight_precision_b) * w_prior_variance
    _fit(model, update_ard, update_active, update_variance)
    model.clear_precompute()
    return model

############
#  CAVIAR  #
############
def make_caviar_command(CAVIAR_PATH, ld_path, z_path, c=2, o='.tmp/caviar'):
    cmd = ' '.join(
        [CAVIAR_PATH,
         '-l', ld_path,
         '-z', z_path,
         '-c', str(c),
         '-o', str(o)])
    return cmd


def run_caviar(B, se, LD):
    """
    run caviar, return list of caviar POST dataframes
    use this downstream for eCAVIAR and finemapping evaluation
    """
    prefix_path = '.tmp/caviar/'
    if not os.path.isdir(prefix_path):
        os.makedirs(prefix_path)

    # seed hack so that tmp files dont get overwritten
    import time
    np.random.seed(int((time.time() * 1e6) % (2**32 -1)))
    prefix = prefix_path + ''.join(np.random.choice(10, 20).astype(str))
    np.random.seed(DSC_SEED)



    z = B / se
    # save caviar summary stats: zscore, LD
    for i, z_i in enumerate(z):
        pd.Series(z_i).to_csv('{}z{}'.format(prefix, i), sep='\t', header=None)

    LD_path = '{}LD'.format(prefix)
    np.savetxt(LD_path, LD, fmt='%.2f')

    #run caviar
    post = []
    results = []
    for i in range(z.shape[0]):
        cmd = make_caviar_command(
            config['caviar_path'], '{}LD'.format(prefix),
            '{}z{}'.format(prefix, i), c=2,
            o='{}_{}'.format(prefix, i))
        print(cmd)
        subprocess.run(cmd, shell=True)

        post = pd.read_csv('{}_{}_post'.format(prefix, i), sep='\t')
        post = post.rename(columns={
            'SNP_ID': 'snp',
            'Prob_in_pCausalSet': 'snp_in_cset_prob',
            'Causal_Post._Prob.': 'pip'
        })
        post = post.set_index('snp')
        try:
            cs = pd.read_csv('{}_{}_set'.format(
                prefix, i), header=None).values.flatten()
        except pd.errors.EmptyDataError:
            cs = np.array([])
        results.append(SimpleNamespace(posterior=post, credible_set=cs))

    # clean up temporary files
    subprocess.run('rm -r {}*'.format(prefix), shell=True)
    return results


def ecaviar_from_caviar(results):
    """
    results a list of caviar results of length l
    return pairwise CLPP for all (l choose 2) pairs of studies
    """
    l = len(results)
    out = {}
    for i, j in combinations(range(l), 2):
        CLPP = results[i].posterior.pip * results[j].posterior.pip
        out[(i, j)] = CLPP.values
    return out

############
#  Coloc  #
###########

def compute_lnabf(beta, se, W=0.15):
    """

    W, prior variance of effect size, set to coloc default for continuous traits
    """
    r = W / (W + se**2)
    z = beta / se
    return np.log(np.sqrt(1-r) * np.exp(z**2 / 2 * r))


def run_coloc(beta, se):
    """
    run coloc pairwise on rows of summary stats
    """
    lnabf = compute_lnabf(beta, se)
    out = {}
    for i, j in combinations(range(lnabf.shape[0]), 2):
        out[(i, j)] = SimpleNamespace(
            **{'pph{}'.format(i): x for i, x in enumerate(
                coloc(lnabf[i], lnabf[j]))})
    return out


def load_model_from_history(X, Y, params):
    model = CAFEHG(X.T, Y, K=10)
    model.__dict__.update(params)
    model._decompress_model()
    return model


############
#  FINEMAP  #
############

def run_finemap(B, se, afreq, LD):
    FINEMAP_PATH = config['finemap_path']
    prefix_path = '.tmp/finemap/'
    if not os.path.isdir(prefix_path):
        os.makedirs(prefix_path)

    # seed hack so that tmp files dont get overwritten
    import time
    np.random.seed(int((time.time() * 1e6) % (2**32 -1)))
    prefix = prefix_path + ''.join(np.random.choice(10, 20).astype(str))

    t, n = B.shape

    pd.DataFrame(LD).to_csv(prefix + 'ld', sep=' ', index=None, header=None)
    for study in range(t):
        d = {
            'rsid': ['rs{}'.format(i + 1) for i in range(n)],
            'chromosome': np.ones(n),
            'position': np.arange(n) + 1,
            'allele1': ['T' for _ in range(n)],
            'allele2': ['A' for _ in range(n)],
            'maf': np.minimum(afreq, np.abs(afreq-1)),
            'beta': B[study],
            'se': se[study]
        }

        d = pd.DataFrame(d)
        d.chromosome = d.chromosome.astype(int)
        d.to_csv(prefix + 'study{}.z'.format(study), sep=' ', index=None)

    master = {
        'z': [prefix + 'study{}.z'.format(i) for i in range(t)],
        'ld': [prefix + 'ld' for _ in range(t)],
        'snp': [prefix + 'study{}.snp'.format(i) for i in range(t)],
        'config': [prefix + 'study{}.config'.format(i) for i in range(t)],
        'cred': [prefix + 'study{}.cred'.format(i) for i in range(t)],
        'log': [prefix + 'study{}.log'.format(i) for i in range(t)],
        'n_samples': [n for _ in range(t)],
    }
    pd.DataFrame(master).to_csv(prefix + 'master', sep=';', index=None)

    cmd = './{} --in-files {}master --sss'.format(FINEMAP_PATH, prefix)
    subprocess.run(cmd, shell=True)

    # get credible sets
    results = []
    for i in range(B.shape[0]):
        p = '{}study{}.snp'.format(prefix, i)
        df = pd.read_csv('../finemap_test/study0.snp', index_col=0, sep=' ')
        df.loc[:, 'snp'] = df.rsid.apply(lambda x: int(x[2:]) - 1)
        df = df.set_index('snp')
        df = df.sort_index()

        results.append(SimpleNamespace(posterior=df))

    # get pips
    subprocess.run('rm -r {}*'.format(prefix), shell=True)
    return cmd
