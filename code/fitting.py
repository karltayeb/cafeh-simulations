from cafeh.cafeh_ss import CAFEH
from cafeh.cafeh_ss_simple import CAFEHSimple
from cafeh.independent_model_ss import CAFEHG
from cafeh.fitting import forward_fit_procedure
from coloc import coloc
import numpy as np
import pandas as pd
import subprocess
import os
from itertools import combinations
from types import SimpleNamespace
import yaml

config = yaml.load(open('config.yml', 'r'))

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

def fit_cafeh_genotype(X, Y, K, p0k):
    cafehg = CAFEHG(X=X, Y=Y, K=K)
    cafehg.prior_activity = np.ones(K) * p0k
    print(cafehg.X.shape)
    print(cafehg.Y.shape)
    print(cafehg.dims)

    forward_fit_procedure(cafehg)
    return cafehg


def fit_cafeh_summary_simple(LD, B, S, K, p0k):
    model = CAFEHSimple(LD=LD, B=B, S=S, K=K)
    model.prior_activity = np.ones(K) * p0k
    forward_fit_procedure(model)
    return model

def fit_cafeh_summary(LD, B, S, K, p0k):
    model = CAFEH(LD=LD, B=B, S=S, K=K)
    model.prior_activity = np.ones(K) * p0k
    forward_fit_procedure(model)
    return model

def fit_susie_genotype(X, Y, K, p0k):
    expected_effects = []
    study_pip = []
    credible_sets = []
    purity = []
    params = []

    for y in Y:
        model = CAFEHG(X=X, Y=y[None], K=K)
        model.prior_activity = np.ones(K) * p0k
        forward_fit_procedure(model)

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

def fit_susie_summary(LD, B, S, K, p0k):
    expected_effects = []
    study_pip = []
    credible_sets = []
    purity = []
    params = []

    for i in range(B.shape[0]):
        model = CAFEH(LD=LD, B=B[[i]], S=S[[i]], K=K)
        model.prior_activity = np.ones(K) * p0k
        forward_fit_procedure(model)
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

    z = B/se
    # save caviar summary stats: zscore, LD
    for i, z_i in enumerate(z):
        pd.Series(z_i).to_csv('{}z{}'.format(prefix, i), sep='\t', header=None)

    LD_path = '{}LD'.format(prefix)
    np.savetxt(LD_path, LD, fmt='%.2f')

    #run caviar
    post = []
    for i in range(z.shape[0]):
        cmd = make_caviar_command(
            config['caviar_path'], '{}LD'.format(prefix),
            '{}z{}'.format(prefix, i), c=2,
            o='{}_{}'.format(prefix, i))
        print(cmd)
        subprocess.run(cmd, shell=True)

    results = []
    for i in range(z.shape[0]):
        post = pd.read_csv('{}_{}_post'.format(prefix, i), sep='\t')
        post = post.rename(columns={
            'SNP_ID': 'snp', 'Prob_in_pCausalSet': 'snp_in_cset_prob', 'Causal_Post._Prob.': 'pip'})
        post = post.set_index('snp')
        try:
            cs = pd.read_csv('{}_{}_set'.format(prefix, i), header=None).values.flatten()
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
        CLPP = results[0].posterior.pip * results[1].posterior.pip
        out[(i, j)] = CLPP.values.max()
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


