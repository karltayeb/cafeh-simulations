from types import SimpleNamespace
import subprocess
import copy
import json
import os
import numpy as np
import pandas as pd

def get_chr(gene):
    gene2chr = json.load(open('data/gene2chr', 'r'))
    return gene2chr.get(gene, None)


def get_tss(gene):
    gene2chr = json.load(open('data/gene2tss', 'r'))
    return gene2chr.get(gene, None)


def make_plink_cmd(gene, save_path):
    tss = get_tss(gene)
    cmd = ' '.join(
        ['plink',
         '--bfile', 'data/1k_genomes/1kg.{}'.format(get_chr(gene)),
         '--chr', get_chr(gene)[3:],
         '--from-bp', str(np.maximum(tss-1e6, 0)),
         '--to-bp', str(tss+1e6),
         '--maf', '0.01',
         '--geno', '0.1',
         '--recode', 'A',
         '--keep-allele-order',
         '--snps-only',
         '--out', save_path])
    return cmd


def load_genotype(gene, subset=None):
    if not os.path.isdir('.tmp/1k_genomes'):
        os.makedirs('.tmp/1k_genomes')
    genotype_path = '.tmp/1k_genomes/{}.raw'.format(gene)
    if not os.path.isfile(genotype_path):
        print('getting genotype')
        subprocess.run(make_plink_cmd(gene, genotype_path[:-4]), shell=True)
    genotype = pd.read_csv(genotype_path, sep=' ').set_index('IID').iloc[:, 5:]
    gentoype = genotype.rename(columns={x: x.split('_')[0] for x in genotype.columns})

    if subset is not None:
        snp_subset = np.random.choice(genotype.shape[1], subset, replace=False)
        genotype = genotype.iloc[:, snp_subset]

    # clean up
    subprocess.run('rm {}*'.format(genotype_path[:-4]), shell=True)
    return genotype


def load_ld(gene, out):
    genotype = load_genotype(gene)
    ld = np.corrcoef(np.nan_to_num(
        genotype.values - np.nanmean(genotype.values, 0)).T)
    return pd.DataFrame(ld, index=genotype.columns, columns=genotype.columns)


def center_mean_impute(genotye):
    imputed = np.nan_to_num(genotye.values - np.nanmean(genotye.values, 0))
    return pd.DataFrame(imputed, genotye.index, genotye.columns)

def linregress(y, X):
    """
    y: m x t expression
    X: m x n genotype
    compute t x n pairwise linear regressions
    reports slopes/standard errors
    """
    diag = np.einsum('ij,ij->i', X.T, X.T)
    betas = y.T @ X / diag
    if np.ndim(y) == 1:
        var = np.var(y[:, None] - betas * X, 0) / diag
    else:
        var = np.array([np.var(y[:, t][:, None] - betas[t] * X, 0) / diag for t in range(y.shape[1])])
    return betas, np.sqrt(var)


def get_cafeh_summary_stats(y, X):
    """
    shape data for CAFEHS
    """
    beta, se = linregress(y, X)
    return {
        'beta': beta,
        'S': np.sqrt((beta**2 / X.shape[0]) + se**2),
        'se': se
    }
