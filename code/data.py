from types import SimpleNamespace
import subprocess
import copy
import json
import yaml
import os
import numpy as np
import pandas as pd


config = yaml.load(open('config.yml', 'r'))
gene_list = np.loadtxt(config['sim_gene_path'], dtype=str)

def get_chr(gene):
    gene2chr = json.load(open(config['gene2chr_path'], 'r'))
    return gene2chr.get(gene, None)


def get_tss(gene):
    gene2chr = json.load(open(config['gene2tss_path'], 'r'))
    return gene2chr.get(gene, None)


def make_plink_cmd(gene, save_path):
    tss = get_tss(gene)
    cmd = ' '.join(
        ['plink',
         '--bfile', config['genotype_binary_path'],
         '--chr', get_chr(gene)[3:],
         '--from-bp', str(np.maximum(tss-1e6, 0)),
         '--to-bp', str(tss+1e6),
         '--maf', '0.01',
         '--geno', '0.1',
         '--recode', 'A',
         '--keep-allele-order',
         '--snps-only',
         '--out', save_path])
    print(cmd)
    return cmd

def make_plink_beagle_cmd(gene, save_path):
    tss = get_tss(gene)
    cmd = ' '.join(
        ['plink',
         '--bfile', config['genotype_binary_path'],
         '--chr', get_chr(gene)[3:],
         '--from-bp', str(np.maximum(tss-1e6, 0)),
         '--to-bp', str(tss+1e6),
         '--maf', '0.01',
         '--geno', '0.1',
         '--recode', 'beagle',
         '--keep-allele-order',
         '--snps-only',
         '--out', save_path])
    print(cmd)
    return cmd

def load_genotype(gene, subset=None, random=False):
    if not os.path.isdir('.tmp/genotype'):
        os.makedirs('.tmp/genotype')
    genotype_path = '.tmp/genotype/{}.raw'.format(gene)
    map_path = '.tmp/genotype/{}.chr-{}.map'.format(gene, get_chr(gene)[3:])
    if not os.path.isfile(genotype_path):
        print('getting genotype')
        subprocess.run(make_plink_cmd(gene, genotype_path[:-4]), shell=True)
        subprocess.run(make_plink_beagle_cmd(gene, genotype_path[:-4]), shell=True)

    genotype = pd.read_csv(genotype_path, sep=' ').set_index('IID').iloc[:, 5:]
    if subset is not None:
        if not random:
            print('use variants near tss')
            tss = get_tss(gene)
            pos = pd.read_csv(map_path, sep='\t').iloc[:, 1].values
            idx = np.argsort(np.abs(pos - tss))[:subset]
        else:
            print('use random subset of variants near tss')
            idx = np.sort(np.random.choice(
                genotype.shape[1], subset
            ))
        genotype = genotype.iloc[:, idx]

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
