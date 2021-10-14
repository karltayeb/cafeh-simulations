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
sv_gene_list = np.loadtxt(config['sim_sv_gene_path'], dtype=str)

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

def load_genotype(gene, subset=None, dense=True):
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
        if dense:
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
    # subprocess.run('rm {}*'.format(genotype_path[:-4]), shell=True)
    return genotype


def make_plink_cmd_sv(gene, save_path):
    tss = get_tss(gene)
    cmd = ' '.join(
        ['plink',
         '--bfile', config['genotype_binary_path'],
         '--chr', get_chr(gene)[3:],
         '--from-bp', str(np.maximum(tss-1e6, 0)),
         '--to-bp', str(tss+1e6),
        '--keep', config['sv_fam_path'],
         '--maf', '0.01',
         '--geno', '0.1',
         '--recode', 'A',
         '--keep-allele-order',
         '--snps-only',
         '--out', save_path])
    print(cmd)
    return cmd

def make_plink_beagle_cmd_sv(gene, save_path):
    tss = get_tss(gene)
    cmd = ' '.join(
        ['plink',
         '--bfile', config['genotype_binary_path'],
         '--chr', get_chr(gene)[3:],
         '--from-bp', str(np.maximum(tss-1e6, 0)),
         '--to-bp', str(tss+1e6),
         '--keep', config['sv_fam_path'],
         '--maf', '0.01',
         '--geno', '0.1',
         '--recode', 'beagle',
         '--keep-allele-order',
         '--snps-only',
         '--out', save_path])
    print(cmd)
    return cmd

def load_genotype_sv(gene, subset=None, dense=True):
    if not os.path.isdir('.tmp/genotype_sv'):
        os.makedirs('.tmp/genotype_sv')
    genotype_path = '.tmp/genotype_sv/{}.raw'.format(gene)
    map_path = '.tmp/genotype_sv/{}.chr-{}.map'.format(gene, get_chr(gene)[3:])
    if not os.path.isfile(genotype_path):
        print('getting genotype')
        subprocess.run(make_plink_cmd_sv(gene, genotype_path[:-4]), shell=True)
        subprocess.run(make_plink_beagle_cmd_sv(gene, genotype_path[:-4]), shell=True)

    # load SNP genotype
    genotype = pd.read_csv(genotype_path, sep=' ').set_index('IID').iloc[:, 5:]
    if subset is not None:
        if dense:
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


    # load structural variants genotype
    gene2sv = json.load(open(config['gene2sv_path'], 'r'))
    sv = pd.read_csv(config['sv_raw_path'], sep='\t')#.loc[:, gene2sv.get('ENSG00000163513.17')]
    sv = sv.set_index('IID').iloc[:, 5:]
    sv.columns = ['_'.join(c.split('_')[:-1]) for c in sv.columns]
    sv = sv.loc[:, gene2sv.get(gene)]

    idx = np.intersect1d(sv.index, genotype.index)

    # clean up
    # subprocess.run('rm {}*'.format(genotype_path[:-4]), shell=True)
    return genotype.loc[idx], sv.loc[idx]

def compute_ldscore(X):
    n_samples = X.shape[0]
    R2 = np.corrcoef(X.T) ** 2
    R2_adj = R2 - (1 - R2) / (n_samples - 2)
    ldscore = R2_adj.sum(1) - np.diag(R2_adj) + 1
    return ldscore

def compute_afreq(G):
    return (G.sum(0) / (~G.isna()).sum(0) / 2).values

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
