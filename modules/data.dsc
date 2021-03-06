# Module for simulating data from GTEx genotype

# Module output
# =============
# $X: genotype m x n
# $Y: expression t x n
# $LD: LD matrix (correlation between genotypes)

full_genotype: data.py + Python(gene=gene_list[DSC_REPLICATE % 200]; G=load_genotype(gene, subset, dense); X=center_mean_impute(G).values; afreq=compute_afreq(G); ldscore=compute_ldscore(X))
  tag: "full"
  subset: None
  dense: True
  $X: X
  $afreq: afreq
  $ldscore: ldscore
  $n_sample: X.shape[0]
  $n_variants: X.shape[1]
  $gene: gene

small_genotype(full_genotype):
  tag: '1k_genotype'
  subset: 1000
  dense: True

small_genotype_random(full_genotype):
  tag: '1k_genotype_random'
  subset: 1000
  dense: False

medium_genotype(full_genotype):
  tag: '2k_genotype'
  subset: 2000

large_genotype(full_genotype):
  tag: '5k_genotype'
  subset: 5000

genotype2ld: Python(LD = numpy.corrcoef(X, rowvar=False))
  @CONF: python_modules = (numpy)
  X: $X
  $LD: LD

individual2summary: data.py + Python(sumstats = get_cafeh_summary_stats(Y.T, X))
  X: $X
  Y: $Y
  $B: sumstats['beta']
  $se: sumstats['se']
  $S: sumstats['S']