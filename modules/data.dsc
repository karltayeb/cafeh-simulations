# Module for simulating data from GTEx genotype

# Module output
# =============
# $X: genotype m x n
# $Y: expression t x n
# $LD: LD matrix (correlation between genotypes)

full_genotype: data.py + Python(X=center_mean_impute(load_genotype(gene_list[DSC_REPLICATE], subset, dense)).values)
  tag: "full"
  subset: None
  dense: True
  $X: X
  $n_sample: X.shape[0]
  $n_variants: X.shape[1]

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

genotype2ld: Python(LD = numpy.corrcoef(X, rowvar=False))
  @CONF: python_modules = (numpy)
  X: $X
  $LD: LD

ld2ldscore: Python(LD_corr = LD - (1-LD)/(X.shape[0] -2); LD_score = LD_corr.sum(1) - np.diag(LD_corr) + 1)
  X: $X
  LD: $LD
  $LD_score: LD_score

individual2summary: data.py + Python(sumstats = get_cafeh_summary_stats(Y.T, X))
  X: $X
  Y: $Y
  $B: sumstats['beta']
  $se: sumstats['se']
  $S: sumstats['S']