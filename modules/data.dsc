# Module for simulating data from GTEx genotype

# Module output
# =============
# $X: genotype m x n
# $Y: expression t x n
# $LD: LD matrix (correlation between genotypes)

full_genotype: data.py + Python(genotype=load_genotype(gene_list[DSC_REPLICATE], subset).values)
  tag: "full"
  subset: None
  $genotype: genotype
  $n_sample: genotype.shape[0]
  $n_variants: genotype.shape[1]

small_genotype(full_genotype):
  tag: 'small_genotype'
  subset: 1000

normalize: data.py + Python(X = center_mean_impute(genotype).values; X / X.std(0) if normalize else X)
  normalize: True, False
  genotype: $genotype
  $X: X

genotype2ld: data.py + Python(LD = np.corrcoef(center_mean_impute(genotype).values, rowvar=False))
  genotype: $genotype
  $LD: LD

individual2summary: data.py + Python(sumstats = get_cafeh_summary_stats(Y.T, X))
  @CONF: python_modules = (numpy)
  X: $X
  Y: $Y
  $B: sumstats['beta']
  $se: sumstats['se']
  $S: sumstats['S']