# Module for simulating data from GTEx genotype

# Module output
# =============
# $X: genotype m x n
# $Y: expression t x n
# $LD: LD matrix (correlation between genotypes)

full_genotype: data.py + Python(X=center_mean_impute(load_genotype(gene, subset)).values)
  tag: "full"
  gene: 'ENSG00000284733.1'
  subset: None
  $X: X
  $n_sample: X.shape[0]
  $n_variants: X.shape[1]

small_genotype(full_genotype):
  subset: 1000

genotype2ld: Python(LD = numpy.corrcoef(X, rowvar=False))
  @CONF: python_modules = (numpy)
  X: $X
  $LD: LD

individual2summary: data.py + Python(sumstats = get_cafeh_summary_stats(Y.T, X))
  @CONF: python_modules = (numpy)
  X: $X
  Y: $Y
  $B: sumstats['B']
  $S: sumstats['S']

random_data: Python(X = numpy.random.random((100, 100)))
  @CONF: python_modules = (numpy)
  $X: X2