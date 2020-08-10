# Module for simulating data from GTEx genotype

# Module output
# =============
# $X: genotype m x n
# $Y: expression t x n
# $LD: LD matrix (correlation between genotypes)

full_genotype: data.py + Python(gene_list = np.loadtxt(gene_list_path, dtype=str); X=center_mean_impute(load_genotype(np.random.choice(gene_list), subset)).values)
  tag: "full"
  gene_list_path: 'data/sim_genes.txt'
  subset: None
  $X: X
  $n_sample: X.shape[0]
  $n_variants: X.shape[1]

small_genotype(full_genotype):
  tag: 'small_genotype'
  subset: 1000

genotype2ld: Python(LD = numpy.corrcoef(X, rowvar=False))
  @CONF: python_modules = (numpy)
  X: $X
  $LD: LD

individual2summary: data.py + Python(sumstats = get_cafeh_summary_stats(Y.T, X))
  @CONF: python_modules = (numpy)
  X: $X
  Y: $Y
  $B: sumstats['beta']
  $se: sumstats['se']
  $S: sumstats['S']

random_data: Python(X = numpy.random.random((100, 100)))
  @CONF: python_modules = (numpy)
  $X: X2