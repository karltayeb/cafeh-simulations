capture: Python(pickle.dump(X, open('./tmp/X', 'wb')); pickle.dump(Y, open('./tmp/Y', 'wb')); a=X)
  @CONF: python_modules = (numpy)
  X: $X
  Y: $Y
  $a: a

fit_cafeh_genotype: fitting.py + Python(model = fit_cafeh_genotype(X=X.T, Y=Y, K=K); params = get_param_dict(model))
  X: $X
  Y: $Y
  K: $K
  $expected_effects: model.expected_effects
  $pip: model.get_pip()
  $study_pip: model.get_study_pip().values
  $pi: model.pi
  $active: model.active
  $credible_sets: model.credible_sets
  $purity: model.purity
  $params: params


fit_cafeh_summary: fitting.py + Python(model = fit_cafeh_summary(LD, B, S, K=K); params = get_param_dict(model))
  LD: $LD
  B: $B
  S: $S
  K: $K
  $expected_effects: model.expected_effects
  $pip: model.get_pip()
  $study_pip: model.get_study_pip().values
  $pi: model.pi
  $active: model.active
  $credible_sets: model.credible_sets
  $purity: model.purity
  $params: params


fit_caviar: fitting.py + Python(caviar_out = run_caviar(B, se, LD))
  LD: $LD
  B: $B
  se: $se
  CAVIAR_PATH: '/Users/karltayeb/Research/software/caviar/CAVIAR-C++/CAVIAR'
  $caviar_out: caviar_out


fit_ecaviar: fitting.py + Python(ecaviar_out = ecaviar_from_caviar(caviar_out))
  caviar_out: $caviar_out
  $ecaviar_out: ecaviar_out


fit_coloc: fitting.py + Python(coloc_out = run_coloc(B, se))
  B: $B
  se: $se
  $coloc_out: coloc_out
