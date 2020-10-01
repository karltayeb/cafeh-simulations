fit_cafeh_genotype: fitting.py + Python(model = fit_cafeh_genotype(X.T, Y, K, p0k, standardize, fit, max_iter=max_iter, update_variance=update_variance); params = get_param_dict(model))
  model: 'cafeh_genotype'
  X: $X
  Y: $Y
  K: $K
  p0k: 0.01
  standardize: True, False
  fit: "weight_ard_active", "weight_active"
  max_iter: 50
  update_variance: False
  $expected_effects: model.expected_effects
  $pip: model.get_pip()
  $study_pip: model.get_study_pip().values
  $pi: model.pi
  $active: model.active
  $credible_sets: model.credible_sets
  $purity: model.purity
  $params: params

fit_susie_genotype: fitting.py + Python(results = fit_susie_genotype(X.T, Y, K, p0k, standardize, fit, max_iter=max_iter, update_variance=update_variance))
  model: 'susie_genotype'
  X: $X
  Y: $Y
  K: 5
  p0k: 1.0
  standardize: True, False
  fit: "weight_ard_active", "weight_active"
  max_iter: 50
  update_variance: False
  $expected_effects: results.expected_effects
  $study_pip: results.study_pip
  $credible_sets: results.credible_sets
  $purity: results.purity
  $params: results.params

fit_susie_genotype_ss(fit_susie_genotype):
  model: 'susie_genotype_ss'
  p0k: 0.01

fit_cafeh_summary: fitting.py + Python(model = fit_cafeh_summary(LD, B, se, S, K, p0k, standardize, fit, max_iter=max_iter); params = get_param_dict(model))
  model: 'cafeh_summary'
  LD: $LD
  B: $B
  se: $se
  S: $S
  K: $K
  p0k: 0.01
  standardize: True, False
  fit: "weight_ard_active", "weight_active"
  max_iter: 50
  update_variance: False
  $expected_effects: model.expected_effects
  $pip: model.get_pip()
  $study_pip: model.get_study_pip().values
  $pi: model.pi
  $active: model.active
  $credible_sets: model.credible_sets
  $purity: model.purity
  $params: params

fit_cafeh_summary_simple: fitting.py + Python(model = fit_cafeh_summary_simple(LD, B, se, S, K, p0k, standardize, fit, max_iter=max_iter); params = get_param_dict(model, compress=False))
  model: 'cafeh_summary_simple'
  LD: $LD
  B: $B
  se: $se
  S: $S
  K: $K
  p0k: 0.01
  standardize: True, False
  fit: "weight_ard_active", "weight_active"
  max_iter: 50
  update_variance: False
  $expected_effects: model.expected_effects
  $pip: model.get_pip()
  $study_pip: model.get_study_pip().values
  $pi: model.pi
  $active: model.active
  $credible_sets: model.credible_sets
  $purity: model.purity
  $params: params

fit_susie_summary: fitting.py + Python(results = fit_susie_summary(LD, B, se, S, K, p0k, standardize, fit, max_iter=max_iter))
  model: 'susie_summary'
  LD: $LD
  B: $B
  se: $se
  S: $S
  K: 5
  p0k: 1.0
  standardize: True, False
  update_variance: False
  fit: "weight_ard_active", "weight_active"
  max_iter: 50
  $expected_effects: results.expected_effects
  $study_pip: results.study_pip
  $credible_sets: results.credible_sets
  $purity: results.purity
  $params: results.params

fit_susie_summary_ss(fit_susie_summary):
  model: 'susie_summary_ss'
  p0k: 0.01

fit_caviar: fitting.py + Python(caviar_out = run_caviar(B, se, LD))
  model: 'caviar'
  LD: $LD
  B: $B
  se: $se
  $caviar_out: caviar_out


fit_ecaviar: fitting.py + Python(ecaviar_out = ecaviar_from_caviar(caviar_out))
  model: 'ecaviar'
  caviar_out: $caviar_out
  $ecaviar_out: ecaviar_out


fit_coloc: fitting.py + Python(coloc_out = run_coloc(B, se))
  model: 'coloc'
  B: $B
  se: $se
  $coloc_out: coloc_out
