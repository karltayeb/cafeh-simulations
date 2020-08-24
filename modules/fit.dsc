fit_cafeh_genotype: fitting.py + Python(model = fit_cafeh_genotype(X_norm.T, Y, K, p0k, fit, update_variance=update_variance); params = get_param_dict(model))
  model: 'cafeh_genotype'
  X_norm: $X_norm
  Y: $Y
  K: $K
  p0k: 0.01, 0.1, 0.5
  fit: "weight_ard_active"
  update_variance: True, False
  $expected_effects: model.expected_effects
  $pip: model.get_pip()
  $study_pip: model.get_study_pip().values
  $pi: model.pi
  $active: model.active
  $credible_sets: model.credible_sets
  $purity: model.purity
  $params: params

fit_susie_genotype: fitting.py + Python(results = fit_susie_genotype(X_norm.T, Y, K, p0k, fit, update_variance=update_variance))
  model: 'susie_genotype'
  X_norm: $X_norm
  Y: $Y
  K: 5
  p0k: 1.0
  fit: "weight_ard_active"
  update_variance: True, False
  $expected_effects: results.expected_effects
  $study_pip: results.study_pip
  $credible_sets: results.credible_sets
  $purity: results.purity
  $params: results.params

fit_susie_genotype_ss(fit_susie_genotype):
  model: 'susie_genotype_ss'
  p0k: 0.01

fit_cafeh_summary: fitting.py + Python(model = fit_cafeh_summary(LD, B, S, K, p0k, fit); params = get_param_dict(model))
  model: 'cafeh_summary'
  LD: $LD
  B: $B
  S: $S
  K: $K
  p0k: 0.01, 0.1, 0.5
  fit: "weight_ard_active"
  update_variance: False
  $expected_effects: model.expected_effects
  $pip: model.get_pip()
  $study_pip: model.get_study_pip().values
  $pi: model.pi
  $active: model.active
  $credible_sets: model.credible_sets
  $purity: model.purity
  $params: params

fit_cafeh_summary_simple: fitting.py + Python(model = fit_cafeh_summary_simple(LD, B, S, K, p0k, fit); params = get_param_dict(model, compress=False))
  model: 'cafeh_summary_simple'
  LD: $LD
  B: $B
  S: $S
  K: $K
  p0k: 0.01, 0.1, 0.5
  fit: "weight_ard_active"
  update_variance: False
  $expected_effects: model.expected_effects
  $pip: model.get_pip()
  $study_pip: model.get_study_pip().values
  $pi: model.pi
  $active: model.active
  $credible_sets: model.credible_sets
  $purity: model.purity
  $params: params


fit_susie_summary: fitting.py + Python(results = fit_susie_summary(LD, B, S, K, p0k, fit))
  model: 'susie_summary'
  LD: $LD
  B: $B
  S: $S
  K: 5
  p0k: 1.0
  update_variance: False
  fit: "weight_ard_active"
  $expected_effects: results.expected_effects
  $study_pip: results.study_pip
  $credible_sets: results.credible_sets
  $purity: results.purity
  $params: results.params

fit_susie_summary_ss(fit_susie_summary):
  model: 'susie_summary_ss'
  p0k: 0.01, 0.1, 0.5

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
