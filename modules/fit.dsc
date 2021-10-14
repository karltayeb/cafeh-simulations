fit_cafeh_genotype: fitting.py + Python(model = fit_cafeh_genotype(X.T, Y, K, p0k, w_prior_variance, standardize, update_ard, update_active, update_variance); params = get_param_dict(model))
  model: 'cafeh_genotype'
  X: $X
  Y: $Y
  K: $K

  p0k: 0.1
  standardize: True, False
  update_ard: True, False
  update_active: True
  update_variance: True, False
  w_prior_variance: 1, 0.1, 0.01, 0.001

  $pip: model.get_pip()
  $study_pip: model.get_study_pip().values
  $pi: model.pi
  $active: model.active
  $fit_model: model

fit_cafeh_genotype_suggested(fit_cafeh_genotype):
  p0k: 0.1
  standardize: True
  update_ard: True
  update_active: True
  update_variance: True
  w_prior_variance: 0.1


fit_cafeh_genotype_snp(fit_cafeh_genotype_suggested):
  X: $X_snp

fit_cafeh_genotype_sv(fit_cafeh_genotype_suggested):
  X: $X_sv

fit_cafeh_genotype_snp_sv(fit_cafeh_genotype_suggested):
  X: $X_snp_sv

fit_cafeh_genotype_conservative(fit_cafeh_genotype):
  model: 'cafeh_genotype_conservative'
  p0k: 0.01
  standardize: True
  update_ard: True
  update_active: True
  update_variance: True
  w_prior_variance: 0.001

fit_cafeh_param_grid(fit_cafeh_genotype):
  model: 'cafeh_genotype_param_grid'
  p0k: 0.001, 0.01, 0.1
  standardize: True
  update_ard: True
  update_active: True
  update_variance: True
  w_prior_variance: 0.0001, 0.001, 0.01, 0.1

fit_cafeh_genotype_fixed_var(fit_cafeh_genotype):
  p0k: 0.1
  standardize: True
  update_ard: True
  update_active: True
  update_variance: False
  w_prior_variance: 0.1


fit_cafeh_genotype_small_w_var(fit_cafeh_genotype):
  p0k: 0.1
  standardize: True
  update_ard: True
  update_active: True
  update_variance: True
  w_prior_variance: 0.001, 0.0001


fit_cafeh_genotype_no_ard(fit_cafeh_genotype):
  model: 'cafeh_genotype_no_ard'
  p0k: 0.1
  standardize: True
  update_ard: False
  update_active: True
  update_variance: True
  w_prior_variance: 0.1


fit_cafeh_genotype_small_p(fit_cafeh_genotype):
  model: 'cafeh_genotype_small_p0k'
  p0k: 0.001
  standardize: True
  update_ard: True
  update_active: True
  update_variance: True
  w_prior_variance: 0.1

fit_susie_genotype: fitting.py + Python(results = fit_susie_genotype(X.T, Y, K, p0k, w_prior_variance, standardize, update_ard, update_active, update_variance))
  model: 'susie_genotype'
  X: $X
  Y: $Y
  K: 10

  p0k: 1.0
  standardize: True
  update_ard: True
  update_active: False
  update_variance: True
  w_prior_variance: 0.2

  $expected_effects: results.expected_effects
  $study_pip: results.study_pip
  $credible_sets: results.credible_sets
  $purity: results.purity
  $params: results.params

fit_susie_genotype_ss(fit_susie_genotype):
  model: 'susie_genotype_ss'
  p0k: 0.1
  update_active: True

fit_cafeh_summary: fitting.py + Python(model = fit_cafeh_summary(LD, B, se, S, K, p0k, w_prior_variance, standardize, update_ard, update_active); params = get_param_dict(model))
  model: 'cafeh_summary'
  LD: $LD
  B: $B
  se: $se
  S: $S
  K: $K
  p0k: 0.1
  standardize: False
  update_ard: True
  update_active: True
  update_variance: False # does nothing if true
  w_prior_variance: 1.0

  $pip: model.get_pip()
  $study_pip: model.get_study_pip().values
  $pi: model.pi
  $active: model.active
  $m: model

fit_cafeh_summary_suggested(fit_cafeh_summary):
  p0k: 0.1
  standardize: True
  update_ard: True
  update_active: True
  update_variance: True

fit_susie_summary: fitting.py + Python(results = fit_susie_summary(LD, B, se, S, K, p0k, w_prior_variance, standardize, update_ard, update_active))
  model: 'susie_summary'
  LD: $LD
  B: $B
  se: $se
  S: $S
  K: 5

  p0k: 1.0
  standardize: False
  update_ard: True
  update_active: True
  update_variance: False # does nothing if true
  w_prior_variance: 1.0

  $expected_effects: results.expected_effects
  $study_pip: results.study_pip
  $credible_sets: results.credible_sets
  $purity: results.purity
  $params: results.params


fit_susie_summary_suggested(fit_susie_summary):
  p0k: 1.0
  standardize: True
  update_ard: True
  update_active: True
  update_variance: True


fit_cafeh_genotype_pairwise: fitting.py + Python(out = fit_cafeh_genotype_pairwise(X.T, Y, K, p0k, w_prior_variance, standardize, update_ard, update_active, update_variance))
  model: 'cafeh_genotype'
  X: $X
  Y: $Y
  K: $K
  p0k: 0.1
  standardize: True
  update_ard: True
  update_active: True
  update_variance: True
  w_prior_variance: 1.0

  $p_coloc: out['p_coloc']

fit_cafeh_summary_simple: fitting.py + Python(model = fit_cafeh_summary_simple(LD, B, se, S, K, p0k, w_prior_variance, standardize, fit, max_iter=max_iter); params = get_param_dict(model, compress=False))
  model: 'cafeh_summary_simple'
  LD: $LD
  B: $B
  se: $se
  S: $S
  K: $K
  p0k: 0.1
  standardize: True, False
  fit: "weight_ard_active", "weight_active"
  max_iter: 50
  update_variance: False
  w_prior_variance: 1.0

  $expected_effects: model.expected_effects
  $pip: model.get_pip()
  $study_pip: model.get_study_pip().values
  $pi: model.pi
  $active: model.active
  $credible_sets: model.credible_sets
  $purity: model.purity
  $params: params

fit_susie_summary_ss(fit_susie_summary):
  model: 'susie_summary_ss'
  p0k: 0.1

fit_caviar: fitting.py + Python(caviar_out = run_caviar(B, se, LD, z_filter))
  model: 'caviar'
  z_filter: 0
  LD: $LD
  B: $B
  se: $se
  $caviar_out: caviar_out

fit_caviar_z2(fit_caviar):
  model: 'caviar_z2'
  z_filter: 2

fit_caviar_z3(fit_caviar):
  model: 'caviar_z3'
  z_filter: 3

fit_ecaviar: fitting.py + Python(ecaviar_out = ecaviar_from_caviar(caviar_out))
  model: 'ecaviar'
  caviar_out: $caviar_out
  $ecaviar_out: ecaviar_out

fit_coloc: fitting.py + Python(coloc_out = run_coloc(B, se))
  model: 'coloc'
  B: $B
  se: $se
  $coloc_out: coloc_out

fit_finemap: fitting.py + Python(finemap_out = run_finemap(B, se, afreq, LD))
  model: 'finemap'
  LD: $LD
  B: $B
  se: $se
  afreq: $afreq
  $finemap_out: finemap_out
