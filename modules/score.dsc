score_coloc_cafeh: scoring.py + Python(results = score_coloc_cafeh(active))
  active: $active
  $p_coloc: results['p_coloc']
  $maxmin: results['maxmin']

cafeh_elbo:
  fit_model: $fit_model
  $elbos: fit_model.elbos

score_coloc_coloc: scoring.py + Python(results = score_coloc_coloc(coloc_out))
  coloc_out: $coloc_out
  $p_coloc: results['p_coloc']

score_coloc_ecaviar: scoring.py + Python(results = score_coloc_ecaviar(ecaviar_out))
  ecaviar_out: $ecaviar_out
  $p_coloc: results['p_coloc']

score_coloc_finemap: scoring.py + Python(results = score_coloc_finemap(finemap_out))
  finemap_out: $finemap_out
  $p_coloc: results['p_coloc']

score_coloc_susie: scoring.py + Python(results = score_coloc_susie(study_pip))
  study_pip: $study_pip
  $p_coloc: results['p_coloc']

score_finemapping_cafeh: scoring.py + Python(results = score_finemapping_cafeh(credible_sets, purity, true_effects))
  credible_sets: $credible_sets
  purity: $purity
  true_effects: $true_effects
  $n_components: results['n_components']
  $n_components_with_causal: results['n_components_with_causal']
  $n_causal_in_cs: results['n_causal_in_cs']
  $n_causal: results['n_causal']
  $n_top_causal: results['n_top_causal']

score_finemapping_caviar: scoring.py + Python(results = score_finemapping_caviar(caviar_out, true_effects))
  caviar_out: $caviar_out
  true_effects: $true_effects
  $study_pip: results.study_pip
  $credible_sets: results.credible_sets
  $n_causal_in_cs: results.n_causal_in_cs
  $n_causal: results.n_causal


score_finemapping_finemap: scoring.py + Python(results = score_finemapping_finemap(finemap_out, true_effects))
  finemap_out: $finemap_out
  true_effects: $true_effects
  $study_pip: results.study_pip

