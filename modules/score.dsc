score_coloc_cafeh: scoring.py + Python(results = score_coloc_cafeh(active, true_coloc, thresh))
  thresh: 0.99
  active: $active
  true_coloc: $true_coloc
  $p_coloc: results['p_coloc']
  $TP: results['true_positive']
  $FP: results['false_positive']
  $TN: results['true_negative']
  $FN: results['false_negative']

score_coloc_coloc: scoring.py + Python(results = score_coloc_coloc(coloc_out, true_coloc, thresh))
  thresh: 0.99
  coloc_out: $coloc_out
  true_coloc: $true_coloc
  $p_coloc: results['p_coloc']
  $TP: results['true_positive']
  $FP: results['false_positive']
  $TN: results['true_negative']
  $FN: results['false_negative']

score_coloc_ecaviar: scoring.py + Python(results = score_coloc_ecaviar(ecaviar_out, true_coloc, thresh))
  thresh: 0.99
  ecaviar_out: $ecaviar_out
  true_coloc: $true_coloc
  $p_coloc: results['p_coloc']
  $TP: results['true_positive']
  $FP: results['false_positive']
  $TN: results['true_negative']
  $FN: results['false_negative']

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



