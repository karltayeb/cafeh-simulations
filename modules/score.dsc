score_coloc_cafeh: scoring.py + Python(results = score_coloc_cafeh(active, sim['true_effects'], thresh))
  thresh: 0.99
  active: $active
  sim: $sim
  $TP: results['true_positive']
  $FP: results['false_positive']
  $TN: results['true_negative']
  $FN: results['false_negative']


score_finemapping_cafeh: scoring.py + Python(results = score_finemapping_cafeh(credible_sets, purity, sim['true_effects']))
  credible_sets: $credible_sets
  purity: $purity
  sim: $sim
  $n_components: results['n_components']
  $n_components_with_causal: results['n_components_with_causal']
  $n_causal_in_cs: results['n_causal_in_cs']
  $n_causal: results['n_causal']
  $n_top_causal: results['n_top_causal']