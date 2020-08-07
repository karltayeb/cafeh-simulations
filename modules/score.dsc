score_coloc_cafeh: scoring.py + Python(results = score_coloc(active, sim['true_effects'], thresh))
  thresh: 0.99
  active: $active
  sim: $sim
  $TP: results['true_positive']
  $FP: results['false_positive']
  $TN: results['true_negative']
  $FN: results['fasle_negative']

  $expected_effects: model.expected_effects
  $pi: model.pi
  $active: model.active
  $params: params
