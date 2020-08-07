capture: Python(pickle.dump(X, open('./tmp/X', 'wb')); pickle.dump(Y, open('./tmp/Y', 'wb')); a=X)
  @CONF: python_modules = (numpy)
  X: $X
  Y: $Y
  $a: a

fit_cafeh_genotype: fitting.py + Python(model = fit_cafeh_genotype(X=X.T, Y=Y, K=10); params = get_param_dict(model))
  X: $X
  Y: $Y
  $expected_effects: model.expected_effects
  $pi: model.pi
  $active: model.active
  $params: params


fit_cafeh_summary: fitting.py + Python(model = fit_cafeh_summary(LD, B, S, K=10); params = get_param_dict(model))
  LD: $LD
  B: $B
  S: $S
  $expected_effects: model.expected_effects
  $pi: model.pi
  $active: model.active
  $params: params