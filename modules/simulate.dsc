#---------------------------------------
# Simulation modules written by Gao Wang
# Used in SuSiE paper
# Implemented in Python
#---------------------------------------

n_causal_per_study_sim: simulation.py \
  + Python(sim=sim_n_causal_per_study(X, n_study, prop_colocalizing, n_causal_per_study, pve, effect_distribution))
  n_study: 2
  prop_colocalizing: 0.5
  n_causal_per_study: 3
  effect_distribution: 'normal'
  pve: 0.05, 0.1, 0.2, 0.4
  X: $X
  $sim: sim
  $Y: sim['expression']
  $K: sim['K']

single_sim(n_causal_per_study_sim):
  n_study: 20
  n_causal_per_study: 3
  pve: 0.1

one_study_sim(n_causal_per_study_sim):
  n_study: 1
  prop_colocalizing: 1.0
  n_causal_per_study: 1, 2, 3, 4, 5

two_study_sim(n_causal_per_study_sim):
  n_study: 2
  prop_colocalizing: 0.5
  n_causal_per_study: 1, 2, 3, 4, 5

increase_study_sim(n_causal_per_study_sim):
  n_study: 2, 5, 10, 20, 50
  n_causal_per_study: 1, 2, 3
  pve: 0.05, 0.1
