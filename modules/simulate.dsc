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
  pve: 0.01, 0.05, 0.1, 0.2
  X: $X
  $sim: sim
  $true_effects: sim['true_effects']
  $true_coloc: sim['true_coloc']
  $Y: sim['expression']
  $K: sim['K']

tissue_specific_sim(n_causal_per_study_sim):
  # want to show that cafeh does not have a lot of false postives as tissues increase
  prop_colocalizing: 0.0
  n_study: 2, 5, 10
  n_causal_per_study: 1, 2, 3

general_sim(n_causal_per_study_sim):
  # demonstrate relative performance at a range of settings
  prop_colocalizing: 0.5
  n_study: 2, 5, 10
  n_causal_per_study: 1, 2, 3, 4

large_n_sim(n_causal_per_study_sim):
  # demonstrate relative performance at a range of settings
  prop_colocalizing: 0.5
  n_study: 20, 50
  n_causal_per_study: 1, 2, 3, 4, 5

global_effect_sim(n_causal_per_study_sim):
  # test at lower pve to show power increases w/ n_study
  prop_colocalizing: 1.0
  n_study: 2, 5, 10
  n_causal_per_study: 1, 2, 3
  pve: 0.005, 0.01, 0.05, 0.1