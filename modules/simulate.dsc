#---------------------------------------
# Simulation modules written by Gao Wang
# Used in SuSiE paper
# Implemented in Python
#---------------------------------------

block_study_sim: simulation.py \
  + Python(sim = sim_block_study(X, afreq, ldscore, n_study, n_blocks, n_causal_per_block, block_p, pve, effect_distribution, *r2_range, *ldscore_range))
  # demonstrate relative performance at a range of settings
  X: $X
  afreq: $afreq
  ldscore: $ldscore

  n_study: 2, 10
  block_p: 0.0
  n_blocks: 2
  n_causal_per_block: 1, 2, 3
  pve: 0.1, 0.05, 0.01
  effect_distribution: 'normal'
  r2_range: (0, 0.8)
  ldscore_range: (-1e10, 1e10)
  $residual_variance: sim['residual_variance']
  $true_effects: sim['true_effects']
  $true_coloc: sim['true_coloc']
  $Y: sim['expression']
  $K: sim['K']

simple_coloc_sim(block_study_sim):
  n_study: 2
  n_blocks: 1
  pve: 0.05, 0.1
  n_causal_per_block: 1, 2, 3

simple_no_coloc_sim(block_study_sim):
  n_study: 2
  n_blocks: 2
  pve: 0.05, 0.1
  n_causal_per_block: 1, 2, 3

simple_sim(block_study_sim):
  n_study: 4
  n_blocks: 2
  pve: 0.01, 0.05, 0.1, 0.2
  n_causal_per_block: 1, 2, 3

allelic_het_sim(block_study_sim):
  n_study: 4
  n_blocks: 2
  pve: 0.2
  n_causal_per_block: 5, 10

increase_tissue_sim(block_study_sim):
  # demonstrate relative performance at a range of settings
  X: $X
  n_study: 4, 8, 16
  block_p: 0.0
  n_blocks: 2
  n_causal_per_block: 1, 2, 3
  pve: 0.05
  effect_distribution: 'normal'
  r2_range: (0, 0.8)
  $residual_variance: sim['residual_variance']
  $true_effects: sim['true_effects']
  $true_coloc: sim['true_coloc']
  $Y: sim['expression']
  $K: sim['K']

r2_between_blocks(block_study_sim):
  # demonstrate relative performance at a range of settings
  X: $X
  n_study: 4
  block_p: 0.0
  n_blocks: 2
  n_causal_per_block: 1
  pve: 0.05, 0.1, 0.2
  r2_range: (0, 0.5), (0.5, 0.7), (0.7, 0.9)


high_r2(block_study_sim):
  # demonstrate relative performance at a range of settings
  X: $X
  n_study: 4
  block_p: 0.0
  n_blocks: 2
  n_causal_per_block: 1
  pve: 0.05, 0.1, 0.2, 0.4
  r2_range: (0.6, 0.7), (0.7, 0.8), (0.8, 0.9), (0.9, 0.1)


increase_r2_sim(block_study_sim):
  # demonstrate relative performance at a range of settings
  rep: 1, 2, 3, 4, 5
  X: $X
  n_study: 2
  block_p: 0.0
  n_blocks: 2, 1
  n_causal_per_block: 1
  pve: 0.1
  r2_range: (0, 0.5), (0.5, 0.6), (0.6, 0.7), (0.7, 0.8), (0.8, 0.9)
