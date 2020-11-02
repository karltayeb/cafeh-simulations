#!/usr/bin/env dsc

%include modules/data
%include modules/simulate
%include modules/fit
%include modules/score

DSC:
  define:
    #simulations: tissue_specific_sim, general_sim, global_effect_sim
    get_genotype: small_genotype, small_genotype_random
    simulations: block_study_sim
    simulate: get_genotype * genotype2ld * simulations * individual2summary 

    # cafeh 
    fit_cafeh: fit_cafeh_genotype #, fit_cafeh_summary, fit_cafeh_summary_simple

    # susie
    fit_susie: (fit_susie_genotype, fit_susie_summary)

    # caviar
    caviar: fit_caviar * (score_finemapping_caviar, fit_ecaviar*score_coloc_ecaviar)

    # coloc
    coloc: fit_coloc * score_coloc_coloc

  run:
    #default: simulate * (fit_cafeh * score_coloc_cafeh, caviar, coloc)
    vary_r2_pipeline: small_genotype * genotype2ld * r2_between_blocks * individual2summary * (fit_cafeh * score_coloc_cafeh, coloc)
    block_sim_pipeline: small_genotype * genotype2ld * block_study_sim * individual2summary * (fit_cafeh_genotype_suggested * score_coloc_cafeh, coloc, caviar, fit_susie_genotype_suggested)
    one_block_pipeline: small_genotype * genotype2ld * one_block_sim * individual2summary * (fit_cafeh_genotype_suggested, fit_susie_genotype_suggested, fit_caviar * score_finemapping_caviar)
    increase_tissue_pipeline: small_genotype * genotype2ld * increase_tissue_sim * individual2summary * (fit_cafeh_genotype_suggested * score_coloc_cafeh, coloc, fit_susie_genotype_suggested)
    #default: simulate * (fit_cafeh * score_coloc_cafeh, fit_susie, caviar, coloc)
    # cafeh_genotype: full_genotype * one_study_sim * fit_cafeh_genotype
    # cafeh_summary: full_genotype * one_study_sim * genotype2ld * individual2summary * fit_cafeh_genotype
    # susie_summary: full_genotype * two_study_sim * genotype2ld * individual2summary * fit_susie_summary
    # cafeh_summary_simple: full_genotype * two_study_sim * genotype2ld * individual2summary * fit_cafeh_summary_simple


  exec_path: code
  output: output/cafeh-simulations


  # dsc cafeh.dsc --target vary_r2_pipeline block_sim_pipeline  --replicate 10 -o output --host marcc.yml

  # dsc cafeh.dsc --target "small_genotype * genotype2ld * block_study_sim * individual2summary * (fit_susie_summary_suggested, fit_cafeh_summary_suggested)" block_sim_pipeline  --truncate -o test_del --host marcc.yml