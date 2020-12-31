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
    ld_sim: ld_sim_coloc, ld_sim_no_coloc
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
    ld_small_pipeline: small_genotype * genotype2ld * ld_sim * individual2summary * (fit_cafeh_genotype_suggested * score_coloc_cafeh, coloc, caviar)
    ld_large_pipeline: full_genotype * ld_sim * individual2summary * (fit_cafeh_genotype_suggested * score_coloc_cafeh, coloc)

  exec_path: code
  output: output/cafeh-simulations


  # dsc cafeh.dsc --target vary_r2_pipeline block_sim_pipeline  --replicate 10 -o output --host marcc.yml

  # dsc cafeh.dsc --target "small_genotype * genotype2ld * block_study_sim * individual2summary * (fit_susie_summary_suggested, fit_cafeh_summary_suggested)" block_sim_pipeline  --truncate -o test_del --host marcc.yml