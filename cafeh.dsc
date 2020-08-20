#!/usr/bin/env dsc

%include modules/data
%include modules/simulate
%include modules/fit
%include modules/score

DSC:
  define:
    #simulations: tissue_specific_sim, general_sim, global_effect_sim
    simulations: general_sim
    simulate: small_genotype * genotype2ld * simulations * individual2summary

    # cafeh 
    fit_cafeh: fit_cafeh_genotype, fit_cafeh_summary, fit_cafeh_summary_simple

    # susie
    fit_susie: (fit_susie_genotype, fit_susie_summary, fit_susie_genotype_ss, fit_susie_summary_ss)

    # caviar
    caviar: fit_caviar * (score_finemapping_caviar, fit_ecaviar*score_coloc_ecaviar)

    # coloc
    coloc: fit_coloc * score_coloc_coloc

  run:
    default: simulate * (fit_cafeh * score_coloc_cafeh, fit_susie, caviar, coloc)
    run_simulations: simulate
    # cafeh_genotype: full_genotype * one_study_sim * fit_cafeh_genotype
    # cafeh_summary: full_genotype * one_study_sim * genotype2ld * individual2summary * fit_cafeh_genotype
    # susie_summary: full_genotype * two_study_sim * genotype2ld * individual2summary * fit_susie_summary
    # cafeh_summary_simple: full_genotype * two_study_sim * genotype2ld * individual2summary * fit_cafeh_summary_simple


  exec_path: code
  output: output/cafeh-simulations