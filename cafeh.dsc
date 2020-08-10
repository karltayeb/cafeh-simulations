#!/usr/bin/env dsc

%include modules/data
%include modules/simulate
%include modules/fit
%include modules/score

DSC:
  define:
    simulate: small_genotype * genotype2ld * (one_study_sim, two_study_sim, increase_study_sim) * individual2summary
  run:
    default: simulate * fit_cafeh_genotype * score_coloc_cafeh
    summary: simulate * fit_cafeh_summary * score_coloc_cafeh
    caviar: small_genotype * genotype2ld * (one_study_sim, two_study_sim) * individual2summary * fit_caviar
    ecaviar: small_genotype * genotype2ld * (one_study_sim, two_study_sim) * individual2summary * fit_caviar * fit_ecaviar
    coloc: small_genotype * genotype2ld * two_study_sim * individual2summary * fit_coloc

    cafeh_genotype: full_genotype * one_study_sim * fit_cafeh_genotype
    cafeh_summary: full_genotype * one_study_sim * genotype2ld * individual2summary * fit_cafeh_genotype
  exec_path: code
  output: output/cafeh-simulations