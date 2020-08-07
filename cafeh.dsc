#!/usr/bin/env dsc

%include modules/data
%include modules/simulate
%include modules/fit

DSC:
  run:
    default: full_genotype
    cafeh_genotype: full_genotype * one_study_sim * fit_cafeh_genotype
    cafeh_summary: full_genotype * one_study_sim * genotype2ld * individual2summary * fit_cafeh_genotype
  exec_path: code
  output: output/cafeh-simulations