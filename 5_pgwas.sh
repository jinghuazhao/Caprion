#!/usr/bin/bash

export ref=/rds/project/jmmh2/rds-jmmh2-post_qc_data/interval/reference_files/genetic/reference_files_genotyped_imputed/

seq 22 | \
xargs -I {} cut -f15 ${ref}/impute_{}_interval.snpstats | \
grep -v MAF | \
sed '/^$/d' | \
Rscript -e 'maf <- scan("stdin"); cutoffs <- cut(maf,c(0,0.001,0.01,0.05,0.1,0.2,0.3,0.4));table(cutoffs);table(cutoffs)/length(maf)'

# Read 87696888 items
# cutoffs
#    (0,0.001] (0.001,0.01]  (0.01,0.05]   (0.05,0.1]    (0.1,0.2]    (0.2,0.3]
#     63805560      9376909      2990789      1414448      1825305      1393653
#    (0.3,0.4]
#      1236541
# cutoffs
#    (0,0.001] (0.001,0.01]  (0.01,0.05]   (0.05,0.1]    (0.1,0.2]    (0.2,0.3]
#   0.72756926   0.10692408   0.03410371   0.01612883   0.02081379   0.01589170
#    (0.3,0.4]
#   0.01410017

# by chromosome
# seq 22 | xargs -I {} bash -c "cut -f15 ${ref}/impute_{}_interval.snpstats > {}.MAF"
