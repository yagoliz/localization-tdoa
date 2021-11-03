library(R.matlab)
library(ggplot2)

# Files for data analysis
folder = "~/desktop/phd/research/localization/core/localization-tdoa/matlab/results"
file = paste(folder,"results_comparison.mat",sep="/")

# Regular file loads and preparation
data = readMat(file)

dab = data$dab
dvbt = data$dvbt
gsm = data$gsm

## Kolmogorov tests: DAB struct
dab_dab_original = dab[[2]][[3]]
dab_dab_modified = dab[[3]][[3]]
dab_dvbt_original = dab[[5]][[3]]
dab_dvbt_modified = dab[[6]][[3]]
dab_lte_original = dab[[8]][[3]]
dab_lte_modified = dab[[9]][[3]]
dab_gsm_original = dab[[11]][[3]]
dab_gsm_modified = dab[[12]][[3]]

ks.test(dab_dab_original, dab_dab_modified)
ks.test(dab_dvbt_original, dab_dvbt_modified)
ks.test(dab_lte_original, dab_lte_modified)
ks.test(dab_gsm_original, dab_gsm_modified)

## Kolmogorov tests: DVB-T struct
dvbt_dab_original = dvbt[[2]][[3]]
dvbt_dab_modified = dvbt[[3]][[3]]
dvbt_dvbt_original = dvbt[[5]][[3]]
dvbt_dvbt_modified = dvbt[[6]][[3]]
dvbt_lte_original = dvbt[[8]][[3]]
dvbt_lte_modified = dvbt[[9]][[3]]
dvbt_gsm_original = dvbt[[11]][[3]]
dvbt_gsm_modified = dvbt[[12]][[3]]

ks.test(dvbt_dab_original, dvbt_dab_modified)
ks.test(dvbt_dvbt_original, dvbt_dvbt_modified)
ks.test(dvbt_lte_original, dvbt_lte_modified)
ks.test(dvbt_gsm_original, dvbt_gsm_modified)

## Kolmogorov tests: GSM struct
gsm_dab_original = gsm[[2]][[3]]
gsm_dab_modified = gsm[[3]][[3]]
gsm_dvbt_original = gsm[[5]][[3]]
gsm_dvbt_modified = gsm[[6]][[3]]
gsm_lte_original = gsm[[8]][[3]]
gsm_lte_modified = gsm[[9]][[3]]
gsm_gsm_original = gsm[[11]][[3]]
gsm_gsm_modified = gsm[[12]][[3]]

ks.test(gsm_dab_original, gsm_dab_modified)
ks.test(gsm_dvbt_original, gsm_dvbt_modified)
ks.test(gsm_lte_original, gsm_lte_modified)
ks.test(gsm_gsm_original, gsm_gsm_modified)