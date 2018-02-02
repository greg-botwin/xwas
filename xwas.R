# load library
library(tidyverse)
library(readxl)

# cedars females only cd
# cedars males only cd
# international females only cd
# international males only cd

<<<<<<< HEAD
# read in all data, filter for ADD, compute: Z, effective sample size, weighted Z and weights squared, log p
=======
# read in all data, filter for ADD, compute: Z, effective sample size, weighted Z and weights squared
>>>>>>> bce594adefa0bd77c0036576072ffe6f8736e088
cedars_cd_female <- read_xlsx("data/Ichip1to6_BBC_chrXandY_clean_FINAL2_metal.CD.assoc.logistic.xlsx", 
                              col_types = c("numeric", "text", "numeric", "text",
                                            "text", "numeric", "numeric", "numeric", 
                                            "numeric", "numeric", "numeric", "numeric")) %>%
  filter(TEST == "ADD") %>%
  mutate(Z = BETA/SE) %>%
  mutate(effective_sample_size = NMISS) %>%
  mutate(weighted_z = Z * effective_sample_size) %>%
  mutate(weights_sq = effective_sample_size^2) %>%
  mutate(log_p = log(P))

cedars_cd_male <- read_xlsx("data/Ichip1to6_BBC_chrXandY_clean_FINAL2_male_metal.CD.assoc.logistic.xlsx", 
                            col_types = c("numeric", "text", "numeric", "text",
                                          "text", "numeric", "numeric", "numeric", 
                                          "numeric", "numeric", "numeric", "numeric")) %>%
  filter(TEST == "ADD") %>%
  select(-CHR, -BP, -A1, -TEST) %>%
  mutate(Z = BETA/SE) %>%
  mutate(effective_sample_size = NMISS/2) %>%
  mutate(weighted_z = Z * effective_sample_size) %>%
  mutate(weights_sq = effective_sample_size^2)%>%
  mutate(log_p = log(P))

int_cd_female <- read_xlsx("data/release5_IIBDgc_XandY_cleanFINAL_female_metal.CD.assoc.logistic.xlsx", 
                           col_types = c("numeric", "text", "numeric", "text",
                                         "text", "numeric", "numeric", "numeric", 
                                         "numeric", "numeric", "numeric", "numeric")) %>%
  filter(TEST == "ADD") %>%
  select(-CHR, -BP, -A1, -TEST) %>%
  mutate(Z = BETA/SE) %>%
  mutate(effective_sample_size = NMISS) %>%
  mutate(weighted_z = Z * effective_sample_size) %>%
  mutate(weights_sq = effective_sample_size^2) %>%
  mutate(log_p = log(P))

int_cd_male <-read_xlsx("data/release5_IIBDgc_XandY_cleanFINAL_male_metal.P2.assoc.logisticCD.xlsx", 
                        col_types = c("numeric", "text", "numeric", "text",
                                      "text", "numeric", "numeric", "numeric", 
                                      "numeric", "numeric", "numeric", "numeric")) %>%
  filter(TEST == "ADD") %>%
  select(-CHR, -BP, -A1, -TEST) %>%
  mutate(Z = BETA/SE) %>%
  mutate(effective_sample_size = NMISS/2) %>%
  mutate(weighted_z = Z * effective_sample_size) %>%
  mutate(weights_sq = effective_sample_size^2) %>%
  mutate(log_p = log(P))

# rename to specific colnames
colnames(cedars_cd_female) <- paste("cedars_cd_female", colnames(cedars_cd_female), sep = "_")
colnames(cedars_cd_male) <- paste("cedars_cd_male", colnames(cedars_cd_male), sep = "_")
colnames(int_cd_female) <- paste("int_cd_female", colnames(int_cd_female), sep = "_")
colnames(int_cd_male) <- paste("int_cd_male", colnames(int_cd_male), sep = "_")

# rename unique SNP to SNP
cedars_cd_female <- rename(cedars_cd_female, SNP = cedars_cd_female_SNP)
cedars_cd_male <- rename(cedars_cd_male, SNP = cedars_cd_male_SNP)
int_cd_female <- rename(int_cd_female, SNP = int_cd_female_SNP)
int_cd_male <- rename(int_cd_male, SNP = int_cd_male_SNP)

# combine all data 
combined_cedars_int_cd <- left_join(cedars_cd_female, cedars_cd_male, by = "SNP")
combined_cedars_int_cd <- left_join(combined_cedars_int_cd, int_cd_female, by = "SNP")
combined_cedars_int_cd <- left_join(combined_cedars_int_cd, int_cd_male, by = "SNP")

# compute p combined fisher
combined_cedars_int_cd <- combined_cedars_int_cd %>%
  mutate(p_comb_fisher = 1-pchisq(-2*(cedars_cd_female_log_p + cedars_cd_male_log_p +
                                             int_cd_female_log_p + int_cd_male_log_p), 8))
# compute p_comb_stouffer
combined_cedars_int_cd <- combined_cedars_int_cd %>%
  mutate(sum_weigthed_z = cedars_cd_female_weighted_z + cedars_cd_male_weighted_z +
           int_cd_female_weighted_z + int_cd_male_weighted_z) %>%
  mutate(sum_weights_sq = cedars_cd_female_weights_sq + cedars_cd_male_weights_sq +
           int_cd_female_weights_sq + int_cd_male_weights_sq) %>%
  mutate(sqrt_sum_weights_sq = sqrt(sum_weights_sq)) %>%
  mutate(p_comb_stouffer = 2*(1-pnorm((sum_weigthed_z/sqrt_sum_weights_sq))))

# write results to csv
combined_cedars_int_cd %>%
<<<<<<< HEAD
  write_csv("results/combined_cedars_int_cd.csv")

combined_cedars_int_cd %>%
  select(SNP, p_comb_fisher, p_comb_stouffer) %>%
  write_csv("results/results.csv")
=======
  write_csv("combined_cedars_int_cd.csv")

combined_cedars_int_cd %>%
  select(SNP, p_comb_fisher, p_comb_stouffer) %>%
  write_csv("results.csv")
>>>>>>> bce594adefa0bd77c0036576072ffe6f8736e088
