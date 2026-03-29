# Best model
# fit H1 to H6 with the model comparision
required_packages <- c("tidyverse", "metafor", "clubSandwich", "MAd", "weightr", "puniform", "metapower", "rempsyc","bruceR")

# Check and install missing packages
install_and_Load <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE) # Load the required packages
  }
}

install_and_Load(required_packages)

# reading the data
MA_dat <- read_csv("Coding Sample - Final Meta Data.csv", skip = 1)
MA_dat$`Fairness/Cooperativeness A Binary` <- as.numeric(MA_dat$`Fairness/Cooperativeness A Binary`)

# Unfair Condition Analysis ----
## H1 The Effects of Violator’s Identity on Punishment ----
# fair_condition <- c("Fair", "Cooperate", "Neutral", "Generous")
# unfair_condition <- c("Unfair", "Defect", "Selfish", "Cheat")
H1_data <- MA_dat %>%
  filter(
    `Fairness/Cooperativeness A Binary` == 0, # unfair condition
    `Person A affiliation` %in% c("Ingroup", "Outgroup")
  ) %>%
  select(
    No., Study, `Person A affiliation`, `Recipient/Player B Affliation`, `Sample ID`, `IV 1 Levels`, `Age (mean)`,
    `Fairness/Cooperativeness of Person A`, `Fairness/Cooperativeness of Person B`,
    `P_Average Invested Amount`, `P_SD invested amount`, `N by condition`, `Total Trials`, `Intext-Citation`,
    Country, `Punishment-Cost Ratio`, `Lab or Online`, `Group Manipulation`, Female, Male, `Game Type`, `Relative Endowments`, `Real or Deception`, `Student Sample`
  ) %>%
  filter(!is.na(`P_Average Invested Amount`)) %>% # Exclude mean is NA
  pivot_wider(
    names_from = `Person A affiliation`,
    values_from = c(`P_Average Invested Amount`, `P_SD invested amount`, `N by condition`,`Total Trials`)
  ) %>% 
  mutate(`P_Average Invested Amount_Outgroup` = as.numeric(`P_Average Invested Amount_Outgroup`),
         `P_Average Invested Amount_Ingroup` = as.numeric(`P_Average Invested Amount_Ingroup`),
         `P_SD invested amount_Outgroup` = as.numeric(`P_SD invested amount_Outgroup`),
         `P_SD invested amount_Ingroup` = as.numeric(`P_SD invested amount_Ingroup`),
         `N by condition_Outgroup` = as.numeric(`N by condition_Outgroup`),
         `N by condition_Ingroup` = as.numeric(`N by condition_Ingroup`)
  ) %>%
  escalc(
    measure = "SMDH", # Hedges' g
    m1i = `P_Average Invested Amount_Outgroup`, sd1i = `P_SD invested amount_Outgroup`, n1i = `N by condition_Outgroup`,
    m2i = `P_Average Invested Amount_Ingroup`,  sd2i = `P_SD invested amount_Ingroup`,  n2i = `N by condition_Ingroup`,
    data = .,
  ) %>% 
  group_by(No., Study) %>%
  mutate(
    effect_n = row_number(),
    effectID = paste(No., Study, effect_n, sep = "_")
  ) %>% 
  ungroup() %>%
  filter(!is.na(yi)) # Exclude effect size is NA

# deal with between data (assign temporary sample ID)
H1_Between_Data <- MA_dat %>% 
  filter(
    `Between or Within?...38` == "Between",
    `Fairness/Cooperativeness A Binary` == 0, # unfair condition
    `Person A affiliation` %in% c("Ingroup", "Outgroup")
  ) %>%
  select(
    No., Study, `Person A affiliation`, `Recipient/Player B Affliation`, `IV 1 Levels`, `Age (mean)`,,
    `Fairness/Cooperativeness of Person A`, `Fairness/Cooperativeness of Person B`,
    `P_Average Invested Amount`, `P_SD invested amount`, `N by condition`, `Total Trials`, `Intext-Citation`,
    Country, `Punishment-Cost Ratio`, `Lab or Online`, `Group Manipulation`, Female, Male, `Game Type`, `Relative Endowments`, `Real or Deception`, `Student Sample`
  ) %>%
  filter(!is.na(`P_Average Invested Amount`)) %>% # Exclude mean is NA
  group_by(No., Study,`Person A affiliation`) %>% 
  mutate(
    ID_n = row_number(),
    `Sample ID` = paste(No., Study, ID_n, sep = "_")
  ) %>% 
  select(-ID_n) %>% 
  ungroup() %>% 
  pivot_wider(
    names_from = `Person A affiliation`,
    values_from = c(`P_Average Invested Amount`, `P_SD invested amount`, `N by condition`,`Total Trials`)
  ) %>% 
  mutate(`P_Average Invested Amount_Outgroup` = as.numeric(`P_Average Invested Amount_Outgroup`),
         `P_Average Invested Amount_Ingroup` = as.numeric(`P_Average Invested Amount_Ingroup`),
         `P_SD invested amount_Outgroup` = as.numeric(`P_SD invested amount_Outgroup`),
         `P_SD invested amount_Ingroup` = as.numeric(`P_SD invested amount_Ingroup`),
         `N by condition_Outgroup` = as.numeric(`N by condition_Outgroup`),
         `N by condition_Ingroup` = as.numeric(`N by condition_Ingroup`)
  ) %>%
  escalc(
    measure = "SMDH", # Hedges' g
    m1i = `P_Average Invested Amount_Outgroup`, sd1i = `P_SD invested amount_Outgroup`, n1i = `N by condition_Outgroup`,
    m2i = `P_Average Invested Amount_Ingroup`,  sd2i = `P_SD invested amount_Ingroup`,  n2i = `N by condition_Ingroup`,
    data = .,
  ) %>% 
  group_by(No., Study) %>%
  mutate(
    effect_n = row_number(),
    effectID = paste(No., Study, effect_n, sep = "_")
  ) %>% 
  ungroup() %>%
  filter(!is.na(yi)) # Exclude effect size is NA

H1_data <- rbind(H1_data, H1_Between_Data)

V_H1 <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H1_data)
H1_model_1 <- rma.mv(yi, V_H1, random = ~ 1 | No. / Study / effectID, data = H1_data, method = "ML") # use ML for better model comparation
H1_model_2 <- rma.mv(yi, V_H1, random = ~ 1 | No. / effectID, data = H1_data, method = "ML")
H1_model_3 <- rma.mv(yi, V_H1, random = ~ 1 | effectID, data = H1_data, method = "ML")

H1_compare_1_2 <- anova(H1_model_1, H1_model_2)
H1_compare_1_3 <- anova(H1_model_1, H1_model_3)

## creat a table of model comparison
H1_model_compare <- data.frame(
  Model = c("Four-level model", 
            "Three-level model", 
            "Two-level model"),
  df = c(H1_compare_1_2$parms.f, 
         H1_compare_1_2$parms.r, 
         H1_compare_1_3$parms.r),
  BIC = c(H1_compare_1_2$fit.stats.f[4], 
          H1_compare_1_2$fit.stats.r[4], 
          H1_compare_1_3$fit.stats.r[4]), 
  LL = c(H1_compare_1_2$fit.stats.f[1], 
         H1_compare_1_2$fit.stats.r[1], 
         H1_compare_1_3$fit.stats.r[1]), 
  LRT = c("", 
          round(H1_compare_1_2$LRT, 3), 
          round(H1_compare_1_3$LRT, 3)), 
  p = c(NA, 
        H1_compare_1_2$pval, 
        H1_compare_1_3$pval))

MC_H1 <- H1_model_compare %>% 
  nice_table(title = "Model Comparison for H1")

# refit model with REML
H1_model_f <- rma.mv(yi, V_H1, random = ~ 1 | No. / Study / effectID, data = H1_data, method = "REML") # use REML for better fitting
robust(H1_model_f, cluster = No., clubSandwich=TRUE)

# checking likelihood profile plots
profile(H1_model_f)

# calculating the overall amount of heterogeneity
sum(H1_model_f$sigma2)

# obtaining I2 value
W <- diag(1/H1_data$vi)
X <- model.matrix(H1_model_f)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
H1_I_t <- 100 * sum(H1_model_f$sigma2) / (sum(H1_model_f$sigma2) + (H1_model_f$k-H1_model_f$p)/sum(diag(P))) # overall I2
H1_I <- 100 * H1_model_f$sigma2 / (sum(H1_model_f$sigma2) + (H1_model_f$k-H1_model_f$p)/sum(diag(P)))  # I2 per levels

## H2 (Deterrence perspective): ----
H2_data <- MA_dat %>%
  filter(
    `Fairness/Cooperativeness A Binary` == 0, # unfair condition
    `Recipient/Player B Affliation` %in% c("Ingroup", "Outgroup")
  ) %>%
  # group_by(No., Study, `Recipient/Player B Affliation`) %>%
  select(
    No., Study, `Person A affiliation`, `Recipient/Player B Affliation`, `Sample ID`, `IV 1 Levels`, `Age (mean)`,,
    `Fairness/Cooperativeness of Person A`, `Fairness/Cooperativeness of Person B`,
    `P_Average Invested Amount`, `P_SD invested amount`, `N by condition`, `Total Trials`, `Intext-Citation`,
    Country, `Punishment-Cost Ratio`, `Lab or Online`, `Group Manipulation`, Female, Male, `Game Type`, `Relative Endowments`, `Real or Deception`, `Student Sample`
  ) %>%
  filter(!is.na(`P_Average Invested Amount`)) %>% # Exclude mean is NA
  pivot_wider(
    names_from = `Recipient/Player B Affliation`,
    values_from = c(`P_Average Invested Amount`, `P_SD invested amount`, `N by condition`,`Total Trials`)
  ) %>% 
  mutate(`P_Average Invested Amount_Outgroup` = as.numeric(`P_Average Invested Amount_Outgroup`),
         `P_Average Invested Amount_Ingroup` = as.numeric(`P_Average Invested Amount_Ingroup`),
         `P_SD invested amount_Outgroup` = as.numeric(`P_SD invested amount_Outgroup`),
         `P_SD invested amount_Ingroup` = as.numeric(`P_SD invested amount_Ingroup`),
         `N by condition_Outgroup` = as.numeric(`N by condition_Outgroup`),
         `N by condition_Ingroup` = as.numeric(`N by condition_Ingroup`)
  ) %>% 
  escalc(
    measure = "SMDH", # Hedges' g
    m1i = `P_Average Invested Amount_Outgroup`, sd1i = `P_SD invested amount_Outgroup`, n1i = `N by condition_Outgroup`,
    m2i = `P_Average Invested Amount_Ingroup`,  sd2i = `P_SD invested amount_Ingroup`,  n2i = `N by condition_Ingroup`,
    data = .,
  ) %>% 
  group_by(No., Study) %>%
  mutate(
    effect_n = row_number(),
    effectID = paste(No., Study, effect_n, sep = "_")
  ) %>% 
  ungroup() %>%
  filter(!is.na(yi)) # Exclude effect size is NA

# deal with between data (assign temporary sample ID)
H2_Between_Data <- MA_dat %>% 
  filter(
    `Between or Within?...40` == "Between",
    `Fairness/Cooperativeness A Binary` == 0, # unfair condition
    `Recipient/Player B Affliation` %in% c("Ingroup", "Outgroup")
  ) %>%
  select(
    No., Study, `Person A affiliation`, `Recipient/Player B Affliation`, `IV 1 Levels`, `Age (mean)`,,
    `Fairness/Cooperativeness of Person A`, `Fairness/Cooperativeness of Person B`,
    `P_Average Invested Amount`, `P_SD invested amount`, `N by condition`, `Total Trials`, `Intext-Citation`,
    Country, `Punishment-Cost Ratio`, `Lab or Online`, `Group Manipulation`, Female, Male, `Game Type`, `Relative Endowments`, `Real or Deception`, `Student Sample`
  ) %>%
  filter(!is.na(`P_Average Invested Amount`)) %>% # Exclude mean is NA
  group_by(No., Study,`Recipient/Player B Affliation`) %>% 
  mutate(
    ID_n = row_number(),
    `Sample ID` = paste(No., Study, ID_n, sep = "_")
  ) %>% 
  select(-ID_n) %>% 
  ungroup() %>% 
  pivot_wider(
    names_from = `Recipient/Player B Affliation`,
    values_from = c(`P_Average Invested Amount`, `P_SD invested amount`, `N by condition`,`Total Trials`)
  ) %>% 
  mutate(`P_Average Invested Amount_Outgroup` = as.numeric(`P_Average Invested Amount_Outgroup`),
         `P_Average Invested Amount_Ingroup` = as.numeric(`P_Average Invested Amount_Ingroup`),
         `P_SD invested amount_Outgroup` = as.numeric(`P_SD invested amount_Outgroup`),
         `P_SD invested amount_Ingroup` = as.numeric(`P_SD invested amount_Ingroup`),
         `N by condition_Outgroup` = as.numeric(`N by condition_Outgroup`),
         `N by condition_Ingroup` = as.numeric(`N by condition_Ingroup`)
  ) %>%
  escalc(
    measure = "SMDH", # Hedges' g
    m1i = `P_Average Invested Amount_Outgroup`, sd1i = `P_SD invested amount_Outgroup`, n1i = `N by condition_Outgroup`,
    m2i = `P_Average Invested Amount_Ingroup`,  sd2i = `P_SD invested amount_Ingroup`,  n2i = `N by condition_Ingroup`,
    data = .,
  ) %>% 
  group_by(No., Study) %>%
  mutate(
    effect_n = row_number(),
    effectID = paste(No., Study, effect_n, sep = "_")
  ) %>% 
  ungroup() %>%
  filter(!is.na(yi)) # Exclude effect size is NA

H2_data <- rbind(H2_data, H2_Between_Data)

V_H2 <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H2_data)
H2_model_1 <- rma.mv(yi, V_H2, random = ~ 1 | No. / Study / effectID, data = H2_data, method = "ML")
H2_model_2 <- rma.mv(yi, V_H2, random = ~ 1 | No. / effectID, data = H2_data, method = "ML")
H2_model_3 <- rma.mv(yi, V_H2, random = ~ 1 | effectID, data = H2_data, method = "ML")

H2_compare_1_2 <- anova(H2_model_1, H2_model_2)
H2_compare_1_3 <- anova(H2_model_1, H2_model_3)
H2_compare_2_3 <- anova(H2_model_2, H2_model_3)

## creat a table of model comparison
H2_model_compare <- data.frame(
  Model = c("Four-level model", 
            "Three-level model", 
            "Two-level model"),
  df = c(H2_compare_1_2$parms.f, 
         H2_compare_1_2$parms.r, 
         H2_compare_1_3$parms.r),
  BIC = c(H2_compare_1_2$fit.stats.f[4], 
          H2_compare_1_2$fit.stats.r[4], 
          H2_compare_1_3$fit.stats.r[4]), 
  LL = c(H2_compare_1_2$fit.stats.f[1], 
         H2_compare_1_2$fit.stats.r[1], 
         H2_compare_1_3$fit.stats.r[1]), 
  LRT = c("", 
          round(H2_compare_1_2$LRT, 3), 
          round(H2_compare_1_3$LRT, 3)), 
  p = c(NA, 
        H2_compare_1_2$pval, 
        H2_compare_1_3$pval))

MC_H2 <- H2_model_compare %>% 
  nice_table(title = "Model Comparison for H2")

# refit model with REML
H2_model_f <- rma.mv(yi, V_H2, random = ~ 1 | effectID, data = H2_data, method = "REML")

robust(H2_model_f, cluster = No., clubSandwich=TRUE)

# checking likelihood profile plots
profile(H2_model_f)

# calculating the overall amount of heterogeneity
sum(H2_model_f$sigma2)

# obtaining I2 value
W <- diag(1/H2_data$vi)
X <- model.matrix(H2_model_f)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
H2_I_t <- 100 * sum(H2_model_f$sigma2) / (sum(H2_model_f$sigma2) + (H2_model_f$k-H2_model_f$p)/sum(diag(P))) # overall I2

# Fair Condition Analysis ################
## H5 The Effects of Violator’s Identity on Punishment ----
H5_data <- MA_dat %>%
  filter(
    `Fairness/Cooperativeness A Binary` == 1, # fair condition
    `Person A affiliation` %in% c("Ingroup", "Outgroup")
  ) %>%
  select(
    No., Study, `Person A affiliation`, `Recipient/Player B Affliation`, `Sample ID`, `IV 1 Levels`, `Age (mean)`,,
    `Fairness/Cooperativeness of Person A`, `Fairness/Cooperativeness of Person B`,
    `P_Average Invested Amount`, `P_SD invested amount`, `N by condition`, `Total Trials`, `Intext-Citation`,
    Country, `Punishment-Cost Ratio`, `Lab or Online`, `Group Manipulation`, Female, Male, `Game Type`, `Relative Endowments`, `Real or Deception`, `Student Sample`
  ) %>%
  filter(!is.na(`P_Average Invested Amount`)) %>% # Exclude mean is NA
  pivot_wider(
    names_from = `Person A affiliation`,
    values_from = c(`P_Average Invested Amount`, `P_SD invested amount`, `N by condition`,`Total Trials`)
  ) %>%
  mutate(`P_Average Invested Amount_Outgroup` = as.numeric(`P_Average Invested Amount_Outgroup`),
         `P_Average Invested Amount_Ingroup` = as.numeric(`P_Average Invested Amount_Ingroup`),
         `P_SD invested amount_Outgroup` = as.numeric(`P_SD invested amount_Outgroup`),
         `P_SD invested amount_Ingroup` = as.numeric(`P_SD invested amount_Ingroup`),
         `N by condition_Outgroup` = as.numeric(`N by condition_Outgroup`),
         `N by condition_Ingroup` = as.numeric(`N by condition_Ingroup`)
  ) %>%
  escalc(
    measure = "SMDH", # Hedges' g
    m1i = `P_Average Invested Amount_Outgroup`, sd1i = `P_SD invested amount_Outgroup`, n1i = `N by condition_Outgroup`,
    m2i = `P_Average Invested Amount_Ingroup`,  sd2i = `P_SD invested amount_Ingroup`,  n2i = `N by condition_Ingroup`,
    data = .,
  ) %>% 
  group_by(No., Study) %>%
  mutate(
    effect_n = row_number(),
    effectID = paste(No., Study, effect_n, sep = "_")
  ) %>% 
  ungroup() %>%
  filter(!is.na(yi)) # Exclude effect size is NA

H5_Between_Data <- MA_dat %>% 
  filter(
    `Between or Within?...38` == "Between",
    `Fairness/Cooperativeness A Binary` == 1, # fair condition
    `Person A affiliation` %in% c("Ingroup", "Outgroup")
  ) %>%
  select(
    No., Study, `Person A affiliation`, `Recipient/Player B Affliation`, `IV 1 Levels`, `Age (mean)`,,
    `Fairness/Cooperativeness of Person A`, `Fairness/Cooperativeness of Person B`,
    `P_Average Invested Amount`, `P_SD invested amount`, `N by condition`, `Total Trials`, `Intext-Citation`,
    Country, `Punishment-Cost Ratio`, `Lab or Online`, `Group Manipulation`, Female, Male, `Game Type`, `Relative Endowments`, `Real or Deception`, `Student Sample`
  ) %>%
  filter(!is.na(`P_Average Invested Amount`)) %>% # Exclude mean is NA
  group_by(No., Study,`Person A affiliation`) %>% 
  mutate(
    ID_n = row_number(),
    `Sample ID` = paste(No., Study, ID_n, sep = "_")
  ) %>% 
  select(-ID_n) %>% 
  ungroup() %>% 
  pivot_wider(
    names_from = `Person A affiliation`,
    values_from = c(`P_Average Invested Amount`, `P_SD invested amount`, `N by condition`,`Total Trials`)
  ) %>% 
  mutate(`P_Average Invested Amount_Outgroup` = as.numeric(`P_Average Invested Amount_Outgroup`),
         `P_Average Invested Amount_Ingroup` = as.numeric(`P_Average Invested Amount_Ingroup`),
         `P_SD invested amount_Outgroup` = as.numeric(`P_SD invested amount_Outgroup`),
         `P_SD invested amount_Ingroup` = as.numeric(`P_SD invested amount_Ingroup`),
         `N by condition_Outgroup` = as.numeric(`N by condition_Outgroup`),
         `N by condition_Ingroup` = as.numeric(`N by condition_Ingroup`)
  ) %>%
  escalc(
    measure = "SMDH", # Hedges' g
    m1i = `P_Average Invested Amount_Outgroup`, sd1i = `P_SD invested amount_Outgroup`, n1i = `N by condition_Outgroup`,
    m2i = `P_Average Invested Amount_Ingroup`,  sd2i = `P_SD invested amount_Ingroup`,  n2i = `N by condition_Ingroup`,
    data = .,
  ) %>% 
  group_by(No., Study) %>%
  mutate(
    effect_n = row_number(),
    effectID = paste(No., Study, effect_n, sep = "_")
  ) %>% 
  ungroup() %>%
  filter(!is.na(yi)) # Exclude effect size is NA

H5_data <- rbind(H5_data, H5_Between_Data)

V_H5 <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H5_data)
H5_model_1 <- rma.mv(yi, V_H5, random = ~ 1 | No. / Study / effectID, data = H5_data, method = "ML")
H5_model_2 <- rma.mv(yi, V_H5, random = ~ 1 | No. / effectID, data = H5_data, method = "ML")
H5_model_3 <- rma.mv(yi, V_H5, random = ~ 1 | effectID, data = H5_data, method = "ML")

H5_compare_1_2 <- anova(H5_model_1, H5_model_2)
H5_compare_1_3 <- anova(H5_model_1, H5_model_3)
H5_compare_2_3 <- anova(H5_model_2, H5_model_3)

## creat a table of model comparison
H5_model_compare <- data.frame(
  Model = c("Four-level model", 
            "Three-level model", 
            "Two-level model"),
  df = c(H5_compare_1_2$parms.f, 
         H5_compare_1_2$parms.r, 
         H5_compare_1_3$parms.r),
  BIC = c(H5_compare_1_2$fit.stats.f[4], 
          H5_compare_1_2$fit.stats.r[4], 
          H5_compare_1_3$fit.stats.r[4]), 
  LL = c(H5_compare_1_2$fit.stats.f[1], 
         H5_compare_1_2$fit.stats.r[1], 
         H5_compare_1_3$fit.stats.r[1]), 
  LRT = c("", 
          round(H5_compare_1_2$LRT, 3), 
          round(H5_compare_1_3$LRT, 3)), 
  p = c(NA, 
        H5_compare_1_2$pval, 
        H5_compare_1_3$pval))

MC_H5 <- H5_model_compare %>% 
  nice_table(title = "Model Comparison for H5")

# refit model with REML
H5_model_f <- rma.mv(yi, V_H5, random = ~ 1 | No. / effectID, data = H5_data, method = "REML")

robust(H5_model_f, cluster = No., clubSandwich=TRUE)

# calculating the overall amount of heterogeneity
sum(H5_model_f$sigma2)

# obtaining I2 value
W <- diag(1/H5_data$vi)
X <- model.matrix(H5_model_f)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
H5_I_t <- 100 * sum(H5_model_f$sigma2) / (sum(H5_model_f$sigma2) + (H5_model_f$k-H5_model_f$p)/sum(diag(P))) # overall I2
H5_I <- 100 * H5_model_f$sigma2 / (sum(H5_model_f$sigma2) + (H5_model_f$k-H5_model_f$p)/sum(diag(P)))  # I2 per levels

## H6 The Effects of Victim’s Identity on Punishment ----
H6_data <- MA_dat %>%
  filter(
    `Fairness/Cooperativeness A Binary` == 1, # fair condition
    `Recipient/Player B Affliation` %in% c("Ingroup", "Outgroup")
  ) %>%
  # group_by(No., Study, `Recipient/Player B Affliation`) %>%
  select(
    No., Study, `Person A affiliation`, `Recipient/Player B Affliation`, `Sample ID`, `IV 1 Levels`, `Age (mean)`,,
    `Fairness/Cooperativeness of Person A`, `Fairness/Cooperativeness of Person B`,
    `P_Average Invested Amount`, `P_SD invested amount`, `N by condition`, `Total Trials`, `Intext-Citation`,
    Country, `Punishment-Cost Ratio`, `Lab or Online`, `Group Manipulation`, Female, Male, `Game Type`, `Relative Endowments`, `Real or Deception`, `Student Sample`
  ) %>%
  filter(!is.na(`P_Average Invested Amount`)) %>% # Exclude mean is NA
  pivot_wider(
    names_from = `Recipient/Player B Affliation`,
    values_from = c(`P_Average Invested Amount`, `P_SD invested amount`, `N by condition`,`Total Trials`)
  ) %>% 
  mutate(`P_Average Invested Amount_Outgroup` = as.numeric(`P_Average Invested Amount_Outgroup`),
         `P_Average Invested Amount_Ingroup` = as.numeric(`P_Average Invested Amount_Ingroup`),
         `P_SD invested amount_Outgroup` = as.numeric(`P_SD invested amount_Outgroup`),
         `P_SD invested amount_Ingroup` = as.numeric(`P_SD invested amount_Ingroup`),
         `N by condition_Outgroup` = as.numeric(`N by condition_Outgroup`),
         `N by condition_Ingroup` = as.numeric(`N by condition_Ingroup`)
  ) %>%
  escalc(
    measure = "SMDH", # Hedges' g
    m1i = `P_Average Invested Amount_Outgroup`, sd1i = `P_SD invested amount_Outgroup`, n1i = `N by condition_Outgroup`,
    m2i = `P_Average Invested Amount_Ingroup`,  sd2i = `P_SD invested amount_Ingroup`,  n2i = `N by condition_Ingroup`,
    data = .,
  ) %>% 
  group_by(No., Study) %>%
  mutate(
    effect_n = row_number(),
    effectID = paste(No., Study, effect_n, sep = "_")
  ) %>% 
  ungroup() %>%
  filter(!is.na(yi)) # Exclude effect size is NA

H6_Between_Data <- MA_dat %>% 
  filter(
    `Between or Within?...40` == "Between",
    `Fairness/Cooperativeness A Binary` == 1, # fair condition
    `Recipient/Player B Affliation` %in% c("Ingroup", "Outgroup")
  ) %>%
  select(
    No., Study, `Person A affiliation`, `Recipient/Player B Affliation`, `IV 1 Levels`, `Age (mean)`,,
    `Fairness/Cooperativeness of Person A`, `Fairness/Cooperativeness of Person B`,
    `P_Average Invested Amount`, `P_SD invested amount`, `N by condition`, `Total Trials`, `Intext-Citation`,
    Country, `Punishment-Cost Ratio`, `Lab or Online`, `Group Manipulation`, Female, Male, `Game Type`, `Relative Endowments`, `Real or Deception`, `Student Sample`
  ) %>%
  filter(!is.na(`P_Average Invested Amount`)) %>% # Exclude mean is NA
  group_by(No., Study,`Recipient/Player B Affliation`) %>% 
  mutate(
    ID_n = row_number(),
    `Sample ID` = paste(No., Study, ID_n, sep = "_")
  ) %>% 
  select(-ID_n) %>% 
  ungroup() %>% 
  pivot_wider(
    names_from = `Recipient/Player B Affliation`,
    values_from = c(`P_Average Invested Amount`, `P_SD invested amount`, `N by condition`,`Total Trials`)
  ) %>% 
  mutate(`P_Average Invested Amount_Outgroup` = as.numeric(`P_Average Invested Amount_Outgroup`),
         `P_Average Invested Amount_Ingroup` = as.numeric(`P_Average Invested Amount_Ingroup`),
         `P_SD invested amount_Outgroup` = as.numeric(`P_SD invested amount_Outgroup`),
         `P_SD invested amount_Ingroup` = as.numeric(`P_SD invested amount_Ingroup`),
         `N by condition_Outgroup` = as.numeric(`N by condition_Outgroup`),
         `N by condition_Ingroup` = as.numeric(`N by condition_Ingroup`)
  ) %>%
  escalc(
    measure = "SMDH", # Hedges' g
    m1i = `P_Average Invested Amount_Outgroup`, sd1i = `P_SD invested amount_Outgroup`, n1i = `N by condition_Outgroup`,
    m2i = `P_Average Invested Amount_Ingroup`,  sd2i = `P_SD invested amount_Ingroup`,  n2i = `N by condition_Ingroup`,
    data = .,
  ) %>% 
  group_by(No., Study) %>%
  mutate(
    effect_n = row_number(),
    effectID = paste(No., Study, effect_n, sep = "_")
  ) %>% 
  ungroup() %>%
  filter(!is.na(yi)) # Exclude effect size is NA

H6_data <- rbind(H6_data, H6_Between_Data)

V_H6 <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H6_data)
H6_model_1 <- rma.mv(yi, V_H6,
                   random = ~ 1 | No. / Study / effectID, data = H6_data, method = "ML")
H6_model_2 <- rma.mv(yi, V_H6,
                     random = ~ 1 | No. / effectID, data = H6_data, method = "ML")
H6_model_3 <- rma.mv(yi, V_H6,
                     random = ~ 1 | effectID, data = H6_data, method = "ML")

H6_compare_1_2 <- anova(H6_model_1, H6_model_2)
H6_compare_1_3 <- anova(H6_model_1, H6_model_3)
H6_compare_2_3 <- anova(H6_model_2, H6_model_3)

## creat a table of model comparison
H6_model_compare <- data.frame(
  Model = c("Four-level model", 
            "Three-level model", 
            "Two-level model"),
  df = c(H6_compare_1_2$parms.f, 
         H6_compare_1_2$parms.r, 
         H6_compare_1_3$parms.r),
  BIC = c(H6_compare_1_2$fit.stats.f[4], 
          H6_compare_1_2$fit.stats.r[4], 
          H6_compare_1_3$fit.stats.r[4]), 
  LL = c(H6_compare_1_2$fit.stats.f[1], 
         H6_compare_1_2$fit.stats.r[1], 
         H6_compare_1_3$fit.stats.r[1]), 
  LRT = c("", 
          round(H6_compare_1_2$LRT, 3), 
          round(H6_compare_1_3$LRT, 3)), 
  p = c(NA, 
        H6_compare_1_2$pval, 
        H6_compare_1_3$pval))

MC_H6 <- H6_model_compare %>% 
  nice_table(title = "Model Comparison for H6")

# refit model with REML
H6_model_f <- rma.mv(yi, V_H6,
                     random = ~ 1 | effectID, data = H6_data, method = "REML")
robust(H6_model_f, cluster = No., clubSandwich=TRUE)

# calculating the overall amount of heterogeneity
sum(H6_model_f$sigma2)

# obtaining I2 value
W <- diag(1/H6_data$vi)
X <- model.matrix(H6_model_f)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
H6_I_t <- 100 * sum(H6_model_f$sigma2) / (sum(H6_model_f$sigma2) + (H6_model_f$k-H6_model_f$p)/sum(diag(P))) # overall I2

## Write table ----
# summary of model comparison
# Save table in folder tab
if (!dir.exists("tab")) {
  dir.create("tab")
}

# summary of best model across H1-H6
Best_table <- data.frame(
  Model = c("H1", 
            "H2", 
            "H5",
            "H6"),
  tau2 = c(sum(H1_model_f$sigma2),
           sum(H2_model_f$sigma2),
           sum(H5_model_f$sigma2),
           sum(H6_model_f$sigma2)),
  `QE(df)` = c(
    paste0(round(H1_model_f$QE,2),"(",H1_model_f$QEdf,")***"),
    paste0(round(H2_model_f$QE,2),"(",H2_model_f$QEdf,")"),
    paste0(round(H5_model_f$QE,2),"(",H5_model_f$QEdf,")*"),
    paste0(round(H6_model_f$QE,2),"(",H6_model_f$QEdf,")")
  ),
  I2_total = c(H1_I_t,
               H2_I_t,
               H5_I_t,
               H6_I_t),
  I2_paper = c(round(H1_I[1],3),
               "",
               round(H5_I[1],3),
               ""),
  I2_study = c(round(H1_I[2],3),
               "","",""),
  I2_in_stu= c(H1_I[3],
               H2_I_t,
               H5_I[2],
               H6_I_t),
  g = c(H1_model_f$b,
        H2_model_f$b,
        H5_model_f$b,
        H6_model_f$b),
  t = c(coef_test(H1_model_f, vcov = "CR2")$tstat,
        coef_test(H2_model_f, vcov = "CR2")$tstat,
        coef_test(H5_model_f, vcov = "CR2")$tstat,
        coef_test(H6_model_f, vcov = "CR2")$tstat),
  df = c(coef_test(H1_model_f, vcov = "CR2")$df_Satt,
         coef_test(H2_model_f, vcov = "CR2")$df_Satt,
         coef_test(H5_model_f, vcov = "CR2")$df_Satt,
         coef_test(H6_model_f, vcov = "CR2")$df_Satt),
  CI_lower = c(conf_int(H1_model_f, vcov = "CR2")$CI_L,
               conf_int(H2_model_f, vcov = "CR2")$CI_L,
               conf_int(H5_model_f, vcov = "CR2")$CI_L,
               conf_int(H6_model_f, vcov = "CR2")$CI_L),
  CI_upper = c(conf_int(H1_model_f, vcov = "CR2")$CI_U,
               conf_int(H2_model_f, vcov = "CR2")$CI_U,
               conf_int(H5_model_f, vcov = "CR2")$CI_U,
               conf_int(H6_model_f, vcov = "CR2")$CI_U),
  p = c(coef_test(H1_model_f, vcov = "CR2")$p_Satt,
        coef_test(H2_model_f, vcov = "CR2")$p_Satt,
        coef_test(H5_model_f, vcov = "CR2")$p_Satt,
        coef_test(H6_model_f, vcov = "CR2")$p_Satt),
  check.names=FALSE
  )

BM <- Best_table %>% 
  nice_table(title = "Summary of Best-Fitting Models for Each Hypothesis")

flextable::save_as_docx(`Table 1` = MC_H1, 
                        `Table 2` = MC_H2, 
                        `Table 3` = MC_H5, 
                        `Table 4` = MC_H6,
                        `Table 5` = BM,
                        path = "tab/S_model_comparison.docx")