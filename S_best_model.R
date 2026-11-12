library(dplyr)
library(tidyr)
library(readr)
library(metafor)
library(clubSandwich)

# reading the data
MA_dat <- read_csv("Coding Sample - Final Meta Data.csv", skip = 1)
MA_dat$`Action_A_Binary` <- as.numeric(MA_dat$`Action_A_Binary`)

# Unfair Condition Analysis ----
## H1 The Effects of Violator’s Identity on Punishment ----
# fair_condition <- c("Fair", "Cooperate", "Neutral", "Generous")
# unfair_condition <- c("Unfair", "Defect", "Selfish", "Cheat")
H1_data <- MA_dat %>%
  filter(
    `Action_A_Binary` == 0, # unfair condition
    `Membership_A` %in% c("Ingroup", "Outgroup")
  ) %>%
  select(
    No., Study, `Membership_A`, `Membership_B`, `Sample_ID`, `IV 1 Levels`, `Age_Mean`,
    `Action_A`, `Action_B`,
    `Punishment_Mean`, `Punishment_SD`, `Sample_Size`, `Total_Trials`, `Intext_Citation`,
    Country, `Cost_Punishment_Ratio`, `Study_Setting`, `Group_Type`, Female, Male, `Game_Type`, `Relative_Endowment`, `Real_Interaction`, `Student_Sample`
  ) %>%
  filter(!is.na(`Punishment_Mean`)) %>% # Exclude mean is NA
  pivot_wider(
    names_from = `Membership_A`,
    values_from = c(`Punishment_Mean`, `Punishment_SD`, `Sample_Size`,`Total_Trials`)
  ) %>% 
  mutate(`Punishment_Mean_Outgroup` = as.numeric(`Punishment_Mean_Outgroup`),
         `Punishment_Mean_Ingroup` = as.numeric(`Punishment_Mean_Ingroup`),
         `Punishment_SD_Outgroup` = as.numeric(`Punishment_SD_Outgroup`),
         `Punishment_SD_Ingroup` = as.numeric(`Punishment_SD_Ingroup`),
         `Sample_Size_Outgroup` = as.numeric(`Sample_Size_Outgroup`),
         `Sample_Size_Ingroup` = as.numeric(`Sample_Size_Ingroup`)
  ) %>%
  escalc(
    measure = "SMDH", # Hedges' g
    m1i = `Punishment_Mean_Outgroup`, sd1i = `Punishment_SD_Outgroup`, n1i = `Sample_Size_Outgroup`,
    m2i = `Punishment_Mean_Ingroup`,  sd2i = `Punishment_SD_Ingroup`,  n2i = `Sample_Size_Ingroup`,
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
    `Group_Manipulation_A` == "Between",
    `Action_A_Binary` == 0, # unfair condition
    `Membership_A` %in% c("Ingroup", "Outgroup")
  ) %>%
  select(
    No., Study, `Membership_A`, `Membership_B`, `IV 1 Levels`, `Age_Mean`,,
    `Action_A`, `Action_B`,
    `Punishment_Mean`, `Punishment_SD`, `Sample_Size`, `Total_Trials`, `Intext_Citation`,
    Country, `Cost_Punishment_Ratio`, `Study_Setting`, `Group_Type`, Female, Male, `Game_Type`, `Relative_Endowment`, `Real_Interaction`, `Student_Sample`
  ) %>%
  filter(!is.na(`Punishment_Mean`)) %>% # Exclude mean is NA
  group_by(No., Study,`Membership_A`) %>% 
  mutate(
    `Sample_ID` = paste(No., Study, `Membership_B`, sep = "_")
  ) %>% 
  ungroup() %>% 
  group_by(`Sample_ID`) %>% 
  mutate(`Age_Mean` = as.numeric(`Age_Mean`),
         `Age_Mean` = weighted.mean(`Age_Mean`, w = `Sample_Size`),
         Female = mean(Female),
         Male = mean(Male)) %>% # have to put the same value for effect sizes calculation
  ungroup() %>% 
  pivot_wider(
    names_from = `Membership_A`,
    values_from = c(`Punishment_Mean`, `Punishment_SD`, `Sample_Size`,`Total_Trials`)
  ) %>% 
  mutate(`Punishment_Mean_Outgroup` = as.numeric(`Punishment_Mean_Outgroup`),
         `Punishment_Mean_Ingroup` = as.numeric(`Punishment_Mean_Ingroup`),
         `Punishment_SD_Outgroup` = as.numeric(`Punishment_SD_Outgroup`),
         `Punishment_SD_Ingroup` = as.numeric(`Punishment_SD_Ingroup`),
         `Sample_Size_Outgroup` = as.numeric(`Sample_Size_Outgroup`),
         `Sample_Size_Ingroup` = as.numeric(`Sample_Size_Ingroup`)
  ) %>%
  escalc(
    measure = "SMDH", # Hedges' g
    m1i = `Punishment_Mean_Outgroup`, sd1i = `Punishment_SD_Outgroup`, n1i = `Sample_Size_Outgroup`,
    m2i = `Punishment_Mean_Ingroup`,  sd2i = `Punishment_SD_Ingroup`,  n2i = `Sample_Size_Ingroup`,
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

V_H1 <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H1_data)
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

# refit model with REML
H1_model_f <- rma.mv(yi, V_H1, random = ~ 1 | No. / Study / effectID, data = H1_data, method = "REML") # use REML for better fitting
# obtaining I2 value
W <- diag(1/H1_data$vi)
X <- model.matrix(H1_model_f)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
H1_I_t <- 100 * sum(H1_model_f$sigma2) / (sum(H1_model_f$sigma2) + (H1_model_f$k-H1_model_f$p)/sum(diag(P))) # overall I2
H1_I <- 100 * H1_model_f$sigma2 / (sum(H1_model_f$sigma2) + (H1_model_f$k-H1_model_f$p)/sum(diag(P)))  # I2 per levels

## H2 (Deterrence perspective): ----
H2_data <- MA_dat %>%
  filter(
    `Action_A_Binary` == 0, # unfair condition
    `Membership_B` %in% c("Ingroup", "Outgroup")
  ) %>%
  # group_by(No., Study, `Membership_B`) %>%
  select(
    No., Study, `Membership_A`, `Membership_B`, `Sample_ID`, `IV 1 Levels`, `Age_Mean`,,
    `Action_A`, `Action_B`,
    `Punishment_Mean`, `Punishment_SD`, `Sample_Size`, `Total_Trials`, `Intext_Citation`,
    Country, `Cost_Punishment_Ratio`, `Study_Setting`, `Group_Type`, Female, Male, `Game_Type`, `Relative_Endowment`, `Real_Interaction`, `Student_Sample`
  ) %>%
  filter(!is.na(`Punishment_Mean`)) %>% # Exclude mean is NA
  pivot_wider(
    names_from = `Membership_B`,
    values_from = c(`Punishment_Mean`, `Punishment_SD`, `Sample_Size`,`Total_Trials`)
  ) %>% 
  mutate(`Punishment_Mean_Outgroup` = as.numeric(`Punishment_Mean_Outgroup`),
         `Punishment_Mean_Ingroup` = as.numeric(`Punishment_Mean_Ingroup`),
         `Punishment_SD_Outgroup` = as.numeric(`Punishment_SD_Outgroup`),
         `Punishment_SD_Ingroup` = as.numeric(`Punishment_SD_Ingroup`),
         `Sample_Size_Outgroup` = as.numeric(`Sample_Size_Outgroup`),
         `Sample_Size_Ingroup` = as.numeric(`Sample_Size_Ingroup`)
  ) %>% 
  escalc(
    measure = "SMDH", # Hedges' g
    m1i = `Punishment_Mean_Outgroup`, sd1i = `Punishment_SD_Outgroup`, n1i = `Sample_Size_Outgroup`,
    m2i = `Punishment_Mean_Ingroup`,  sd2i = `Punishment_SD_Ingroup`,  n2i = `Sample_Size_Ingroup`,
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
    `Group_Manipulation_B` == "Between",
    `Action_A_Binary` == 0, # unfair condition
    `Membership_B` %in% c("Ingroup", "Outgroup")
  ) %>%
  select(
    No., Study, `Membership_A`, `Membership_B`, `IV 1 Levels`, `Age_Mean`,,
    `Action_A`, `Action_B`,
    `Punishment_Mean`, `Punishment_SD`, `Sample_Size`, `Total_Trials`, `Intext_Citation`,
    Country, `Cost_Punishment_Ratio`, `Study_Setting`, `Group_Type`, Female, Male, `Game_Type`, `Relative_Endowment`, `Real_Interaction`, `Student_Sample`
  ) %>%
  filter(!is.na(`Punishment_Mean`)) %>% # Exclude mean is NA
  group_by(No., Study,`Membership_B`) %>% 
  mutate(
    `Sample_ID` = paste(No., Study, `Membership_A`, sep = "_")
  ) %>% 
  ungroup() %>% 
  group_by(`Sample_ID`) %>% 
  mutate(`Age_Mean` = as.numeric(`Age_Mean`),
         `Age_Mean` = weighted.mean(`Age_Mean`, w = `Sample_Size`),
         Female = mean(Female),
         Male = mean(Male)) %>% # have to put the same value for effect sizes calculation
  ungroup() %>% 
  pivot_wider(
    names_from = `Membership_B`,
    values_from = c(`Punishment_Mean`, `Punishment_SD`, `Sample_Size`,`Total_Trials`)
  ) %>% 
  mutate(`Punishment_Mean_Outgroup` = as.numeric(`Punishment_Mean_Outgroup`),
         `Punishment_Mean_Ingroup` = as.numeric(`Punishment_Mean_Ingroup`),
         `Punishment_SD_Outgroup` = as.numeric(`Punishment_SD_Outgroup`),
         `Punishment_SD_Ingroup` = as.numeric(`Punishment_SD_Ingroup`),
         `Sample_Size_Outgroup` = as.numeric(`Sample_Size_Outgroup`),
         `Sample_Size_Ingroup` = as.numeric(`Sample_Size_Ingroup`)
  ) %>%
  escalc(
    measure = "SMDH", # Hedges' g
    m1i = `Punishment_Mean_Outgroup`, sd1i = `Punishment_SD_Outgroup`, n1i = `Sample_Size_Outgroup`,
    m2i = `Punishment_Mean_Ingroup`,  sd2i = `Punishment_SD_Ingroup`,  n2i = `Sample_Size_Ingroup`,
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

V_H2 <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H2_data)
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

# refit model with REML
H2_model_f <- rma.mv(yi, V_H2, random = ~ 1 | effectID, data = H2_data, method = "REML")

# obtaining I2 value
W <- diag(1/H2_data$vi)
X <- model.matrix(H2_model_f)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
H2_I_t <- 100 * sum(H2_model_f$sigma2) / (sum(H2_model_f$sigma2) + (H2_model_f$k-H2_model_f$p)/sum(diag(P))) # overall I2

# Fair Condition Analysis ################
## H5 The Effects of Violator’s Identity on Punishment ----
H5_data <- MA_dat %>%
  filter(
    `Action_A_Binary` == 1, # fair condition
    `Membership_A` %in% c("Ingroup", "Outgroup")
  ) %>%
  select(
    No., Study, `Membership_A`, `Membership_B`, `Sample_ID`, `IV 1 Levels`, `Age_Mean`,,
    `Action_A`, `Action_B`,
    `Punishment_Mean`, `Punishment_SD`, `Sample_Size`, `Total_Trials`, `Intext_Citation`,
    Country, `Cost_Punishment_Ratio`, `Study_Setting`, `Group_Type`, Female, Male, `Game_Type`, `Relative_Endowment`, `Real_Interaction`, `Student_Sample`
  ) %>%
  filter(!is.na(`Punishment_Mean`)) %>% # Exclude mean is NA
  pivot_wider(
    names_from = `Membership_A`,
    values_from = c(`Punishment_Mean`, `Punishment_SD`, `Sample_Size`,`Total_Trials`)
  ) %>%
  mutate(`Punishment_Mean_Outgroup` = as.numeric(`Punishment_Mean_Outgroup`),
         `Punishment_Mean_Ingroup` = as.numeric(`Punishment_Mean_Ingroup`),
         `Punishment_SD_Outgroup` = as.numeric(`Punishment_SD_Outgroup`),
         `Punishment_SD_Ingroup` = as.numeric(`Punishment_SD_Ingroup`),
         `Sample_Size_Outgroup` = as.numeric(`Sample_Size_Outgroup`),
         `Sample_Size_Ingroup` = as.numeric(`Sample_Size_Ingroup`)
  ) %>%
  escalc(
    measure = "SMDH", # Hedges' g
    m1i = `Punishment_Mean_Outgroup`, sd1i = `Punishment_SD_Outgroup`, n1i = `Sample_Size_Outgroup`,
    m2i = `Punishment_Mean_Ingroup`,  sd2i = `Punishment_SD_Ingroup`,  n2i = `Sample_Size_Ingroup`,
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
    `Group_Manipulation_A` == "Between",
    `Action_A_Binary` == 1, # fair condition
    `Membership_A` %in% c("Ingroup", "Outgroup")
  ) %>%
  select(
    No., Study, `Membership_A`, `Membership_B`, `IV 1 Levels`, `Age_Mean`,,
    `Action_A`, `Action_B`,
    `Punishment_Mean`, `Punishment_SD`, `Sample_Size`, `Total_Trials`, `Intext_Citation`,
    Country, `Cost_Punishment_Ratio`, `Study_Setting`, `Group_Type`, Female, Male, `Game_Type`, `Relative_Endowment`, `Real_Interaction`, `Student_Sample`
  ) %>%
  filter(!is.na(`Punishment_Mean`)) %>% # Exclude mean is NA
  group_by(No., Study,`Membership_A`) %>% 
  mutate(
    `Sample_ID` = paste(No., Study, `Membership_B`, sep = "_")
  ) %>% 
  ungroup() %>% 
  group_by(`Sample_ID`) %>% 
  mutate(`Age_Mean` = as.numeric(`Age_Mean`),
         `Age_Mean` = weighted.mean(`Age_Mean`, w = `Sample_Size`),
         Female = mean(Female),
         Male = mean(Male)) %>% # have to put the same value for effect sizes calculation
  ungroup() %>% 
  pivot_wider(
    names_from = `Membership_A`,
    values_from = c(`Punishment_Mean`, `Punishment_SD`, `Sample_Size`,`Total_Trials`)
  ) %>% 
  mutate(`Punishment_Mean_Outgroup` = as.numeric(`Punishment_Mean_Outgroup`),
         `Punishment_Mean_Ingroup` = as.numeric(`Punishment_Mean_Ingroup`),
         `Punishment_SD_Outgroup` = as.numeric(`Punishment_SD_Outgroup`),
         `Punishment_SD_Ingroup` = as.numeric(`Punishment_SD_Ingroup`),
         `Sample_Size_Outgroup` = as.numeric(`Sample_Size_Outgroup`),
         `Sample_Size_Ingroup` = as.numeric(`Sample_Size_Ingroup`)
  ) %>%
  escalc(
    measure = "SMDH", # Hedges' g
    m1i = `Punishment_Mean_Outgroup`, sd1i = `Punishment_SD_Outgroup`, n1i = `Sample_Size_Outgroup`,
    m2i = `Punishment_Mean_Ingroup`,  sd2i = `Punishment_SD_Ingroup`,  n2i = `Sample_Size_Ingroup`,
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

V_H5 <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H5_data)
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

# refit model with REML
H5_model_f <- rma.mv(yi, V_H5, random = ~ 1 | No. / effectID, data = H5_data, method = "REML")

# obtaining I2 value
W <- diag(1/H5_data$vi)
X <- model.matrix(H5_model_f)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
H5_I_t <- 100 * sum(H5_model_f$sigma2) / (sum(H5_model_f$sigma2) + (H5_model_f$k-H5_model_f$p)/sum(diag(P))) # overall I2
H5_I <- 100 * H5_model_f$sigma2 / (sum(H5_model_f$sigma2) + (H5_model_f$k-H5_model_f$p)/sum(diag(P)))  # I2 per levels

## H6 The Effects of Victim’s Identity on Punishment ----
H6_data <- MA_dat %>%
  filter(
    `Action_A_Binary` == 1, # fair condition
    `Membership_B` %in% c("Ingroup", "Outgroup")
  ) %>%
  # group_by(No., Study, `Membership_B`) %>%
  select(
    No., Study, `Membership_A`, `Membership_B`, `Sample_ID`, `IV 1 Levels`, `Age_Mean`,,
    `Action_A`, `Action_B`,
    `Punishment_Mean`, `Punishment_SD`, `Sample_Size`, `Total_Trials`, `Intext_Citation`,
    Country, `Cost_Punishment_Ratio`, `Study_Setting`, `Group_Type`, Female, Male, `Game_Type`, `Relative_Endowment`, `Real_Interaction`, `Student_Sample`
  ) %>%
  filter(!is.na(`Punishment_Mean`)) %>% # Exclude mean is NA
  pivot_wider(
    names_from = `Membership_B`,
    values_from = c(`Punishment_Mean`, `Punishment_SD`, `Sample_Size`,`Total_Trials`)
  ) %>% 
  mutate(`Punishment_Mean_Outgroup` = as.numeric(`Punishment_Mean_Outgroup`),
         `Punishment_Mean_Ingroup` = as.numeric(`Punishment_Mean_Ingroup`),
         `Punishment_SD_Outgroup` = as.numeric(`Punishment_SD_Outgroup`),
         `Punishment_SD_Ingroup` = as.numeric(`Punishment_SD_Ingroup`),
         `Sample_Size_Outgroup` = as.numeric(`Sample_Size_Outgroup`),
         `Sample_Size_Ingroup` = as.numeric(`Sample_Size_Ingroup`)
  ) %>%
  escalc(
    measure = "SMDH", # Hedges' g
    m1i = `Punishment_Mean_Outgroup`, sd1i = `Punishment_SD_Outgroup`, n1i = `Sample_Size_Outgroup`,
    m2i = `Punishment_Mean_Ingroup`,  sd2i = `Punishment_SD_Ingroup`,  n2i = `Sample_Size_Ingroup`,
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
    `Group_Manipulation_B` == "Between",
    `Action_A_Binary` == 1, # fair condition
    `Membership_B` %in% c("Ingroup", "Outgroup")
  ) %>%
  select(
    No., Study, `Membership_A`, `Membership_B`, `IV 1 Levels`, `Age_Mean`,,
    `Action_A`, `Action_B`,
    `Punishment_Mean`, `Punishment_SD`, `Sample_Size`, `Total_Trials`, `Intext_Citation`,
    Country, `Cost_Punishment_Ratio`, `Study_Setting`, `Group_Type`, Female, Male, `Game_Type`, `Relative_Endowment`, `Real_Interaction`, `Student_Sample`
  ) %>%
  filter(!is.na(`Punishment_Mean`)) %>% # Exclude mean is NA
  group_by(No., Study,`Membership_B`) %>% 
  mutate(
    `Sample_ID` = paste(No., Study, `Membership_A`, sep = "_")
  ) %>% 
  ungroup() %>%
  group_by(`Sample_ID`) %>% 
  mutate(`Age_Mean` = as.numeric(`Age_Mean`),
         `Age_Mean` = weighted.mean(`Age_Mean`, w = `Sample_Size`),
         Female = mean(Female),
         Male = mean(Male)) %>% # have to put the same value for effect sizes calculation
  ungroup() %>% 
  pivot_wider(
    names_from = `Membership_B`,
    values_from = c(`Punishment_Mean`, `Punishment_SD`, `Sample_Size`,`Total_Trials`)
  ) %>% 
  mutate(`Punishment_Mean_Outgroup` = as.numeric(`Punishment_Mean_Outgroup`),
         `Punishment_Mean_Ingroup` = as.numeric(`Punishment_Mean_Ingroup`),
         `Punishment_SD_Outgroup` = as.numeric(`Punishment_SD_Outgroup`),
         `Punishment_SD_Ingroup` = as.numeric(`Punishment_SD_Ingroup`),
         `Sample_Size_Outgroup` = as.numeric(`Sample_Size_Outgroup`),
         `Sample_Size_Ingroup` = as.numeric(`Sample_Size_Ingroup`)
  ) %>%
  escalc(
    measure = "SMDH", # Hedges' g
    m1i = `Punishment_Mean_Outgroup`, sd1i = `Punishment_SD_Outgroup`, n1i = `Sample_Size_Outgroup`,
    m2i = `Punishment_Mean_Ingroup`,  sd2i = `Punishment_SD_Ingroup`,  n2i = `Sample_Size_Ingroup`,
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

V_H6 <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H6_data)
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

# refit model with REML
H6_model_f <- rma.mv(yi, V_H6,
                     random = ~ 1 | effectID, data = H6_data, method = "REML")
# obtaining I2 value
W <- diag(1/H6_data$vi)
X <- model.matrix(H6_model_f)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
H6_I_t <- 100 * sum(H6_model_f$sigma2) / (sum(H6_model_f$sigma2) + (H6_model_f$k-H6_model_f$p)/sum(diag(P))) # overall I2

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
  `QE(df)` = vapply(list(H1_model_f, H2_model_f, H5_model_f, H6_model_f), function(model) {
    stars <- symnum(model$QEp, cutpoints = c(0, .001, .01, .05, 1),
                    symbols = c("***", "**", "*", ""))
    paste0(round(model$QE, 2), "(", model$k - model$p, ")", stars)
  }, character(1)),
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

save(H1_model_compare, H2_model_compare, H5_model_compare, H6_model_compare,
     Best_table, file = "best_model.rds", compress = "xz")
