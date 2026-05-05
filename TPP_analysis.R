# Notes:       Requires 'Coding Sample - Final Meta Data.csv' in the folder.


# if rma.mv failed to converge, add control=list(rel.tol=1e-9) in the function
required_packages <- c("tidyverse", "metafor", "clubSandwich", "MAd", "weightr", "cowplot", "metapower", "rempsyc")

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

# read the data
MA_dat <- read_csv("Coding Sample - Final Meta Data.csv", skip = 1)
MA_dat$`Action_A_Binary` <- as.numeric(MA_dat$`Action_A_Binary`)

# Main Analysis ================================================================
## H1 The Effect of Focal Player Membership on Altruistic Punishment --------------
### Data Preparation ----
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

# deal with between data (assign temporary Sample_ID)
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

### Model Fitting ----
H1_model <- rma.mv(yi, V_H1, random = ~ 1 | No. / Study / effectID, data = H1_data, method = "REML")
robust(H1_model, cluster = No., clubSandwich=TRUE) %>% predict()
robust(H1_model, cluster = No., clubSandwich=TRUE) %>% confint()

# checking likelihood profile plots
profile(H1_model)

# calculating the overall amount of heterogeneity
sum(H1_model$sigma2)

# obtaining I2 value
W <- diag(1/H1_data$vi)
X <- model.matrix(H1_model)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(H1_model$sigma2) / (sum(H1_model$sigma2) + (H1_model$k-H1_model$p)/sum(diag(P))) # overall I2
100 * H1_model$sigma2 / (sum(H1_model$sigma2) + (H1_model$k-H1_model$p)/sum(diag(P)))  # I2 per levels

# forest plot
forest(H1_model, slab = paste(H1_data$Intext_Citation, "S", H1_data$Study, sep = " "))
# Save fig in folder fig
if (!dir.exists("fig")) {
  dir.create("fig")
}
# saving forest plot to a file
jpeg("fig/H1_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H1_model, slab = paste(H1_data$Intext_Citation, "S", H1_data$Study, sep = " "))
dev.off()

# funnel plot
metafor::funnel(H1_model)

# contour-enhanced funnel plot
jpeg("fig/H1_funnel.jpeg", height = 6, width = 6, units = "in", res = 300)
metafor::funnel(H1_model, refline=0, level=c(95, 99), 
                shade=c("white", "grey"))
dev.off()


## H2 The Effect of Partner Player Membership on Altruistic Punishment ----
### Data Preparation ----
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

# deal with between data (assign temporary Sample_ID)
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

### Model Fitting ----
H2_model <- rma(yi, vi, random = ~ 1 | effectID, data = H2_data, method = "REML", test = "knha")
summary(H2_model) # don't use CR2 as sigma is 0 
H2_model %>% predict()
H2_model %>% confint()

# checking likelihood profile plots
profile(H2_model)

# forest plot
H2_f <- forest(H2_model, slab = paste(H2_data$Intext_Citation, "S", H2_data$Study, sep = " "))
# saving forest plot to a file
jpeg("fig/H2_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H2_model, slab = paste(H2_data$Intext_Citation, "S", H2_data$Study, sep = " "))
dev.off()

# funnel plot
metafor::funnel(H2_model)

# contour-enhanced funnel plot
jpeg("fig/H2_funnel.jpeg", height = 6, width = 6, units = "in", res = 300)
metafor::funnel(H2_model, refline=0, level=c(95, 99), 
                shade=c("white", "grey"))
dev.off()

## Interaction Effect of Group Membership on Altruistic Punishment----
# We applied two methods, moderator method and regression method
### The moderator method ----
H3_model <- H1_data %>%
  rma.mv(yi, V_H1,
         mods = ~ Membership_B - 1,
         random = ~ 1 | No. / Study / effectID,
         data = .,
         method = "REML",
  )
robust(H3_model, cluster = No., clubSandwich=TRUE)

Wald_test(H3_model, constrain_equal(1:3), vcov = "CR2")

# asign color for ingroup, outgroup, and neutral
H3_col <- ifelse(H1_data$Membership_B == "Ingroup", "#004488",
                 ifelse(H1_data$Membership_B == "Outgroup", "#BB5566", "gray"))
jpeg("fig/H3_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H3_model, slab = paste(H1_data$Intext_Citation, "S", H1_data$Study, sep = " "),
       colout = H3_col)
dev.off()

#### Robustness check 1: Exclude studies with neutral recipient/player B affiliation ----
H3_data <- H1_data %>%
  filter(Membership_B != "Neutral")

V_H3 <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H3_data)

H3_model2 <-H3_data  %>%
  rma.mv(yi, V_H3,
         mods = ~ Membership_B - 1,
         random = ~ 1 | No. / Study / effectID,
         data = .,
         method = "REML",
  )
robust(H3_model2, cluster = No., clubSandwich=TRUE)

Wald_test(H3_model2, constraints = constrain_equal(1:2), vcov = "CR2")

#### Robustness check 2: Using the H2_data ----
H3_model3 <- H2_data %>%
  rma.mv(yi, vi, # V cannot convergence
         mods = ~ Membership_A - 1,
         random = ~ 1 | No. / Study / effectID,
         data = .,
         method = "REML",
  )
robust(H3_model3, cluster = No., clubSandwich=TRUE)

Wald_test(H3_model3, constrain_equal(1:3), vcov = "CR2")

#### Robustness check 3: Using H2_data &Exclude studies with neutral recipient/player B affiliation ----
H3_data2 <- H2_data %>%
  filter(Membership_A != "Neutral")

V_H3_2 <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H3_data2)

H3_model4 <-H3_data2  %>%
  rma.mv(yi, V_H3_2,
         mods = ~ Membership_A - 1,
         random = ~ 1 | No. / Study / effectID,
         data = .,
         method = "REML",
  )
robust(H3_model4, cluster = No., clubSandwich=TRUE)

Wald_test(H3_model4, constrain_equal(1:2), vcov = "CR2")

### The regression method ----
#### Data Preparation ----
H3_data_r <- MA_dat %>%
  filter(
    `Action_A_Binary` == 0, # unfair condition
    `Membership_A` %in% c("Ingroup", "Outgroup"),
    `Membership_B` %in% c("Ingroup", "Outgroup"),
    !is.na(`Punishment_Mean_STD`),
    !is.na(`Punishment_SD_STD`),
    !is.na(`Total_Trials`)
  ) %>%
  mutate(
    `Total_Trials` = as.numeric(`Total_Trials`),
    `Person A` = factor(`Membership_A`, levels = c("Ingroup", "Outgroup")),
    `Person B` = factor(`Membership_B`, levels = c("Ingroup", "Outgroup")),
    P_SV = as.numeric(`Punishment_SD_STD`)^2 * `Total_Trials` / ( `Total_Trials` - 1), # Sampling variance of punishment amount
    effect_n = row_number(),
    effectID = paste(No., Study, effect_n, sep = "_"),
    P_log = car::logit(`Punishment_Mean_STD`),
    SV_log = `Punishment_SD`^2 / (`Punishment_Mean`^2 * `Total_Trials` * (1-`Punishment_Mean_STD`)^2)
  )

V_H3  <- vcalc(SV_log, cluster = `Sample_ID`, obs = effectID, rho = 0.5, data = H3_data_r)

H4_data_r <- MA_dat %>%
  filter(
    `Action_A_Binary` == 1, # fair condition
    `Membership_A` %in% c("Ingroup", "Outgroup"),
    `Membership_B` %in% c("Ingroup", "Outgroup"),
    !is.na(`Punishment_Mean_STD`),
    !is.na(`Punishment_SD_STD`),
    !is.na(`Total_Trials`)
  ) %>%
  mutate(
    `Total_Trials` = as.numeric(`Total_Trials`),
    `Person A` = factor(`Membership_A`, levels = c("Ingroup", "Outgroup")),
    `Person B` = factor(`Membership_B`, levels = c("Ingroup", "Outgroup")),
    P_SV = as.numeric(`Punishment_SD_STD`)^2 * `Total_Trials` / ( `Total_Trials` - 1), # Sampling variance of punishment amount
    effect_n = row_number(),
    effectID = paste(No., Study, effect_n, sep = "_"),
    P_log = car::logit(`Punishment_Mean_STD`),
    SV_log = `Punishment_SD`^2 / (`Punishment_Mean`^2 * `Total_Trials` * (1-`Punishment_Mean_STD`)^2)
  )

V_H4  <- vcalc(SV_log, cluster = `Sample_ID`, obs = effectID, rho = 0.5, data = H4_data_r)

#### Model Fitting ----
H3_model_r <- H3_data_r %>% # unfair
  rma.mv(yi = P_log, V = V_H3,
         mods = ~ `Person A` * `Person B`,
         random = ~ 1 | No. / Study / effectID,  data = .,
         method = "REML",
  )
robust(H3_model_r, cluster = No., clubSandwich=TRUE)

H4_model_r <- H4_data_r %>% # fair
  rma.mv(yi = P_log, V = V_H4,
         mods = ~ `Person A` * `Person B`,
         random = ~ 1 | No. / Study / effectID,  data = .,
         method = "REML",
  )
robust(H4_model_r, cluster = No., clubSandwich=TRUE)

#### Estimated Average Effects ----
H3_data_d <- MA_dat %>%
  filter(
    `Action_A_Binary` == 0, # unfair condition
    !is.na(`Punishment_Mean_STD`),
    !is.na(`Punishment_SD_STD`),
    !is.na(`Total_Trials`)
  ) %>%
  mutate(
    `Total_Trials` = as.numeric(`Total_Trials`),
    `Person A` = factor(`Membership_A`, levels = c("Ingroup", "Neutral", "Outgroup")),
    `Person B` = factor(`Membership_B`, levels = c("Ingroup", "Neutral", "Outgroup")),
    P_SV = as.numeric(`Punishment_SD_STD`)^2 * `Total_Trials` / ( `Total_Trials` - 1), # Sampling variance of punishment amount
    effect_n = row_number(),
    effectID = paste(No., Study, effect_n, sep = "_"),
    P_log = car::logit(`Punishment_Mean_STD`),
    SV_log = `Punishment_SD`^2 / (`Punishment_Mean`^2 * `Total_Trials` * (1-`Punishment_Mean_STD`)^2)
  )

V_H3_d  <- vcalc(SV_log, cluster = `Sample_ID`, obs = effectID, rho = 0.5, data = H3_data_d)
H3_model_d <- H3_data_d %>%
  rma.mv(yi = P_log, V = V_H3_d,
         mods = ~ `Person A` * `Person B`,
         random = ~ 1 | No. / Study / effectID,  data = .,
         method = "REML",
  )
robust(H3_model_d, cluster = No., clubSandwich=TRUE)

H3_mean_M <- H3_model_d %>% 
  model.matrix() %>% 
  colMeans(.)[-1]
H3_overall <- predict(H3_model_d, newmods = H3_mean_M) %>% 
  as.data.frame() %>% 
  lapply(boot::inv.logit) %>% as.data.frame()
H3_overall

results_mat_name <- expand.grid(`Person A` = c("Ingroup", "Neutral", "Outgroup"), `Person B` = c("Ingroup", "Neutral", "Outgroup"))
results_mat <- model.matrix(~ `Person A` * `Person B`, results_mat_name)[,-1]
results_dat <- predict(H3_model_d, newmods = results_mat) %>% 
  cbind(results_mat_name) %>% 
  mutate(pred = boot::inv.logit(pred))

heat_all <- results_dat %>% 
  ggplot(aes(x = `Person A`, y = `Person B`, fill = pred)) +
  geom_tile() +
  geom_text(aes(label = paste0(format(round(pred * 100,2), nsmall = 2), "%"))) +
  scale_fill_gradient(low = "white", high = "#BB5566",
                      limits = c(0, 0.5)) +
  theme_cowplot() +
  labs(fill = "Punishment", x = "Focal Player Group Membership",
       y = "Partner Player Group Membership")
ggsave("fig/Summary_unfair_heatmap.jpeg", heat_all, height = 6, width = 7)

H4_data_d <- MA_dat %>% # not appropriate, but let's say it is H4 temporary
  filter(
    `Action_A_Binary` == 1, # fair condition
    !is.na(`Punishment_Mean_STD`),
    !is.na(`Punishment_SD_STD`),
    !is.na(`Total_Trials`)
  ) %>%
  mutate(
    `Total_Trials` = as.numeric(`Total_Trials`),
    `Person A` = factor(`Membership_A`, levels = c("Ingroup", "Neutral", "Outgroup")),
    `Person B` = factor(`Membership_B`, levels = c("Ingroup", "Neutral", "Outgroup")),
    P_SV = as.numeric(`Punishment_SD_STD`)^2 * `Total_Trials` / ( `Total_Trials` - 1), # Sampling variance of punishment amount
    effect_n = row_number(),
    effectID = paste(No., Study, effect_n, sep = "_"),
    P_log = car::logit(`Punishment_Mean_STD`),
    SV_log = `Punishment_SD`^2 / (`Punishment_Mean`^2 * `Total_Trials` * (1-`Punishment_Mean_STD`)^2),
    SV_log = ifelse(is.na(SV_log), 1/`Total_Trials`, SV_log)
  )

V_H4_d  <- vcalc(SV_log, cluster = `Sample_ID`, obs = effectID, rho = 0.5, data = H4_data_d)
H4_model_d <- H4_data_d %>%
  rma.mv(yi = P_log, V = V_H4_d, # V close to singular
         mods = ~ `Person A` * `Person B`,
         random = ~ 1 | No. / Study / effectID,  data = .,
         method = "REML",
  )
summary(H4_model_d)

H4_mean_M <- H4_model_d %>% 
  model.matrix() %>% 
  colMeans(.)[-1]
H4_overall <- predict(H4_model_d, newmods = H4_mean_M) %>% 
  as.data.frame() %>% 
  lapply(boot::inv.logit) %>% as.data.frame() %>% 
  mutate(across(everything(), ~ (.x - 0.025) / (1 - 2 * 0.025)))
H4_overall


results_mat2 <- results_mat[-8,-7]
results_dat2 <- predict(H4_model_d, newmods = results_mat2) %>% 
  cbind(results_mat_name[-8,]) %>% 
  right_join(results_mat_name) %>% 
  mutate(pred = boot::inv.logit(pred),
         pred = (pred - 0.025) / (1 - 2 * 0.025) # Note that due to 0 in the data, car::logit() remapped data, so we need invert it as well
  )

heat_all_2 <- results_dat2 %>% 
  ggplot(aes(x = `Person A`, y = `Person B`, fill = pred)) +
  geom_tile() +
  geom_text(aes(label = ifelse(is.na(pred),"NA",
                               paste0(format(round(pred * 100,2), nsmall = 2), "%"))
  )) +
  scale_fill_gradient(low = "white", high = "#BB5566",
                      limits = c(0, 0.5)) +
  theme_cowplot() +
  labs(fill = "Punishment", x = "Focal Player Group Membership",
       y = "Partner Player Group Membership")
ggsave("fig/Summary_fair_heatmap.jpeg", heat_all_2, height = 6, width = 7)

heat_all_all <- plot_grid(heat_all, heat_all_2, labels = c("a", "b"))
ggsave("fig/Summary_all_heatmap.jpeg", heat_all_all, height = 6, width = 14)

## H5 The Effect of Focal Player Membership on Antisocial Punishment --------------
### Data Preparation ----
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

### Model Fitting ----
H5_model <- rma.mv(yi, V_H5, random = ~ 1 | No. / effectID, data = H5_data, method = "REML")
robust(H5_model, cluster = No., clubSandwich=TRUE) %>% predict()
robust(H5_model, cluster = No., clubSandwich=TRUE) %>% confint()

# calculating the overall amount of heterogeneity
sum(H5_model$sigma2)

# obtaining I2 value
W <- diag(1/H5_data$vi)
X <- model.matrix(H5_model)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(H5_model$sigma2) / (sum(H5_model$sigma2) + (H5_model$k-H5_model$p)/sum(diag(P))) # overall I2
100 * H5_model$sigma2 / (sum(H5_model$sigma2) + (H5_model$k-H5_model$p)/sum(diag(P)))  # I2 per levels

# forest plot
H5_f <- forest(H5_model, slab = paste(H5_data$Intext_Citation, "S", H5_data$Study, sep = " "))
# saving forest plot to a file
jpeg("fig/H5_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H5_model, slab = paste(H5_data$Intext_Citation, "S", H5_data$Study, sep = " "))
dev.off()

# funnel plot
metafor::funnel(H5_model)

# contour-enhanced funnel plot
jpeg("fig/H5_funnel.jpeg", height = 6, width = 6, units = "in", res = 300)
metafor::funnel(H5_model, refline=0, level=c(95, 99), 
                shade=c("white", "grey"))
dev.off()

## H6 The Effects of Victim’s Identity on Punishment ----
### Data Preparation ----
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

### Model Fitting ----
H6_model <- rma(yi, vi, random = ~ 1 | effectID, data = H6_data, method = "REML", test = "knha")
summary(H6_model) # don't use CR2 as sigma is 0 
H6_model %>% predict()

# forest plot
H6_f <- forest(H6_model, slab = paste(H6_data$Intext_Citation, "S", H6_data$Study, sep = " "))
# saving forest plot to a file
jpeg("fig/H6_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H6_model, slab = paste(H6_data$Intext_Citation, "S", H6_data$Study, sep = " "))
dev.off()

# funnel plot
metafor::funnel(H6_model)

# contour-enhanced funnel plot
jpeg("fig/H6_funnel.jpeg", height = 6, width = 6, units = "in", res = 300)
metafor::funnel(H6_model, refline=0, level=c(95, 99), 
                shade=c("white", "grey"))
dev.off()

# H4 ----
H4_model <- H5_data %>%
  rma.mv(yi, V_H5,
         mods = ~ Membership_B - 1,
         random = ~ 1 | No. / effectID,
         data = .,
         method = "REML",
  )
robust(H4_model, cluster = No., clubSandwich=TRUE)

Wald_test(H4_model, constrain_equal(1:3), vcov = "CR2")

# asign color for ingroup, outgroup, and neutral
H4_col <- ifelse(H5_data$Membership_B == "Ingroup", "#004488",
                 ifelse(H5_data$Membership_B == "Outgroup", "#BB5566", "gray"))
jpeg("fig/H4_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H4_model, slab = paste(H5_data$Intext_Citation, "S", H5_data$Study, sep = " "),
       colout = H4_col)
dev.off()


# Single moderator analyses ################
## The Effects of Focal Player Membership on Punishment moderator analyses -----
H1_moderator_table <- tibble()
### The culture ----
H1_country_dat <- H1_data %>% 
  filter(!is.na(Country)) %>% 
  filter(Country != "UK/US") # only one, and that paper reported UK & US separately in other study
V_H1_country <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H1_country_dat)
H1_country_model <- rma.mv(yi, V_H1_country,
                           mods = ~ Country - 1,
                           random = ~ 1 | No. / Study / effectID, data = H1_country_dat, method = "REML")
summary(H1_country_model)

Wald_test(H1_country_model, constraints = constrain_equal(1:7), vcov = "CR2", test = "EDT") # HTZ returns NA
coef_test(H1_country_model, vcov = "CR2")
conf_int(H1_country_model, vcov = "CR2")

# saving forest plot to a file
H1_country_col <- case_match(
  H1_country_dat$Country,
  "China"  ~ "#CC6677",
  "Turkey" ~ "#DDCC77",
  "United States"     ~ "#999933",
  "United Kingdom"  ~ "#117733",
  "Argentina" ~ "#44AA99",
  "Switzerland" ~ "#88CCEE",
  "Papua New Guinea" ~ "#332288"
)

jpeg("fig/H1_country_Forest.jpeg", height = 12, width = 10, units = "in", res = 300)
forest(H1_country_model, slab = paste(H1_country_dat$Intext_Citation, "S", H1_country_dat$Study, sep = " "),
       colout = H1_country_col, ilab = Country, ilab.lab = "Country", ilab.xpos = -3.5, xlim = c(-10, 6.5), alim = c(-2, 4))
dev.off()

# saving statistics
H1_moderator_table <- left_join(
  coef_test(H1_country_model, vcov = "CR2"),
  conf_int(H1_country_model, vcov = "CR2")
)

### The cost ----
H1_cost_dat <- H1_data %>% 
  mutate(Cost_Punishment_Ratio = as.numeric(Cost_Punishment_Ratio)) %>% 
  filter(!is.na(Cost_Punishment_Ratio))
V_H1_cost <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H1_cost_dat)
H1_cost_model <- rma.mv(yi, V_H1_cost,
                        mods = ~ Cost_Punishment_Ratio,
                        random = ~ 1 | No. / Study / effectID, data = H1_cost_dat, method = "REML")
summary(H1_cost_model)
coef_test(H1_cost_model, vcov = "CR2")

# saving statistics
H1_moderator_table <- rbind(H1_moderator_table,
                            slice_tail( # cost ratio has intercept
                              left_join(
                                coef_test(H1_cost_model, vcov = "CR2"),
                                conf_int(H1_cost_model, vcov = "CR2")
                              )
                            ))

### Study_Setting ----
H1_lab_dat <- H1_data %>% 
  filter(!is.na(Study_Setting))
V_H1_lab <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H1_lab_dat)
H1_lab_model <- rma.mv(yi, V_H1_lab,
                       mods = ~ Study_Setting - 1,
                       random = ~ 1 | No. / Study / effectID, data = H1_lab_dat, method = "REML")
summary(H1_lab_model)
Wald_test(H1_lab_model, constraints = constrain_equal(2:3), vcov = "CR2") # we didn't test field, cause it has only 2 effects, 1 paper
coef_test(H1_lab_model, vcov = "CR2")
conf_int(H1_lab_model, vcov = "CR2")

# saving forest plot to a file
H1_lab_col <- case_match(
  H1_lab_dat$Study_Setting,
  "Online" ~ "#004488",
  "Lab"    ~ "#BB5566",
  "Field"  ~ "#DDAA33"
)
jpeg("fig/H1_lab_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H1_lab_model, slab = paste(H1_lab_dat$Intext_Citation, "S", H1_lab_dat$Study, sep = " "),
       colout = H1_lab_col)
dev.off()

# saving statistics
H1_moderator_table <- rbind(H1_moderator_table,
                            left_join(
                              coef_test(H1_lab_model, vcov = "CR2"),
                              conf_int(H1_lab_model, vcov = "CR2")
                            ))

### Group_Type ----
H1_group_dat <- H1_data %>% 
  filter(!is.na(Group_Type))
V_H1_group <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H1_group_dat)
H1_group_model <- rma.mv(yi, V_H1_group,
                         mods = ~ Group_Type - 1,
                         random = ~ 1 | No. / Study / effectID, data = H1_group_dat, method = "REML")
summary(H1_group_model)
Wald_test(H1_group_model, constraints = constrain_equal(1:2), vcov = "CR2")
coef_test(H1_group_model, vcov = "CR2")
conf_int(H1_group_model, vcov = "CR2")

# saving statistics
H1_moderator_table <- rbind(H1_moderator_table,
                            left_join(
                              coef_test(H1_group_model, vcov = "CR2"),
                              conf_int(H1_group_model, vcov = "CR2")
                            ))

### Game_Type ----
H1_game_dat <- H1_data %>% 
  filter(!is.na(Game_Type))
V_H1_game <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H1_game_dat)
H1_game_model <- rma.mv(yi, V_H1_game,
                        mods = ~ Game_Type - 1,
                        random = ~ 1 | No. / Study / effectID, data = H1_game_dat, method = "REML")
summary(H1_game_model)
Wald_test(H1_game_model, constraints = constrain_equal(1:4), vcov = "CR2")
coef_test(H1_game_model, vcov = "CR2")
conf_int(H1_game_model, vcov = "CR2")

# saving forest plot to a file
H1_gam_col <- case_match(
  H1_game_dat$Game_Type,
  "Dictator Game" ~ "#004488",
  "Common Resource Dilemma"    ~ "#BB5566",
  "Prisoner's Dilemma"  ~ "#DDAA33"
)
jpeg("fig/H1_game_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H1_game_model, slab = paste(H1_game_dat$Intext_Citation, "S", H1_game_dat$Study, sep = " "),
       colout = H1_gam_col)
dev.off()

# saving statistics
H1_moderator_table <- rbind(H1_moderator_table,
                            left_join(
                              coef_test(H1_game_model, vcov = "CR2"),
                              conf_int(H1_game_model, vcov = "CR2")
                            ))

### Sex Ratio ----
H1_sex_dat <- H1_data %>% 
  mutate(male.ratio = Male / (Male + Female))
V_H1_sex <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H1_sex_dat)
H1_sex_model <- rma.mv(yi, V_H1_sex,
                       mods = ~ male.ratio,
                       random = ~ 1 | No. / Study / effectID, data = H1_sex_dat, method = "REML")
summary(H1_sex_model)
coef_test(H1_sex_model, vcov = "CR2")
conf_int(H1_sex_model, vcov = "CR2")

# saving statistics
H1_moderator_table <- rbind(H1_moderator_table,
                            slice_tail( # sex ratio has intercept
                              left_join(
                                coef_test(H1_sex_model, vcov = "CR2"),
                                conf_int(H1_sex_model, vcov = "CR2")
                              )
                            ))

### Age ----
H1_age_dat <- H1_data %>%
  mutate(Age_Mean = as.numeric(Age_Mean)) %>% 
  filter(!is.na(Age_Mean))
V_H1_age <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H1_age_dat)
H1_age_model <- rma.mv(yi, V_H1_age,
                       mods = ~ Age_Mean,
                       random = ~ 1 | No. / Study / effectID, data = H1_age_dat, method = "REML")
summary(H1_age_model)
coef_test(H1_age_model, vcov = "CR2")
conf_int(H1_age_model, vcov = "CR2")

# saving statistics
H1_moderator_table <- rbind(H1_moderator_table,
                            slice_tail( # age has intercept
                              left_join(
                                coef_test(H1_age_model, vcov = "CR2"),
                                conf_int(H1_age_model, vcov = "CR2")
                              )
                            ))

### Relative endowment ----
H1_endow_dat <- H1_data %>%
  mutate(Relative_Endowment = as.numeric(Relative_Endowment)) %>% 
  filter(!is.na(Relative_Endowment))
V_H1_endow <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H1_endow_dat)
H1_endow_model <- rma.mv(yi, V_H1_endow,
                         mods = ~ Relative_Endowment,
                         random = ~ 1 | No. / Study / effectID, data = H1_endow_dat, method = "REML")
summary(H1_endow_model)
coef_test(H1_endow_model, vcov = "CR2")
conf_int(H1_endow_model, vcov = "CR2")

# saving statistics
H1_moderator_table <- rbind(H1_moderator_table,
                            slice_tail( # Relative endowment has intercept
                              left_join(
                                coef_test(H1_endow_model, vcov = "CR2"),
                                conf_int(H1_endow_model, vcov = "CR2")
                              )
                            ))

### Real Interaction ----
H1_interact_dat <- H1_data %>% 
  filter(!is.na(Real_Interaction))
V_H1_interact <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H1_interact_dat)
H1_interact_model <- rma.mv(yi, V_H1_interact,
                            mods = ~ Real_Interaction - 1,
                            random = ~ 1 | No. / Study / effectID, data = H1_interact_dat, method = "REML")
summary(H1_interact_model)
Wald_test(H1_interact_model, constraints = constrain_equal(1:2), vcov = "CR2")
coef_test(H1_interact_model, vcov = "CR2")
conf_int(H1_interact_model, vcov = "CR2")

# saving statistics
H1_moderator_table <- rbind(H1_moderator_table,
                            left_join(
                              coef_test(H1_interact_model, vcov = "CR2"),
                              conf_int(H1_interact_model, vcov = "CR2")
                            ))

### Student Sample ----
H1_stu_dat <- H1_data %>% 
  filter(!is.na(Student_Sample))
V_H1_stu <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H1_stu_dat)
H1_stu_model <- rma.mv(yi, V_H1_stu,
                       mods = ~ Student_Sample - 1,
                       random = ~ 1 | No. / Study / effectID, data = H1_stu_dat, method = "REML")
summary(H1_stu_model)
Wald_test(H1_stu_model, constraints = constrain_equal(1:2), vcov = "CR2")
coef_test(H1_stu_model, vcov = "CR2")
conf_int(H1_stu_model, vcov = "CR2")

# saving statistics
H1_moderator_table <- rbind(H1_moderator_table,
                            left_join(
                              coef_test(H1_stu_model, vcov = "CR2"),
                              conf_int(H1_stu_model, vcov = "CR2")
                            ))

### Write Table ----
# Save table in folder tab
if (!dir.exists("tab")) {
  dir.create("tab")
}
H1_moderator_table <- H1_moderator_table %>% 
  select(-df) %>% 
  rename(
    Moderator = "Coef",
    g = "beta",
    t = "tstat",
    df = "df_Satt",
    p = "p_Satt",
    CI_lower = "CI_L",
    CI_upper = "CI_U"
  ) 
H1_moderator_table %>% 
  nice_table() %>% 
  flextable::save_as_docx(., path = "tab/H1_moderator.docx")

## The Effect of Partner Player Membership on Altruistic Punishment ------------
H2_moderator_table <- tibble()
### The culture ----
H2_country_dat <- H2_data %>% 
  filter(!is.na(Country))
V_H2_country <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H2_country_dat)
H2_country_model_omnibus <- rma(yi, vi,
                                mods = ~ Country,
                                random = ~ 1 | effectID, data = H2_country_dat, method = "REML", test = "knha")
summary(H2_country_model_omnibus)

H2_country_model <- rma(yi, vi,
                        mods = ~ Country - 1,
                        random = ~ 1 | effectID, data = H2_country_dat, method = "REML", test = "knha")
summary(H2_country_model)

# saving forest plot to a file
H2_country_col <- case_match(
  H2_country_dat$Country,
  "China"  ~ "#CC6677",
  "Turkey" ~ "#DDCC77",
  "Argentina" ~ "#44AA99",
  "Switzerland" ~ "#88CCEE",
  "Papua New Guinea" ~ "#332288",
  "United States"     ~ "#999933"
)

jpeg("fig/H2_country_Forest.jpeg", height = 12, width = 10, units = "in", res = 300)
forest(H2_country_model, slab = paste(H2_country_dat$Intext_Citation, "S", H2_country_dat$Study, sep = " "),
       colout = H2_country_col, ilab = Country, ilab.lab = "Country", ilab.xpos = -3.5, xlim = c(-8, 3), alim = c(-2, 1))
dev.off()

# saving statistics
H2_moderator_table <- as.data.frame(coef(summary(H2_country_model)))

### The cost ----
H2_cost_dat <- H2_data %>%
  mutate(Cost_Punishment_Ratio = as.numeric(Cost_Punishment_Ratio)) %>%
  filter(!is.na(Cost_Punishment_Ratio))
V_H2_cost <- vcalc(H2_cost_dat$vi, cluster=H2_cost_dat$Sample_ID, rho = 0.5)
H2_cost_model <- rma(yi, vi,
                     mods = ~ Cost_Punishment_Ratio,
                     random = ~ 1 | effectID, data = H2_cost_dat, method = "REML", test = "knha")
summary(H2_cost_model)

H2_moderator_table <- rbind(H2_moderator_table, 
                            slice_tail( # sex ratio has intercept
                              as.data.frame(coef(summary(H2_cost_model))))
)

### Study_Setting ----
H2_lab_dat <- H2_data %>% 
  filter(!is.na(Study_Setting))
V_H2_lab <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H2_lab_dat)
H2_lab_model <- rma(yi, vi,
                    mods = ~ Study_Setting - 1,
                    random = ~ 1 | effectID, data = H2_lab_dat, method = "REML", test = "knha")
summary(H2_lab_model)
anova(H2_lab_model,
      X = rbind(c(0, 1, -1)))

# saving forest plot to a file
H2_lab_col <- case_match(
  H2_lab_dat$Study_Setting,
  "Online" ~ "#004488",
  "Lab"    ~ "#BB5566",
  "Field"  ~ "#DDAA33"
)
jpeg("fig/H2_lab_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H2_lab_model, slab = paste(H2_lab_dat$Intext_Citation, "S", H2_lab_dat$Study, sep = " "),
       colout = H2_lab_col)
dev.off()
# saving statistics

H2_moderator_table <- rbind(H2_moderator_table,
                            as.data.frame(coef(summary(H2_lab_model))))

### Group_Type ----
H2_group_dat <- H2_data %>% 
  filter(!is.na(Group_Type))
V_H2_group <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H2_group_dat)
H2_group_model <- rma(yi, vi,
                      mods = ~ Group_Type - 1,
                      random = ~ 1 | effectID, data = H2_group_dat, method = "REML", test = "knha")
summary(H2_group_model)
anova(H2_group_model,
      X = rbind(c(1, -1)))

# saving statistics
H2_moderator_table <- rbind(H2_moderator_table,
                            as.data.frame(coef(summary(H2_group_model))))

### Game_Type ----
H2_game_dat <- H2_data %>% 
  filter(!is.na(Game_Type))
V_H2_game <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H2_game_dat)
H2_game_model <- rma(yi, vi,
                     mods = ~ Game_Type - 1,
                     random = ~ 1 | effectID, data = H2_game_dat, method = "REML", test = "knha")
summary(H2_game_model)
anova(H2_game_model,
      X = rbind(c(1, -1)))

# saving statistics
H2_moderator_table <- rbind(H2_moderator_table,
                            as.data.frame(coef(summary(H2_game_model))))

### Sex Ratio ----
H2_sex_dat <- H2_data %>% 
  mutate(male.ratio = Male / (Male + Female))
V_H2_sex <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H2_sex_dat)
H2_sex_model <- rma(yi, vi,
                    mods = ~ male.ratio,
                    random = ~ 1 | effectID, data = H2_sex_dat, method = "REML", test = "knha")
summary(H2_sex_model)

# saving statistics
H2_moderator_table <- rbind(H2_moderator_table,
                            slice_tail( # sex ratio has intercept
                              as.data.frame(coef(summary(H2_sex_model))))
)

### Age ----
H2_age_dat <- H2_data %>% 
  filter(!is.na(Age_Mean)) %>%
  mutate(Age_Mean = as.numeric(Age_Mean))
V_H2_age <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H2_age_dat)
H2_age_model <- rma(yi, vi,
                    mods = ~ Age_Mean,
                    random = ~ 1 | effectID, data = H2_age_dat, method = "REML", test = "knha")
summary(H2_age_model)

# saving statistics
H2_moderator_table <- rbind(H2_moderator_table,
                            slice_tail( # age has intercept
                              as.data.frame(coef(summary(H2_age_model))))
)

### relative endowment ----
H2_endow_dat <- H2_data %>% 
  filter(!is.na(Relative_Endowment)) %>%
  mutate(Relative_Endowment = as.numeric(Relative_Endowment))
V_H2_endow <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H2_endow_dat)
H2_endow_model <- rma(yi, vi,
                      mods = ~ Relative_Endowment,
                      random = ~ 1 | effectID, data = H2_endow_dat, method = "REML", test = "knha")
summary(H2_endow_model)

# saving statistics
H2_moderator_table <- rbind(H2_moderator_table,
                            slice_tail( # Relative endowment has intercept
                              as.data.frame(coef(summary(H2_endow_model))))
)

### Real Interaction ----
H2_interact_dat <- H2_data %>% 
  filter(!is.na(Real_Interaction))
V_H2_interact <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H2_interact_dat)
H2_interact_model <- rma(yi, vi,
                         mods = ~ Real_Interaction - 1,
                         random = ~ 1 | effectID, data = H2_interact_dat, method = "REML", test = "knha")
summary(H2_interact_model)
anova(H2_interact_model,
      X = rbind(c(1, -1)))

# saving statistics
H2_moderator_table <- rbind(H2_moderator_table,
                            as.data.frame(coef(summary(H2_interact_model))))

### Student Sample ----
H2_stu_dat <- H2_data %>% 
  filter(!is.na(Student_Sample))
V_H2_stu <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H2_stu_dat)
H2_stu_model <- rma(yi, vi,
                    mods = ~ Student_Sample - 1,
                    random = ~ 1 | effectID, data = H2_stu_dat, method = "REML", test = "knha")
summary(H2_stu_model) # NOTICE：this moderator is completely correlated with Real Interaction
anova(H2_stu_model,
      X = rbind(c(1, -1)))

# saving statistics
H2_moderator_table <- rbind(H2_moderator_table,
                            as.data.frame(coef(summary(H2_stu_model))))

### Write Table ----
H2_moderator_table <- H2_moderator_table %>% 
  rownames_to_column("Moderator") %>% 
  rename(
    g = "estimate",
    t = "tval",
    p = "pval",
    CI_lower = "ci.lb",
    CI_upper = "ci.ub"
  ) 
H2_moderator_table %>% 
  nice_table() %>% 
  flextable::save_as_docx(., path = "tab/H2_moderator.docx")


## The Effect of Focal Player on Antisocial Punishment ----
H5_moderator_table <- tibble()
### The culture ----
H5_country_dat <- H5_data %>% 
  filter(!is.na(Country)) 
V_H5_country <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H5_country_dat)
H5_country_model <- rma.mv(yi, V_H5_country, 
                           mods = ~ Country - 1,
                           random = ~ 1 | No. / effectID, data = H5_country_dat, method = "REML")
summary(H5_country_model)
Wald_test(H5_country_model, constraints = constrain_equal(1:6), vcov = "CR2", test = "EDT") # HTZ returns NA
coef_test(H5_country_model, vcov = "CR2")
conf_int(H5_country_model, vcov = "CR2")

# saving forest plot to a file
H5_country_col <- case_match(
  H5_country_dat$Country,
  "China"  ~ "#CC6677",
  "Turkey" ~ "#DDCC77",
  "United States"     ~ "#999933",
  "United Kingdom"  ~ "#117733",
  "Argentina" ~ "#44AA99",
  "Switzerland" ~ "#88CCEE",
  "Papua New Guinea" ~ "#332288",
  "Italy" ~ "#882255"
)
jpeg("fig/H5_country_Forest.jpeg", height = 12, width = 10, units = "in", res = 300)
forest(H5_country_model, slab = paste(H5_country_dat$Intext_Citation, "S", H5_country_dat$Study, sep = " "),
       colout = H5_country_col, ilab = Country, ilab.lab = "Country", ilab.xpos = -3.8, xlim = c(-9, 4), alim = c(-2.5, 2))
dev.off()

# saving statistics
H5_moderator_table <- left_join(
  coef_test(H5_country_model, vcov = "CR2"),
  conf_int(H5_country_model, vcov = "CR2")
)

### The cost ----
H5_cost_dat <- H5_data %>% 
  mutate(Cost_Punishment_Ratio = as.numeric(Cost_Punishment_Ratio)) %>% 
  filter(!is.na(Cost_Punishment_Ratio))
V_H5_cost <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H5_cost_dat)
H5_cost_model <- rma.mv(yi, V_H5_cost,
                        mods = ~ Cost_Punishment_Ratio,
                        random = ~ 1 | No. / effectID, data = H5_cost_dat, method = "REML")
summary(H5_cost_model)
coef_test(H5_cost_model, vcov = "CR2")
conf_int(H5_cost_model, vcov = "CR2")

# saving statistics
H5_moderator_table <- rbind(H5_moderator_table,
                            slice_tail( # cost ratio has intercept
                              left_join(
                                coef_test(H5_cost_model, vcov = "CR2"),
                                conf_int(H5_cost_model, vcov = "CR2")
                              )
                            ))

### Study_Setting ----
H5_lab_dat <- H5_data %>% 
  filter(!is.na(Study_Setting))
V_H5_lab <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H5_lab_dat)
H5_lab_model <- rma.mv(yi, V_H5_lab,
                       mods = ~ Study_Setting - 1,
                       random = ~ 1 | No. / effectID, data = H5_lab_dat, method = "REML")
summary(H5_lab_model)
Wald_test(H5_lab_model, constraints = constrain_equal(2:3), vcov = "CR2") # we didn't test field, cause it has only 2 effects, 1 paper
coef_test(H5_lab_model, vcov = "CR2")
conf_int(H5_lab_model, vcov = "CR2")

# saving statistics
H5_moderator_table <- rbind(H5_moderator_table,
                            left_join(
                              coef_test(H5_lab_model, vcov = "CR2"),
                              conf_int(H5_lab_model, vcov = "CR2")
                            ))

### Group_Type ----
H5_group_dat <- H5_data %>% 
  filter(!is.na(Group_Type))
V_H5_group <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H5_group_dat)
H5_group_model <- rma.mv(yi, V_H5_group,
                         mods = ~ Group_Type - 1,
                         random = ~ 1 | No. / effectID, data = H5_group_dat, method = "REML")
summary(H5_group_model)
Wald_test(H5_group_model, constraints = constrain_equal(1:2), vcov = "CR2")
coef_test(H5_group_model, vcov = "CR2")
conf_int(H5_group_model, vcov = "CR2")

# saving forest plot to a file
H5_group_col <- ifelse(H5_group_dat$Group_Type == "Natural", "#004488", "#BB5566")
jpeg("fig/H5_group_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H5_group_model, slab = paste(H5_group_dat$Intext_Citation, "S", H5_group_dat$Study, sep = " "),
       colout = H5_group_col)
dev.off()

# saving statistics
H5_moderator_table <- rbind(H5_moderator_table,
                            left_join(
                              coef_test(H5_group_model, vcov = "CR2"),
                              conf_int(H5_group_model, vcov = "CR2")
                            ))

### Game_Type ----
H5_game_dat <- H5_data %>% 
  filter(!is.na(Game_Type))
V_H5_game <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H5_game_dat)
H5_game_model <- rma.mv(yi, V_H5_game,
                        mods = ~ Game_Type - 1,
                        random = ~ 1 | No. / effectID, data = H5_game_dat, method = "REML")
summary(H5_game_model)
Wald_test(H5_game_model, constraints = constrain_equal(1:3), vcov = "CR2")

# saving statistics
H5_moderator_table <- rbind(H5_moderator_table,
                            left_join(
                              coef_test(H5_game_model, vcov = "CR2"),
                              conf_int(H5_game_model, vcov = "CR2")
                            ))

### Sex Ratio ----
H5_sex_dat <- H5_data %>% 
  mutate(male.ratio = Male / (Male + Female))
V_H5_sex <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H5_sex_dat)
H5_sex_model <- rma.mv(yi, V_H5_sex,
                       mods = ~ male.ratio,
                       random = ~ 1 | No. / effectID, data = H5_sex_dat, method = "REML")
summary(H5_sex_model)
coef_test(H5_sex_model, vcov = "CR2")

# saving statistics
H5_moderator_table <- rbind(H5_moderator_table,
                            slice_tail( # sex ratio has intercept
                              left_join(
                                coef_test(H5_sex_model, vcov = "CR2"),
                                conf_int(H5_sex_model, vcov = "CR2")
                              )
                            ))

### Age ----
H5_age_dat <- H5_data %>%
  mutate(Age_Mean = as.numeric(Age_Mean)) %>% 
  filter(!is.na(Age_Mean))
V_H5_age <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H5_age_dat)
H5_age_model <- rma.mv(yi, V_H5_age,
                       mods = ~ Age_Mean,
                       random = ~ 1 | No. / effectID, data = H5_age_dat, method = "REML")
summary(H5_age_model)
coef_test(H5_age_model, vcov = "CR2")
conf_int(H5_age_model, vcov = "CR2")

# saving statistics
H5_moderator_table <- rbind(H5_moderator_table,
                            slice_tail( # age has intercept
                              left_join(
                                coef_test(H5_age_model, vcov = "CR2"),
                                conf_int(H5_age_model, vcov = "CR2")
                              )
                            ))

### relative endowment ----
H5_endow_dat <- H5_data %>%
  mutate(Relative_Endowment = as.numeric(Relative_Endowment)) %>% 
  filter(!is.na(Relative_Endowment))
V_H5_endow <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H5_endow_dat)
H5_endow_model <- rma.mv(yi, V_H5_endow,
                         mods = ~ Relative_Endowment,
                         random = ~ 1 | No. / effectID, data = H5_endow_dat, method = "REML")
summary(H5_endow_model)
coef_test(H5_endow_model, vcov = "CR2")

# saving statistics
H5_moderator_table <- rbind(H5_moderator_table,
                            slice_tail( # Relative endowment has intercept
                              left_join(
                                coef_test(H5_endow_model, vcov = "CR2"),
                                conf_int(H5_endow_model, vcov = "CR2")
                              )
                            ))

### Real Interaction ----
H5_interact_dat <- H5_data %>% 
  filter(!is.na(Real_Interaction))
V_H5_interact <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H5_interact_dat)
H5_interact_model <- rma.mv(yi, V_H5_interact,
                            mods = ~ Real_Interaction - 1,
                            random = ~ 1 | No. / effectID, data = H5_interact_dat, method = "REML")
summary(H5_interact_model)
Wald_test(H5_interact_model, constraints = constrain_equal(1:2), vcov = "CR2")
coef_test(H5_interact_model, vcov = "CR2")

# saving statistics
H5_moderator_table <- rbind(H5_moderator_table,
                            left_join(
                              coef_test(H5_interact_model, vcov = "CR2"),
                              conf_int(H5_interact_model, vcov = "CR2")
                            ))

### Student Sample ----
H5_stu_dat <- H5_data %>% 
  filter(!is.na(Student_Sample))
V_H5_stu <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H5_stu_dat)
H5_stu_model <- rma.mv(yi, V_H5_stu,
                       mods = ~ Student_Sample - 1,
                       random = ~ 1 | No. / effectID, data = H5_stu_dat, method = "REML")
summary(H5_stu_model)
Wald_test(H5_stu_model, constraints = constrain_equal(1:2), vcov = "CR2")
coef_test(H5_stu_model, vcov = "CR2")
conf_int(H5_stu_model, vcov = "CR2")

# saving statistics
H5_moderator_table <- rbind(H5_moderator_table,
                            left_join(
                              coef_test(H5_stu_model, vcov = "CR2"),
                              conf_int(H5_stu_model, vcov = "CR2")
                            ))

### Write Table ----
# Save table in folder tab
H5_moderator_table <- H5_moderator_table %>% 
  select(-df) %>% 
  rename(
    Moderator = "Coef",
    g = "beta",
    t = "tstat",
    df = "df_Satt",
    p = "p_Satt",
    CI_lower = "CI_L",
    CI_upper = "CI_U"
  ) 
H5_moderator_table %>% 
  nice_table() %>% 
  flextable::save_as_docx(., path = "tab/H5_moderator.docx")



## Effect of Partner Player Membership on Antisocial Punishment ----
H6_moderator_table <- tibble()
### The culture ----
H6_country_dat <- H6_data %>% 
  filter(!is.na(Country)) 
V_H6_country <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H6_country_dat)
H6_country_model_omnibus <- rma.mv(yi, vi,
                                   mods = ~ Country - 1,
                                   random = ~ 1 | effectID, data = H6_country_dat, method = "REML", test = "knha")
summary(H6_country_model_omnibus)

H6_country_model <- rma.mv(yi, vi,
                           mods = ~ Country - 1,
                           random = ~ 1 | effectID, data = H6_country_dat, method = "REML", test = "knha")
summary(H6_country_model)

# saving statistics
H6_moderator_table <- as.data.frame(coef(summary(H6_country_model)))

### The cost ----
H6_cost_dat <- H6_data %>% 
  mutate(Cost_Punishment_Ratio = as.numeric(Cost_Punishment_Ratio)) %>% 
  filter(!is.na(Cost_Punishment_Ratio))
V_H6_cost <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H6_cost_dat)
H6_cost_model <- rma(yi, vi,
                     mods = ~ Cost_Punishment_Ratio,
                     random = ~ 1 | effectID, data = H6_cost_dat, method = "REML", test = "knha")
summary(H6_cost_model)

# saving statistics
H6_moderator_table <- rbind(H6_moderator_table,
                            slice_tail( # cost ratio has intercept
                              as.data.frame(coef(summary(H6_cost_model)))
                            ))

### Study_Setting ----
H6_lab_dat <- H6_data %>% 
  filter(!is.na(Study_Setting))
V_H6_lab <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H6_lab_dat)
H6_lab_model <- rma(yi, vi,
                    mods = ~ Study_Setting - 1,
                    random = ~ 1 | effectID, data = H6_lab_dat, method = "REML", test = "knha")
summary(H6_lab_model)
anova(H6_lab_model,
      X = rbind(c(0, 1, -1)))

# saving statistics
H6_moderator_table <- rbind(H6_moderator_table,
                            as.data.frame(coef(summary(H6_lab_model))))

### Group_Type ----
# Not appliable, all data is Natural group
# H6_group_dat <- H6_data %>%
#   filter(!is.na(Group_Type))
# V_H6_group <- vcalc(H6_group_dat$vi, cluster=H6_group_dat$Sample_ID, rho = 0.5)
# H6_group_model <- rma.mv(yi, V_H6_group,
#                          mods = ~ Group_Type,
#                          random = ~ 1 | effectID, data = H6_group_dat, method = "REML")
# summary(H6_group_model)

### Game_Type ----
H6_game_dat <- H6_data %>% 
  filter(!is.na(Game_Type))
V_H6_game <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H6_game_dat)
H6_game_model <- rma(yi, vi, 
                     mods = ~ Game_Type - 1,
                     random = ~ 1 | effectID, data = H6_game_dat, method = "REML", test = "knha")
summary(H6_game_model)
anova(H6_game_model,
      X = rbind(c(1, -1)))

# saving statistics
H6_moderator_table <- rbind(H6_moderator_table,
                            as.data.frame(coef(summary(H6_game_model))))

### Sex Ratio ----
H6_sex_dat <- H6_data %>% 
  mutate(male.ratio = Male / (Male + Female))
V_H6_sex <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H6_sex_dat)
H6_sex_model <- rma(yi, vi,
                    mods = ~ male.ratio,
                    random = ~ 1 | effectID, data = H6_sex_dat, method = "REML", test = "knha")
summary(H6_sex_model)

# saving statistics
H6_moderator_table <- rbind(H6_moderator_table,
                            slice_tail( # sex ratio has intercept
                              as.data.frame(coef(summary(H6_sex_model))))
)

### Age ----
H6_age_dat <- H6_data %>% 
  filter(!is.na(Age_Mean)) %>%
  mutate(Age_Mean = as.numeric(Age_Mean))
V_H6_age <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H6_age_dat)
H6_age_model <- rma(yi, vi,
                    mods = ~ Age_Mean,
                    random = ~ 1 | effectID, data = H6_age_dat, method = "REML", test = "knha")
summary(H6_age_model)

# saving statistics
H6_moderator_table <- rbind(H6_moderator_table,
                            slice_tail( # age has intercept
                              as.data.frame(coef(summary(H6_age_model))))
)

### relative endowment ----
H6_endow_dat <- H6_data %>% 
  filter(!is.na(Relative_Endowment)) %>%
  mutate(Relative_Endowment = as.numeric(Relative_Endowment))
V_H6_endow <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H6_endow_dat)
H6_endow_model <- rma(yi, vi,
                      mods = ~ Relative_Endowment,
                      random = ~ 1 | effectID, data = H6_endow_dat, method = "REML", test = "knha")
summary(H6_endow_model)

# saving statistics
H6_moderator_table <- rbind(H6_moderator_table,
                            slice_tail( # Relative endowment has intercept
                              as.data.frame(coef(summary(H6_endow_model))))
)

### Real Interaction ----
H6_interact_dat <- H6_data %>% 
  filter(!is.na(Real_Interaction))
V_H6_interact <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H6_interact_dat)
H6_interact_model <- rma(yi, vi,
                         mods = ~ Real_Interaction - 1,
                         random = ~ 1 | effectID, data = H6_interact_dat, method = "REML", test = "knha")
summary(H6_interact_model)
anova(H6_interact_model,
      X = rbind(c(1, -1)))

# saving statistics
H6_moderator_table <- rbind(H6_moderator_table,
                            as.data.frame(coef(summary(H6_interact_model))))

### Student Sample ----
H6_stu_dat <- H6_data %>% 
  filter(!is.na(Student_Sample))
V_H6_stu <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H6_stu_dat)
H6_stu_model <- rma(yi, vi,
                    mods = ~ Student_Sample - 1,
                    random = ~ 1 | effectID, data = H6_stu_dat, method = "REML", test = "knha")
summary(H6_stu_model)
anova(H6_stu_model,
      X = rbind(c(1, -1)))

# saving statistics
H6_moderator_table <- rbind(H6_moderator_table,
                            as.data.frame(coef(summary(H6_stu_model))))

### Write Table ----
H6_moderator_table <- H6_moderator_table %>% 
  rownames_to_column("Moderator") %>% 
  rename(
    g = "estimate",
    t = "tval",
    p = "pval",
    CI_lower = "ci.lb",
    CI_upper = "ci.ub"
  ) 
H6_moderator_table %>% 
  # mutate(across(where(is.numeric), ~ round(., 4))) %>% 
  nice_table() %>% 
  flextable::save_as_docx(., path = "tab/H6_moderator.docx")



# Sensitivity analysis for the value of correlation of effect sizes obtained on the same group of participants ########################
H1_sensitivity <- tibble()
H5_sensitivity <- tibble()
# H1_data <- H1_data %>% 
#   mutate(StudyID = paste0(No.,"_",Study))
for (r in 1:10) {
  R = r/10 # real correlation coefficient
  ## H1
  V_H1_s <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = R, data = H1_data)
  H1_model_s <- rma.mv(yi, V_H1_s, random = ~ 1 | No. / Study / effectID, data = H1_data, method = "REML",
                       control=list(rel.tol=1e-9)) # by default is 1e-10
  # summary(H1_model_s)
  # saving statistics
  H1_sensitivity <- rbind(H1_sensitivity,
                          left_join(
                            coef_test(H1_model_s, vcov = "CR2"),
                            conf_int(H1_model_s, vcov = "CR2")
                          )
  )
  ## H5
  V_H5_s <- vcalc(vi, cluster = Sample_ID, obs = effectID, rho = R, data = H5_data)
  H5_model_s <- rma.mv(yi, V_H5_s, random = ~ 1 | No. / effectID, data = H5_data, method = "REML",
                       control=list(rel.tol=1e-9))
  
  # saving statistics
  H5_sensitivity <- rbind(H5_sensitivity,
                          left_join(
                            coef_test(H5_model_s, vcov = "CR2"),
                            conf_int(H5_model_s, vcov = "CR2")
                          )
  )
}

## write table ----
H1_sensitivity[,1] <- seq(from = 0.1, to = 1, by = 0.1)
H5_sensitivity[,1] <- seq(from = 0.1, to = 1, by = 0.1)

H1_sensitivity_short <- H1_sensitivity %>% 
  select(Coef, beta, p_Satt) %>% 
  rename(
    r = "Coef",
    H1.g = "beta",
    H1.p = "p_Satt"
  )

H5_sensitivity_short <- H5_sensitivity %>% 
  select(beta, p_Satt) %>% 
  rename(
    H5.g = "beta",
    H5.p = "p_Satt"
  )


cbind(H1_sensitivity_short, H5_sensitivity_short) %>% 
  nice_table(., separate.header = TRUE) %>% 
  flextable::save_as_docx(., path = "tab/sensitivity.docx")

# Publication Bias ########################
# we use the precision-effect test (PET) and the precision-effect estimate with standard error (PEESE)
## H1 ----
PET_H1 <- rma.mv(
  yi, V_H1,
  mods = ~ sqrt(vi),
  random = ~ 1 | No. / Study / effectID,
  data = H1_data,
  method = "REML"
)
summary(PET_H1)
coef_test(PET_H1, vcov = "CR2") # intercept is not significant, use PET

PETPEESE_CR2 <- slice_head(conf_int(PET_H1, vcov = "CR2"))

## H2 ----
PET_H2 <- rma(yi, vi, mods = ~ sqrt(vi), data = H2_data, method = "REML", test = "knha")
summary(PET_H2) # intercept is significant, use PEESE
PEESE_H2 <- rma(yi, vi, mods = ~ vi, data = H2_data, method = "REML", test = "knha")
summary(PEESE_H2)

PETPEESE_one <- slice_head(
  as.data.frame(coef(summary(PEESE_H2)))
  )

## H5 ----
PET_H5 <- rma.mv(
  yi, V_H5, 
  mods = ~ sqrt(vi),
  random = ~ 1 | No. / effectID, 
  data = H5_data, method = "REML")
summary(PET_H5)
coef_test(PET_H5, vcov = "CR2") # intercept is significantly positive, use PEESE

PEESE_H5 <- rma.mv(
  yi, V_H5,
  mods = ~ vi,
  random = ~ 1 | No. / effectID,
  data = H5_data,
  method = "REML"
)
summary(PEESE_H5)
coef_test(PEESE_H5, vcov = "CR2")

PETPEESE_CR2 <- rbind(PETPEESE_CR2,
                      slice_head(conf_int(PEESE_H5, vcov = "CR2"))
)

## H6 ----
PET_H6 <- rma(yi, vi, mods = ~ sqrt(vi), data = H6_data, method = "REML", test = "knha")
summary(PET_H6) # intercept is not significant, use PET

## Write Table ----
PETPEESE_one <- rbind(PETPEESE_one,
                      slice_head(as.data.frame(coef(summary(PET_H6))))
)

PETPEESE_CR2 <- PETPEESE_CR2 %>% 
  select(beta,CI_L,CI_U) %>% 
  rename(
    PETPEESE.g = "beta",
    PETPEESE.CI_lower = "CI_L",
    PETPEESE.CI_upper = "CI_U"
  ) %>% 
  add_column(
    Model = c("H1","H5")
  )

PETPEESE_one <- PETPEESE_one %>% 
  select(estimate, ci.lb, ci.ub) %>% 
  rename(
    PETPEESE.g = "estimate",
    PETPEESE.CI_lower = "ci.lb",
    PETPEESE.CI_upper = "ci.ub"
  ) %>% 
  add_column(
    Model = c("H2","H6")
  )
PETPEESE_table <- rbind(PETPEESE_CR2,
                        PETPEESE_one) %>% 
  relocate(Model)

main_CR2 <- rbind(
  conf_int(H1_model, vcov = "CR2"),
  conf_int(H5_model, vcov = "CR2")
) %>% 
  select(beta,CI_L,CI_U) %>% 
  rename(
    Original.g = "beta",
    Original.CI_lower = "CI_L",
    Original.CI_upper = "CI_U"
  ) %>% 
  add_column(
    Model = c("H1","H5")
  ) 

main_one <- rbind(
  as.data.frame(coef(summary(H2_model))),
  as.data.frame(coef(summary(H6_model)))
) %>% 
  select(estimate, ci.lb, ci.ub) %>% 
  rename(
    Original.g = "estimate",
    Original.CI_lower = "ci.lb",
    Original.CI_upper = "ci.ub"
  ) %>% 
  add_column(
    Model = c("H2","H6")
  )

main_table <- rbind(
  main_CR2, main_one
) %>% relocate(Model)

PETPEESE_table <- left_join(PETPEESE_table,
                            main_table)
  
PETPEESE_table %>% 
  nice_table(separate.header = TRUE,
             note = c(
               "For the multilevel models (i.e., H1 and H5), cluster-robust (CR2) standard errors were applied. Models highlighted in the table showed a significant intercept in the PET specification; accordingly, the corresponding PEESE estimates are reported as the bias-adjusted effects."
             )) %>% flextable::bg(bg = "gray", i = c(2,3)) %>% 
  flextable::save_as_docx(., path = "tab/publication_bias.docx")

# H4 and inequality aversion ----
ineq_data_fair <- MA_dat %>% 
  filter(
    `Action_A_Binary` == 1,
    `Membership_A` %in% c("Ingroup", "Outgroup")
  )
plot(ineq_data_fair$Relative_Endowment, ineq_data_fair$Punishment_Mean)

H4_endow_dat <- H4_data_r %>%
  filter(!is.na(Relative_Endowment))
V_H4_endow <- vcalc(SV_log, cluster = Sample_ID, obs = effectID, rho = 0.5, data = H4_endow_dat)
H4_endow_model <- H4_endow_dat %>%
  rma.mv(yi = P_log, V = V_H4_endow,
         mods = ~ Relative_Endowment - 1,
         random = ~ 1 | No. / effectID,
         data = .,
         method = "REML",
  )
robust(H4_endow_model, cluster = No., clubSandwich=TRUE)
Wald_test(H4_endow_model, constraints = constrain_equal(1:3), vcov = "CR2")


save.image(file = 'TPP_meta.rds', compress = "xz")

