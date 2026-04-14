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

# reading the data
MA_dat <- read_csv("Coding Sample - Final Meta Data.csv", skip = 1)
MA_dat$`Fairness/Cooperativeness A Binary` <- as.numeric(MA_dat$`Fairness/Cooperativeness A Binary`)

# Unfair Condition Analysis ----
## H1 The Effects of Violator’s Identity on Punishment ----
# fair_condition <- c("Fair", "Cooperate", "Neutral", "Generous")
# unfair_condition <- c("Unfair", "Defect", "Selfish", "Cheat")
### Data Preparation ----
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

### Model Fitting ----
H1_model <- rma.mv(yi, V_H1, random = ~ 1 | No. / Study / effectID, data = H1_data, method = "REML")
robust(H1_model, cluster = No., clubSandwich=TRUE)

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
forest(H1_model, slab = paste(H1_data$Intext.Citation, "S", H1_data$Study, sep = " "))
# Save fig in folder fig
if (!dir.exists("fig")) {
  dir.create("fig")
}
# saving forest plot to a file
jpeg("fig/H1_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H1_model, slab = paste(H1_data$Intext.Citation, "S", H1_data$Study, sep = " "))
dev.off()

# funnel plot
metafor::funnel(H1_model)

# contour-enhanced funnel plot
jpeg("fig/H1_funnel.jpeg", height = 6, width = 6, units = "in", res = 300)
metafor::funnel(H1_model, refline=0, level=c(95, 99), 
                shade=c("white", "grey"))
dev.off()

## H2 (Deterrence perspective): ----
### Data Preparation ----
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

### Model Fitting ----
H2_model <- rma.mv(yi, V_H2, random = ~ 1 | No. / Study / effectID, data = H2_data, method = "REML")
robust(H2_model, cluster = No., clubSandwich=TRUE)

# checking likelihood profile plots
profile(H2_model)

# calculating the overall amount of heterogeneity
sum(H2_model$sigma2)

# obtaining I2 value
W <- diag(1/H2_data$vi)
X <- model.matrix(H2_model)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(H2_model$sigma2) / (sum(H2_model$sigma2) + (H2_model$k-H2_model$p)/sum(diag(P))) # overall I2

# forest plot
H2_f <- forest(H2_model, slab = paste(H2_data$Intext.Citation, "S", H2_data$Study, sep = " "))
# saving forest plot to a file
jpeg("fig/H2_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H2_model, slab = paste(H2_data$Intext.Citation, "S", H2_data$Study, sep = " "))
dev.off()

# funnel plot
metafor::funnel(H2_model)

# contour-enhanced funnel plot
jpeg("fig/H2_funnel.jpeg", height = 6, width = 6, units = "in", res = 300)
metafor::funnel(H2_model, refline=0, level=c(95, 99), 
                shade=c("white", "grey"))
dev.off()

## H3 Interaction Effect between Violator’s Identity and Victim Salience ----
### The moderator method ----
H3_model <- H1_data %>%
  rma.mv(yi, V_H1,
  mods = ~ Recipient.Player.B.Affliation - 1,
  random = ~ 1 | No. / Study / effectID,
  data = .,
  method = "REML",
)
robust(H3_model, cluster = No., clubSandwich=TRUE)

Wald_test(H3_model, constrain_equal(1:3), vcov = "CR2")

# asign color for ingroup, outgroup, and neutral
H3_col <- ifelse(H1_data$Recipient.Player.B.Affliation == "Ingroup", "#004488",
                 ifelse(H1_data$Recipient.Player.B.Affliation == "Outgroup", "#BB5566", "gray"))
jpeg("fig/H3_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H3_model, slab = paste(H1_data$Intext.Citation, "S", H1_data$Study, sep = " "),
       colout = H3_col)
dev.off()

#### Robustness check 1: Exclude studies with neutral recipient/player B affiliation ----
H3_data <- H1_data %>%
  filter(Recipient.Player.B.Affliation != "Neutral")

V_H3 <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H3_data)

H3_model2 <-H3_data  %>%
  rma.mv(yi, V_H3,
         mods = ~ Recipient.Player.B.Affliation - 1,
         random = ~ 1 | No. / Study / effectID,
         data = .,
         method = "REML",
)
robust(H3_model2, cluster = No., clubSandwich=TRUE)

Wald_test(H3_model2, constraints = constrain_equal(1:2), vcov = "CR2")

#### Robustness check 2: Using the H2_data ----
H3_model3 <- H2_data %>%
  rma.mv(yi, vi, # V cannot convergence
         mods = ~ Person.A.affiliation - 1,
         random = ~ 1 | No. / Study / effectID,
         data = .,
         method = "REML",
)
robust(H3_model3, cluster = No., clubSandwich=TRUE)

Wald_test(H3_model3, constrain_equal(1:3), vcov = "CR2")

#### Robustness check 3: Using H2_data &Exclude studies with neutral recipient/player B affiliation ----
H3_data2 <- H2_data %>%
  filter(Person.A.affiliation != "Neutral")

V_H3_2 <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H3_data2)

H3_model4 <-H3_data2  %>%
  rma.mv(yi, V_H3_2,
         mods = ~ Person.A.affiliation - 1,
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
    `Fairness/Cooperativeness A Binary` == 0, # unfair condition
    `Person A affiliation` %in% c("Ingroup", "Outgroup"),
    `Recipient/Player B Affliation` %in% c("Ingroup", "Outgroup"),
    !is.na(`P_Average Invested Amount_Standardized`),
    !is.na(`P_SD invested amount _Standardized`),
    !is.na(`Total Trials`)
  ) %>%
  mutate(
    `Total Trials` = as.numeric(`Total Trials`),
    `Person A` = factor(`Person A affiliation`, levels = c("Ingroup", "Outgroup")),
    `Person B` = factor(`Recipient/Player B Affliation`, levels = c("Ingroup", "Outgroup")),
    P_SV = as.numeric(`P_SD invested amount _Standardized`)^2 * `Total Trials` / ( `Total Trials` - 1), # Sampling variance of punishment amount
    effect_n = row_number(),
    effectID = paste(No., Study, effect_n, sep = "_"),
    P_log = car::logit(`P_Average Invested Amount_Standardized`),
    SV_log = `P_SD invested amount`^2 / (`P_Average Invested Amount`^2 * `Total Trials` * (1-`P_Average Invested Amount_Standardized`)^2)
  )

V_H3  <- vcalc(SV_log, cluster = `Sample ID`, obs = effectID, rho = 0.5, data = H3_data_r)

#### Model Fitting ----
H3_model_r <- H3_data_r %>%
  rma.mv(yi = P_log, V = V_H3, # V close to singular
         mods = ~ `Person A` * `Person B`,
         random = ~ 1 | No. / Study / effectID,  data = .,
         method = "REML",
)
robust(H3_model_r, cluster = No., clubSandwich=TRUE)

#### Estimated Average Effects ----
H3_data_d <- MA_dat %>%
  filter(
    `Fairness/Cooperativeness A Binary` == 0, # unfair condition
    !is.na(`P_Average Invested Amount_Standardized`),
    !is.na(`P_SD invested amount _Standardized`),
    !is.na(`Total Trials`)
  ) %>%
  mutate(
    `Total Trials` = as.numeric(`Total Trials`),
    `Person A` = factor(`Person A affiliation`, levels = c("Ingroup", "Neutral", "Outgroup")),
    `Person B` = factor(`Recipient/Player B Affliation`, levels = c("Ingroup", "Neutral", "Outgroup")),
    P_SV = as.numeric(`P_SD invested amount _Standardized`)^2 * `Total Trials` / ( `Total Trials` - 1), # Sampling variance of punishment amount
    effect_n = row_number(),
    effectID = paste(No., Study, effect_n, sep = "_"),
    P_log = car::logit(`P_Average Invested Amount_Standardized`),
    SV_log = `P_SD invested amount`^2 / (`P_Average Invested Amount`^2 * `Total Trials` * (1-`P_Average Invested Amount_Standardized`)^2)
  )

V_H3_d  <- vcalc(SV_log, cluster = `Sample ID`, obs = effectID, rho = 0.5, data = H3_data_d)
H3_model_d <- H3_data_d %>%
  rma.mv(yi = P_log, V = V_H3_d,
         mods = ~ `Person A` * `Person B`,
         random = ~ 1 | No. / Study / effectID,  data = .,
         method = "REML",
  )
robust(H3_model_d, cluster = No., clubSandwich=TRUE)

results_mat_name <- expand.grid(`Person A` = c("Ingroup", "Neutral", "Outgroup"), `Person B` = c("Ingroup", "Neutral", "Outgroup"))
results_mat <- model.matrix(~ `Person A` * `Person B`, results_mat_name)[,-1]
results_dat <- predict(H3_model_d, newmods = results_mat) %>% 
  cbind(results_mat_name) %>% 
  mutate(pred = boot::inv.logit(pred))

heat_all <- results_dat %>% 
  ggplot(aes(x = `Person A`, y = `Person B`, fill = pred)) +
  geom_tile() +
  geom_text(aes(label = paste0(round(pred * 100,1), "%"))) +
  scale_fill_gradient(low = "white", high = "#BB5566",
                      limits = c(0, 0.5)) +
  theme_cowplot() +
  labs(fill = "Punishment", x = "Focal Player Group Membership",
       y = "Partner Player Group Membership")
ggsave("fig/Summary_unfair_heatmap.jpeg", heat_all, height = 6, width = 7)

H4_data_d <- MA_dat %>% # not appropriate, but let's say it is H4 temporary
  filter(
    `Fairness/Cooperativeness A Binary` == 1, # fair condition
    !is.na(`P_Average Invested Amount_Standardized`),
    !is.na(`P_SD invested amount _Standardized`),
    !is.na(`Total Trials`)
  ) %>%
  mutate(
    `Total Trials` = as.numeric(`Total Trials`),
    `Person A` = factor(`Person A affiliation`, levels = c("Ingroup", "Neutral", "Outgroup")),
    `Person B` = factor(`Recipient/Player B Affliation`, levels = c("Ingroup", "Neutral", "Outgroup")),
    P_SV = as.numeric(`P_SD invested amount _Standardized`)^2 * `Total Trials` / ( `Total Trials` - 1), # Sampling variance of punishment amount
    effect_n = row_number(),
    effectID = paste(No., Study, effect_n, sep = "_"),
    P_log = car::logit(`P_Average Invested Amount_Standardized`),
    SV_log = `P_SD invested amount`^2 / (`P_Average Invested Amount`^2 * `Total Trials` * (1-`P_Average Invested Amount_Standardized`)^2),
    SV_log = ifelse(is.na(SV_log), 1/`Total Trials`, SV_log)
  )

V_H4_d  <- vcalc(SV_log, cluster = `Sample ID`, obs = effectID, rho = 0.5, data = H4_data_d)
H4_model_d <- H4_data_d %>%
  rma.mv(yi = P_log, V = V_H4_d, # V close to singular
         mods = ~ `Person A` * `Person B`,
         random = ~ 1 | No. / Study / effectID,  data = .,
         method = "REML",
  )
summary(H4_model_d)
results_mat2 <- results_mat[-c(5,8),-c(5,7)]
results_dat2 <- predict(H4_model_d, newmods = results_mat2) %>% 
  cbind(results_mat_name[-c(5,8),]) %>% 
  right_join(results_mat_name) %>% 
  mutate(pred = boot::inv.logit(pred),
         pred = (pred - 0.025) / (1 - 2 * 0.025) # Note that due to 0 in the data, car::logit() remapped data, so we need invert it as well
         )

heat_all_2 <- results_dat2 %>% 
  ggplot(aes(x = `Person A`, y = `Person B`, fill = pred)) +
  geom_tile() +
  geom_text(aes(label = ifelse(is.na(pred),"NA",
                               paste0(round(pred * 100,1), "%"))
                )) +
  scale_fill_gradient(low = "white", high = "#BB5566",
                      limits = c(0, 0.5)) +
  theme_cowplot() +
  labs(fill = "Punishment", x = "Focal Player Group Membership",
       y = "Partner Player Group Membership")
ggsave("fig/Summary_fair_heatmap.jpeg", heat_all_2, height = 6, width = 7)

heat_all_all <- plot_grid(heat_all, heat_all_2, labels = c("a", "b"))
ggsave("fig/Summary_all_heatmap.jpeg", heat_all_all, height = 6, width = 14)

# Fair Condition Analysis ################
## H5 The Effects of Violator’s Identity on Punishment ----
### Data Preparation ----
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

### Model Fitting ----
H5_model <- rma.mv(yi, V_H5, random = ~ 1 | No. / Study / effectID, data = H5_data, method = "REML")
robust(H5_model, cluster = No., clubSandwich=TRUE)

# calculating the overall amount of heterogeneity
sum(H5_model$sigma2)

# obtaining I2 value
W <- diag(1/H5_data$vi)
X <- model.matrix(H5_model)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(H5_model$sigma2) / (sum(H5_model$sigma2) + (H5_model$k-H5_model$p)/sum(diag(P))) # overall I2
100 * H5_model$sigma2 / (sum(H5_model$sigma2) + (H5_model$k-H5_model$p)/sum(diag(P)))  # I2 per levels

# forest plot
H5_f <- forest(H5_model, slab = paste(H5_data$Intext.Citation, "S", H5_data$Study, sep = " "))
# saving forest plot to a file
jpeg("fig/H5_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H5_model, slab = paste(H5_data$Intext.Citation, "S", H5_data$Study, sep = " "))
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

### Model Fitting ----
H6_model <- rma.mv(yi, V_H6,
                   random = ~ 1 | No. / Study / effectID, data = H6_data, method = "REML")
robust(H6_model, cluster = No., clubSandwich=TRUE)

# forest plot
H6_f <- forest(H6_model, slab = paste(H6_data$Intext.Citation, "S", H6_data$Study, sep = " "))
# saving forest plot to a file
jpeg("fig/H6_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H6_model, slab = paste(H6_data$Intext.Citation, "S", H6_data$Study, sep = " "))
dev.off()

# calculating the overall amount of heterogeneity
sum(H6_model$sigma2)

# obtaining I2 value
W <- diag(1/H6_data$vi)
X <- model.matrix(H6_model)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(H6_model$sigma2) / (sum(H6_model$sigma2) + (H6_model$k-H6_model$p)/sum(diag(P))) # overall I2

# funnel plot
metafor::funnel(H6_model)

# contour-enhanced funnel plot
jpeg("fig/H6_funnel.jpeg", height = 6, width = 6, units = "in", res = 300)
metafor::funnel(H6_model, refline=0, level=c(95, 99), 
                shade=c("white", "grey"))
dev.off()


# Single moderator analyses ################
## H1 The Effects of Violator’s Identity on Punishment moderator analyses ----
H1_moderator_table <- tibble()
### The culture ----
H1_country_dat <- H1_data %>% 
  filter(!is.na(Country))
V_H1_country <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H1_country_dat)
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
  "US"     ~ "#999933",
  "UK/US"  ~ "#117733",
  "Argentina" ~ "#44AA99",
  "Switzerland" ~ "#88CCEE",
  "Papua New Guinea" ~ "#332288"
)

jpeg("fig/H1_country_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H1_country_model, slab = paste(H1_country_dat$Intext.Citation, "S", H1_country_dat$Study, sep = " "),
       colout = H1_country_col)
dev.off()

# saving statistics
H1_moderator_table <- left_join(
  coef_test(H1_country_model, vcov = "CR2"),
  conf_int(H1_country_model, vcov = "CR2")
)

### The cost ----
H1_cost_dat <- H1_data %>% 
  mutate(Punishment.Cost.Ratio = as.numeric(Punishment.Cost.Ratio)) %>% 
  filter(!is.na(Punishment.Cost.Ratio))
V_H1_cost <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H1_cost_dat)
H1_cost_model <- rma.mv(yi, V_H1_cost,
                           mods = ~ Punishment.Cost.Ratio,
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

### Lab or online ----
H1_lab_dat <- H1_data %>% 
  filter(!is.na(Lab.or.Online))
V_H1_lab <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H1_lab_dat)
H1_lab_model <- rma.mv(yi, V_H1_lab,
                        mods = ~ Lab.or.Online - 1,
                        random = ~ 1 | No. / Study / effectID, data = H1_lab_dat, method = "REML")
summary(H1_lab_model)
Wald_test(H1_lab_model, constraints = constrain_equal(1:2), vcov = "CR2")
coef_test(H1_lab_model, vcov = "CR2")
conf_int(H1_lab_model, vcov = "CR2")

# saving forest plot to a file
H1_lab_col <- ifelse(H1_lab_dat$Lab.or.Online == "Online", "#004488", "#BB5566")
jpeg("fig/H1_lab_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H1_lab_model, slab = paste(H1_lab_dat$Intext.Citation, "S", H1_lab_dat$Study, sep = " "),
       colout = H1_lab_col)
dev.off()

# saving statistics
H1_moderator_table <- rbind(H1_moderator_table,
                            left_join(
                                coef_test(H1_lab_model, vcov = "CR2"),
                                conf_int(H1_lab_model, vcov = "CR2")
                                ))

### Group Manipulation ----
H1_group_dat <- H1_data %>% 
  filter(!is.na(Group.Manipulation))
V_H1_group <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H1_group_dat)
H1_group_model <- rma.mv(yi, V_H1_group,
                       mods = ~ Group.Manipulation - 1,
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

### Game Type ----
H1_game_dat <- H1_data %>% 
  filter(!is.na(Game.Type))
V_H1_game <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H1_game_dat)
H1_game_model <- rma.mv(yi, V_H1_game,
                         mods = ~ Game.Type - 1,
                         random = ~ 1 | No. / Study / effectID, data = H1_game_dat, method = "REML")
summary(H1_game_model)
Wald_test(H1_game_model, constraints = constrain_equal(1:4), vcov = "CR2")
coef_test(H1_game_model, vcov = "CR2")
conf_int(H1_game_model, vcov = "CR2")

# saving statistics
H1_moderator_table <- rbind(H1_moderator_table,
                            left_join(
                              coef_test(H1_game_model, vcov = "CR2"),
                              conf_int(H1_game_model, vcov = "CR2")
                            ))

### Sex Ratio ----
H1_sex_dat <- H1_data %>% 
  mutate(male.ratio = Male / (Male + Female))
V_H1_sex <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H1_sex_dat)
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
  mutate(Age..mean. = as.numeric(Age..mean.)) %>% 
  filter(!is.na(Age..mean.))
V_H1_age <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H1_age_dat)
H1_age_model <- rma.mv(yi, V_H1_age,
                         mods = ~ Age..mean.,
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
  mutate(Relative.Endowments = as.numeric(Relative.Endowments)) %>% 
  filter(!is.na(Relative.Endowments))
V_H1_endow <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H1_endow_dat)
H1_endow_model <- rma.mv(yi, V_H1_endow,
                       mods = ~ Relative.Endowments,
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
  filter(!is.na(Real.or.Deception))
V_H1_interact <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H1_interact_dat)
H1_interact_model <- rma.mv(yi, V_H1_interact,
                        mods = ~ Real.or.Deception - 1,
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
  filter(!is.na(Student.Sample))
V_H1_stu <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H1_stu_dat)
H1_stu_model <- rma.mv(yi, V_H1_stu,
                       mods = ~ Student.Sample - 1,
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

## H2 (Deterrence perspective): ----
H2_moderator_table <- tibble()
### The culture ----
H2_country_dat <- H2_data %>% 
  filter(!is.na(Country))
V_H2_country <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H2_country_dat)
H2_country_model <- rma.mv(yi, V_H2_country,
                           mods = ~ Country - 1,
                           random = ~ 1 | No. / Study / effectID, data = H2_country_dat, method = "REML")
summary(H2_country_model)
Wald_test(H2_country_model, constraints = constrain_equal(1:5), vcov = "CR2", test = "EDT") # HTZ returns NA
coef_test(H2_country_model, vcov = "CR2")
conf_int(H2_country_model, vcov = "CR2")

# saving forest plot to a file
H2_country_col <- case_match(
  H2_country_dat$Country,
  "China"  ~ "#CC6677",
  "Turkey" ~ "#DDCC77",
  "Argentina" ~ "#44AA99",
  "Switzerland" ~ "#88CCEE",
  "Papua New Guinea" ~ "#332288"
)

jpeg("fig/H2_country_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H2_country_model, slab = paste(H2_country_dat$Intext.Citation, "S", H2_country_dat$Study, sep = " "),
       colout = H2_country_col)
dev.off()

# saving statistics
H2_moderator_table <- left_join(
  coef_test(H2_country_model, vcov = "CR2"),
  conf_int(H2_country_model, vcov = "CR2")
)

### The cost ----
# Not appliable, all data is 3
# H2_cost_dat <- H2_data %>% 
#   mutate(Punishment.Cost.Ratio = as.numeric(Punishment.Cost.Ratio)) %>% 
#   filter(!is.na(Punishment.Cost.Ratio))
# V_H2_cost <- vcalc(H2_cost_dat$vi, cluster=H2_cost_dat$Sample.ID, rho = 0.5)
# H2_cost_model <- rma.mv(yi, V_H2_cost,
#                         mods = ~ Punishment.Cost.Ratio,
#                         random = ~ 1 | No. / Study / effectID, data = H2_cost_dat, method = "REML")
# summary(H2_cost_model)
# coef_test(H2_cost_model, vcov = "CR2")
# conf_int(H2_cost_model, vcov = "CR2")

### Lab or online ----
H2_lab_dat <- H2_data %>% 
  filter(!is.na(Lab.or.Online))
V_H2_lab <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H2_lab_dat)
H2_lab_model <- rma.mv(yi, V_H2_lab,
                       mods = ~ Lab.or.Online - 1,
                       random = ~ 1 | No. / Study / effectID, data = H2_lab_dat, method = "REML")
summary(H2_lab_model)
Wald_test(H2_lab_model, constraints = constrain_equal(1:2), vcov = "CR2")
coef_test(H2_lab_model, vcov = "CR2")
conf_int(H2_lab_model, vcov = "CR2")

# saving forest plot to a file
H2_lab_col <- ifelse(H2_lab_dat$Lab.or.Online == "Online", "#004488", "#BB5566")
jpeg("fig/H2_lab_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H2_lab_model, slab = paste(H2_lab_dat$Intext.Citation, "S", H2_lab_dat$Study, sep = " "),
       colout = H2_lab_col)
dev.off()
# saving statistics

H2_moderator_table <- rbind(H2_moderator_table,
                            left_join(
                              coef_test(H2_lab_model, vcov = "CR2"),
                              conf_int(H2_lab_model, vcov = "CR2")
                            ))

### Group Manipulation ----
H2_group_dat <- H2_data %>% 
  filter(!is.na(Group.Manipulation))
V_H2_group <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H2_group_dat)
H2_group_model <- rma.mv(yi, V_H2_group,
                         mods = ~ Group.Manipulation - 1,
                         random = ~ 1 | No. / Study / effectID, data = H2_group_dat, method = "REML")
summary(H2_group_model)
Wald_test(H2_group_model, constraints = constrain_equal(1:2), vcov = "CR2")
coef_test(H2_group_model, vcov = "CR2")
conf_int(H2_group_model, vcov = "CR2")

# saving statistics
H2_moderator_table <- rbind(H2_moderator_table,
                            left_join(
                              coef_test(H2_group_model, vcov = "CR2"),
                              conf_int(H2_group_model, vcov = "CR2")
                            ))

### Game Type ----
H2_game_dat <- H2_data %>% 
  filter(!is.na(Game.Type))
V_H2_game <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H2_game_dat)
H2_game_model <- rma.mv(yi, V_H2_game,
                        mods = ~ Game.Type - 1,
                        random = ~ 1 | No. / Study / effectID, data = H2_game_dat, method = "REML")
summary(H2_game_model)
Wald_test(H2_game_model, constraints = constrain_equal(1:2), vcov = "CR2")
coef_test(H2_game_model, vcov = "CR2")

# saving statistics
H2_moderator_table <- rbind(H2_moderator_table,
                            left_join(
                              coef_test(H2_game_model, vcov = "CR2"),
                              conf_int(H2_game_model, vcov = "CR2")
                            ))

### Sex Ratio ----
H2_sex_dat <- H2_data %>% 
  mutate(male.ratio = Male / (Male + Female))
V_H2_sex <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H2_sex_dat)
H2_sex_model <- rma.mv(yi, V_H2_sex,
                       mods = ~ male.ratio,
                       random = ~ 1 | No. / Study / effectID, data = H2_sex_dat, method = "REML")
summary(H2_sex_model)
coef_test(H2_sex_model, vcov = "CR2")

# saving statistics
H2_moderator_table <- rbind(H2_moderator_table,
                            slice_tail( # sex ratio has intercept
                              left_join(
                                coef_test(H2_sex_model, vcov = "CR2"),
                                conf_int(H2_sex_model, vcov = "CR2")
                              )
                            ))

### Age ----
H2_age_dat <- H2_data %>% 
  filter(!is.na(Age..mean.)) %>%
  mutate(Age..mean. = as.numeric(Age..mean.))
V_H2_age <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H2_age_dat)
H2_age_model <- rma.mv(yi, V_H2_age,
                       mods = ~ Age..mean.,
                       random = ~ 1 | No. / Study / effectID, data = H2_age_dat, method = "REML")
summary(H2_age_model)
coef_test(H2_age_model, vcov = "CR2")
conf_int(H2_age_model, vcov = "CR2")

# saving statistics
H2_moderator_table <- rbind(H2_moderator_table,
                            slice_tail( # age has intercept
                              left_join(
                                coef_test(H2_age_model, vcov = "CR2"),
                                conf_int(H2_age_model, vcov = "CR2")
                              )
                            ))

### relative endowment ----
H2_endow_dat <- H2_data %>% 
  filter(!is.na(Relative.Endowments)) %>%
  mutate(Relative.Endowments = as.numeric(Relative.Endowments))
V_H2_endow <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H2_endow_dat)
H2_endow_model <- rma.mv(yi, V_H2_endow,
                         mods = ~ Relative.Endowments,
                         random = ~ 1 | No. / Study / effectID, data = H2_endow_dat, method = "REML")
summary(H2_endow_model)
coef_test(H2_endow_model, vcov = "CR2")

# saving statistics
H2_moderator_table <- rbind(H2_moderator_table,
                            slice_tail( # Relative endowment has intercept
                              left_join(
                                coef_test(H2_endow_model, vcov = "CR2"),
                                conf_int(H2_endow_model, vcov = "CR2")
                              )
                            ))

### Real Interaction ----
H2_interact_dat <- H2_data %>% 
  filter(!is.na(Real.or.Deception))
V_H2_interact <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H2_interact_dat)
H2_interact_model <- rma.mv(yi, V_H2_interact,
                            mods = ~ Real.or.Deception - 1,
                            random = ~ 1 | No. / Study / effectID, data = H2_interact_dat, method = "REML")
summary(H2_interact_model)
Wald_test(H2_interact_model, constraints = constrain_equal(1:2), vcov = "CR2")
coef_test(H2_interact_model, vcov = "CR2")

# saving statistics
H2_moderator_table <- rbind(H2_moderator_table,
                            left_join(
                              coef_test(H2_interact_model, vcov = "CR2"),
                              conf_int(H2_interact_model, vcov = "CR2")
                            ))

### Student Sample ----
H2_stu_dat <- H2_data %>% 
  filter(!is.na(Student.Sample))
V_H2_stu <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H2_stu_dat)
H2_stu_model <- rma.mv(yi, V_H2_stu,
                       mods = ~ Student.Sample - 1,
                       random = ~ 1 | No. / Study / effectID, data = H2_stu_dat, method = "REML")
summary(H2_stu_model)
Wald_test(H2_stu_model, constraints = constrain_equal(1:2), vcov = "CR2")
coef_test(H2_stu_model, vcov = "CR2")
conf_int(H2_stu_model, vcov = "CR2")

# saving statistics
H2_moderator_table <- rbind(H2_moderator_table,
                            left_join(
                              coef_test(H2_stu_model, vcov = "CR2"),
                              conf_int(H2_stu_model, vcov = "CR2")
                            ))

### Write Table ----
H2_moderator_table <- H2_moderator_table %>% 
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
H2_moderator_table %>% 
  nice_table() %>% 
  flextable::save_as_docx(., path = "tab/H2_moderator.docx")

# Combine H1 with H2
# add prefix before joining:
H1_H2_moderator_table <- full_join(
  H1_moderator_table %>% rename_with(.,
                                     ~ paste0("H1.", .x),
                                     -Moderator)
  ,
  H2_moderator_table %>% rename_with(.,
                                    ~ paste0("H2.", .x),
                                    -Moderator)
  ,
  by = "Moderator")

# Make nice_table
H1_H2_moderator_table %>% 
nice_table(., separate.header = TRUE) %>% 
  flextable::save_as_docx(., path = "tab/H1&2_moderator.docx")


## H5 The Effects of Violator’s Identity on Punishment When Fair moderator analyses ----
H5_moderator_table <- tibble()
### The culture ----
H5_country_dat <- H5_data %>% 
  filter(!is.na(Country)) 
V_H5_country <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H5_country_dat)
H5_country_model <- rma.mv(yi, V_H5_country, 
                           mods = ~ Country - 1,
                           random = ~ 1 | No./ Study / effectID, data = H5_country_dat, method = "REML")
summary(H5_country_model)
Wald_test(H5_country_model, constraints = constrain_equal(1:6), vcov = "CR2", test = "EDT") # HTZ returns NA
coef_test(H5_country_model, vcov = "CR2")
conf_int(H5_country_model, vcov = "CR2")

# saving forest plot to a file
H5_country_col <- case_match(
    H5_country_dat$Country,
    "China"  ~ "#CC6677",
    "Turkey" ~ "#DDCC77",
    "US"     ~ "#999933",
    "UK/US"  ~ "#117733",
    "Argentina" ~ "#44AA99",
    "Switzerland" ~ "#88CCEE",
    "Papua New Guinea" ~ "#332288",
    "Italy" ~ "#882255"
  )
jpeg("fig/H5_country_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H5_country_model, slab = paste(H5_country_dat$Intext.Citation, "S", H5_country_dat$Study, sep = " "),
       colout = H5_country_col)
dev.off()

# saving statistics
H5_moderator_table <- left_join(
  coef_test(H5_country_model, vcov = "CR2"),
  conf_int(H5_country_model, vcov = "CR2")
)

### The cost ----
H5_cost_dat <- H5_data %>% 
  mutate(Punishment.Cost.Ratio = as.numeric(Punishment.Cost.Ratio)) %>% 
  filter(!is.na(Punishment.Cost.Ratio))
V_H5_cost <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H5_cost_dat)
H5_cost_model <- rma.mv(yi, V_H5_cost,
                        mods = ~ Punishment.Cost.Ratio,
                        random = ~ 1 | No. / Study / effectID, data = H5_cost_dat, method = "REML")
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

### Lab or online ----
H5_lab_dat <- H5_data %>% 
  filter(!is.na(Lab.or.Online))
V_H5_lab <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H5_lab_dat)
H5_lab_model <- rma.mv(yi, V_H5_lab,
                       mods = ~ Lab.or.Online - 1,
                       random = ~ 1 | No. / Study / effectID, data = H5_lab_dat, method = "REML")
summary(H5_lab_model)
Wald_test(H5_lab_model, constraints = constrain_equal(1:2), vcov = "CR2")
coef_test(H5_lab_model, vcov = "CR2")
conf_int(H5_lab_model, vcov = "CR2")

# saving statistics
H5_moderator_table <- rbind(H5_moderator_table,
                            left_join(
                              coef_test(H5_lab_model, vcov = "CR2"),
                              conf_int(H5_lab_model, vcov = "CR2")
                            ))

### Group Manipulation ----
H5_group_dat <- H5_data %>% 
  filter(!is.na(Group.Manipulation))
V_H5_group <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H5_group_dat)
H5_group_model <- rma.mv(yi, V_H5_group,
                         mods = ~ Group.Manipulation - 1,
                         random = ~ 1 | No. / Study / effectID, data = H5_group_dat, method = "REML")
summary(H5_group_model)
Wald_test(H5_group_model, constraints = constrain_equal(1:2), vcov = "CR2")
coef_test(H5_group_model, vcov = "CR2")
conf_int(H5_group_model, vcov = "CR2")

# saving forest plot to a file
H5_group_col <- ifelse(H5_group_dat$Group.Manipulation == "Natural", "#004488", "#BB5566")
jpeg("fig/H5_group_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H5_group_model, slab = paste(H5_group_dat$Intext.Citation, "S", H5_group_dat$Study, sep = " "),
       colout = H5_group_col)
dev.off()

# saving statistics
H5_moderator_table <- rbind(H5_moderator_table,
                            left_join(
                              coef_test(H5_group_model, vcov = "CR2"),
                              conf_int(H5_group_model, vcov = "CR2")
                            ))

### Game Type ----
H5_game_dat <- H5_data %>% 
  filter(!is.na(Game.Type))
V_H5_game <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H5_game_dat)
H5_game_model <- rma.mv(yi, V_H5_game,
                        mods = ~ Game.Type - 1,
                        random = ~ 1 | No. / Study / effectID, data = H5_game_dat, method = "REML")
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
V_H5_sex <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H5_sex_dat)
H5_sex_model <- rma.mv(yi, V_H5_sex,
                       mods = ~ male.ratio,
                       random = ~ 1 | No. / Study / effectID, data = H5_sex_dat, method = "REML")
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
  mutate(Age..mean. = as.numeric(Age..mean.)) %>% 
  filter(!is.na(Age..mean.))
V_H5_age <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H5_age_dat)
H5_age_model <- rma.mv(yi, V_H5_age,
                       mods = ~ Age..mean.,
                       random = ~ 1 | No. / Study / effectID, data = H5_age_dat, method = "REML")
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
  mutate(Relative.Endowments = as.numeric(Relative.Endowments)) %>% 
  filter(!is.na(Relative.Endowments))
V_H5_endow <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H5_endow_dat)
H5_endow_model <- rma.mv(yi, V_H5_endow,
                         mods = ~ Relative.Endowments,
                         random = ~ 1 | No. / Study / effectID, data = H5_endow_dat, method = "REML")
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
  filter(!is.na(Real.or.Deception))
V_H5_interact <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H5_interact_dat)
H5_interact_model <- rma.mv(yi, V_H5_interact,
                            mods = ~ Real.or.Deception - 1,
                            random = ~ 1 | No. / Study / effectID, data = H5_interact_dat, method = "REML")
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
  filter(!is.na(Student.Sample))
V_H5_stu <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H5_stu_dat)
H5_stu_model <- rma.mv(yi, V_H5_stu,
                       mods = ~ Student.Sample - 1,
                       random = ~ 1 | No. / Study / effectID, data = H5_stu_dat, method = "REML")
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



## H6 The Effects of Victim’s Identity on Punishment When Fair moderator analyses ----
H6_moderator_table <- tibble()
### The culture ----
H6_country_dat <- H6_data %>% 
  filter(!is.na(Country)) 
V_H6_country <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H6_country_dat)
H6_country_model <- rma.mv(yi, vi, # V cannot convergence
                           mods = ~ Country - 1,
                           random = ~ 1 | No./ Study / effectID, data = H6_country_dat, method = "REML")
summary(H6_country_model)
# Wald_test(H6_country_model, constraints = constrain_equal(1:5), vcov = "CR2", test = "chi-sq") # V is singular
coef_test(H6_country_model, vcov = "CR2")
conf_int(H6_country_model, vcov = "CR2")

# saving statistics
H6_moderator_table <- left_join(
  coef_test(H6_country_model, vcov = "CR2"),
  conf_int(H6_country_model, vcov = "CR2")
)

### The cost ----
H6_cost_dat <- H6_data %>% 
  mutate(Punishment.Cost.Ratio = as.numeric(Punishment.Cost.Ratio)) %>% 
  filter(!is.na(Punishment.Cost.Ratio))
V_H6_cost <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H6_cost_dat)
H6_cost_model <- rma.mv(yi, V_H6_cost,
                        mods = ~ Punishment.Cost.Ratio,
                        random = ~ 1 | No. / Study / effectID, data = H6_cost_dat, method = "REML")
summary(H6_cost_model)
coef_test(H6_cost_model, vcov = "CR2")
conf_int(H6_cost_model, vcov = "CR2")

# saving statistics
H6_moderator_table <- rbind(H6_moderator_table,
                            slice_tail( # cost ratio has intercept
                              left_join(
                                coef_test(H6_cost_model, vcov = "CR2"),
                                conf_int(H6_cost_model, vcov = "CR2")
                              )
                            ))

### Lab or online ----
H6_lab_dat <- H6_data %>% 
  filter(!is.na(Lab.or.Online))
V_H6_lab <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H6_lab_dat)
H6_lab_model <- rma.mv(yi, V_H6_lab,
                       mods = ~ Lab.or.Online - 1,
                       random = ~ 1 | No. / Study / effectID, data = H6_lab_dat, method = "REML",
                       control=list(rel.tol=1e-9))
summary(H6_lab_model)
Wald_test(H6_lab_model, constraints = constrain_equal(1:2), vcov = "CR2")
coef_test(H6_lab_model, vcov = "CR2")
conf_int(H6_lab_model, vcov = "CR2")

# saving statistics
H6_moderator_table <- rbind(H6_moderator_table,
                            left_join(
                              coef_test(H6_lab_model, vcov = "CR2"),
                              conf_int(H6_lab_model, vcov = "CR2")
                            ))

### Group Manipulation ----
# Not appliable, all data is Natural group
# H6_group_dat <- H6_data %>%
#   filter(!is.na(Group.Manipulation))
# V_H6_group <- vcalc(H6_group_dat$vi, cluster=H6_group_dat$Sample.ID, rho = 0.5)
# H6_group_model <- rma.mv(yi, V_H6_group,
#                          mods = ~ Group.Manipulation,
#                          random = ~ 1 | No. / Study / effectID, data = H6_group_dat, method = "REML")
# summary(H6_group_model)

### Game Type ----
H6_game_dat <- H6_data %>% 
  filter(!is.na(Game.Type))
V_H6_game <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H6_game_dat)
H6_game_model <- rma.mv(yi, vi, # V cannot convergence
                        mods = ~ Game.Type - 1,
                        random = ~ 1 | No. / Study / effectID, data = H6_game_dat, method = "REML")
summary(H6_game_model)
Wald_test(H6_game_model, constraints = constrain_equal(1:2), vcov = "CR2")
coef_test(H6_game_model, vcov = "CR2")

# saving statistics
H6_moderator_table <- rbind(H6_moderator_table,
                            left_join(
                              coef_test(H6_game_model, vcov = "CR2"),
                              conf_int(H6_game_model, vcov = "CR2")
                            ))

### Sex Ratio ----
H6_sex_dat <- H6_data %>% 
  mutate(male.ratio = Male / (Male + Female))
V_H6_sex <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H6_sex_dat)
H6_sex_model <- rma.mv(yi, vi, # V cannot convergence
                       mods = ~ male.ratio,
                       random = ~ 1 | No. / Study / effectID, data = H6_sex_dat, method = "REML")
summary(H6_sex_model)
coef_test(H6_sex_model, vcov = "CR2")

# saving statistics
H6_moderator_table <- rbind(H6_moderator_table,
                            slice_tail( # sex ratio has intercept
                              left_join(
                                coef_test(H6_sex_model, vcov = "CR2"),
                                conf_int(H6_sex_model, vcov = "CR2")
                              )
                            ))

### Age ----
H6_age_dat <- H6_data %>% 
  filter(!is.na(Age..mean.)) %>%
  mutate(Age..mean. = as.numeric(Age..mean.))
V_H6_age <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H6_age_dat)
H6_age_model <- rma.mv(yi, V_H6_age,
                       mods = ~ Age..mean.,
                       random = ~ 1 | No. / Study / effectID, data = H6_age_dat, method = "REML")
summary(H6_age_model)
coef_test(H6_age_model, vcov = "CR2")
conf_int(H6_age_model, vcov = "CR2")

# saving statistics
H6_moderator_table <- rbind(H6_moderator_table,
                            slice_tail( # age has intercept
                              left_join(
                                coef_test(H6_age_model, vcov = "CR2"),
                                conf_int(H6_age_model, vcov = "CR2")
                              )
                            ))

### relative endowment ----
H6_endow_dat <- H6_data %>% 
  filter(!is.na(Relative.Endowments)) %>%
  mutate(Relative.Endowments = as.numeric(Relative.Endowments))
V_H6_endow <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H6_endow_dat)
H6_endow_model <- rma.mv(yi, V_H6_endow,
                         mods = ~ Relative.Endowments,
                         random = ~ 1 | No. / Study / effectID, data = H6_endow_dat, method = "REML")
summary(H6_endow_model)
coef_test(H6_endow_model, vcov = "CR2")

# saving statistics
H6_moderator_table <- rbind(H6_moderator_table,
                            slice_tail( # Relative endowment has intercept
                              left_join(
                                coef_test(H6_endow_model, vcov = "CR2"),
                                conf_int(H6_endow_model, vcov = "CR2")
                              )
                            ))

### Real Interaction ----
H6_interact_dat <- H6_data %>% 
  filter(!is.na(Real.or.Deception))
V_H6_interact <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H6_interact_dat)
H6_interact_model <- rma.mv(yi, vi, # V cannot convergence
                            mods = ~ Real.or.Deception - 1,
                            random = ~ 1 | No. / Study / effectID, data = H6_interact_dat, method = "REML")
summary(H6_interact_model)
Wald_test(H6_interact_model, constraints = constrain_equal(1:2), vcov = "CR2")
coef_test(H6_interact_model, vcov = "CR2")

# saving statistics
H6_moderator_table <- rbind(H6_moderator_table,
                            left_join(
                              coef_test(H6_interact_model, vcov = "CR2"),
                              conf_int(H6_interact_model, vcov = "CR2")
                            ))

### Student Sample ----
H6_stu_dat <- H6_data %>% 
  filter(!is.na(Student.Sample))
V_H6_stu <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = 0.5, data = H6_stu_dat)
H6_stu_model <- rma.mv(yi, V_H6_stu,
                       mods = ~ Student.Sample - 1,
                       random = ~ 1 | No. / Study / effectID, data = H6_stu_dat, method = "REML")
summary(H6_stu_model)
Wald_test(H6_stu_model, constraints = constrain_equal(1:2), vcov = "CR2")
coef_test(H6_stu_model, vcov = "CR2")
conf_int(H6_stu_model, vcov = "CR2")

# saving statistics
H6_moderator_table <- rbind(H6_moderator_table,
                            left_join(
                              coef_test(H6_stu_model, vcov = "CR2"),
                              conf_int(H6_stu_model, vcov = "CR2")
                            ))

### Write Table ----
H6_moderator_table <- H6_moderator_table %>% 
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
H6_moderator_table %>% 
  nice_table() %>% 
  flextable::save_as_docx(., path = "tab/H6_moderator.docx")

# Combine H5 with H6
# add prefix before joining:
H5_H6_moderator_table <- full_join(
  H5_moderator_table %>% rename_with(.,
                                     ~ paste0("H5.", .x),
                                     -Moderator)
  ,
  H6_moderator_table %>% rename_with(.,
                                     ~ paste0("H6.", .x),
                                     -Moderator)
  ,
  by = "Moderator")

# Make nice_table
H5_H6_moderator_table %>% 
  nice_table(., separate.header = TRUE) %>% 
  flextable::save_as_docx(., path = "tab/H5&6_moderator.docx")



# Sensitivity analysis for the value of correlation of effect sizes obtained on the same group of participants ########################
H1_sensitivity <- tibble()
H2_sensitivity <- tibble()
H5_sensitivity <- tibble()
H6_sensitivity <- tibble()
# H1_data <- H1_data %>% 
#   mutate(StudyID = paste0(No.,"_",Study))
for (r in 1:10) {
  R = r/10 # real correlation coefficient
  ## H1
  V_H1_s <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = R, data = H1_data)
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
  ## H2
  V_H2_s <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = R, data = H2_data)
  H2_model_s <- rma.mv(yi, V_H2_s, random = ~ 1 | No. / Study / effectID, data = H2_data, method = "REML",
                       control=list(rel.tol=1e-9))
  
  # saving statistics
  H2_sensitivity <- rbind(H2_sensitivity,
                          left_join(
                            coef_test(H2_model_s, vcov = "CR2"),
                            conf_int(H2_model_s, vcov = "CR2")
                          )
  )
  ## H5
  V_H5_s <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = R, data = H5_data)
  H5_model_s <- rma.mv(yi, V_H5_s, random = ~ 1 | No. / Study / effectID, data = H5_data, method = "REML",
                       control=list(rel.tol=1e-9))
  
  # saving statistics
  H5_sensitivity <- rbind(H5_sensitivity,
                          left_join(
                            coef_test(H5_model_s, vcov = "CR2"),
                            conf_int(H5_model_s, vcov = "CR2")
                          )
  )
  ## H6
  V_H6_s <- vcalc(vi, cluster = Sample.ID, obs = effectID, rho = R, data = H6_data)
  H6_model_s <- rma.mv(yi, V_H6_s, random = ~ 1 | No. / Study / effectID, data = H6_data, method = "REML",
                       control=list(rel.tol=1e-9))
  
  # saving statistics
  H6_sensitivity <- rbind(H6_sensitivity,
                          left_join(
                            coef_test(H6_model_s, vcov = "CR2"),
                            conf_int(H6_model_s, vcov = "CR2")
                          )
  )
}

## write table ----
H1_sensitivity[,1] <- seq(from = 0.1, to = 1, by = 0.1)
H2_sensitivity[,1] <- seq(from = 0.1, to = 1, by = 0.1)
H5_sensitivity[,1] <- seq(from = 0.1, to = 1, by = 0.1)
H6_sensitivity[,1] <- seq(from = 0.1, to = 1, by = 0.1)

H1_sensitivity_short <- H1_sensitivity %>% 
  select(Coef, beta, p_Satt) %>% 
  rename(
    r = "Coef",
    H1.g = "beta",
    H1.p = "p_Satt"
  )
H2_sensitivity_short <- H2_sensitivity %>% 
  select(beta, p_Satt) %>% 
  rename(
    H2.g = "beta",
    H2.p = "p_Satt"
  )
H5_sensitivity_short <- H5_sensitivity %>% 
  select(beta, p_Satt) %>% 
  rename(
    H5.g = "beta",
    H5.p = "p_Satt"
  )
H6_sensitivity_short <- H6_sensitivity %>% 
  select(beta, p_Satt) %>% 
  rename(
    H6.g = "beta",
    H6.p = "p_Satt"
  )

cbind(H1_sensitivity_short, H2_sensitivity_short, H5_sensitivity_short, H6_sensitivity_short) %>% 
  nice_table(., separate.header = TRUE) %>% 
  flextable::save_as_docx(., path = "tab/sensitivity.docx")
