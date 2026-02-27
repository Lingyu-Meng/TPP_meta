required_packages <- c("tidyverse", "metafor", "clubSandwich", "MAd", "weightr", "puniform", "metapower")

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
    No., Study, `Person A affiliation`, `Recipient/Player B Affliation`, `Sample ID`,
    `Fairness/Cooperativeness of Person A`, `Fairness/Cooperativeness of Person B`,
    `P_Average Invested Amount`, `P_SD invested amount`, `N by condition`, `Total Trials`, `Intext-Citation`,
    Country, `Punishment-Cost Ratio`, `Lab or Online`, `Group Manipulation`, Female, Male, `Game Type`
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
    `Between or Within?...37` == "Between",
    `Fairness/Cooperativeness A Binary` == 0, # unfair condition
    `Person A affiliation` %in% c("Ingroup", "Outgroup")
  ) %>%
  select(
    No., Study, `Person A affiliation`, `Recipient/Player B Affliation`,
    `Fairness/Cooperativeness of Person A`, `Fairness/Cooperativeness of Person B`,
    `P_Average Invested Amount`, `P_SD invested amount`, `N by condition`, `Total Trials`, `Intext-Citation`,
    Country, `Punishment-Cost Ratio`, `Lab or Online`, `Group Manipulation`, Female, Male, `Game Type`
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

V_H1 <- vcalc(H1_data$vi, cluster=H1_data$Sample.ID, rho = 0.5)
H1_model <- rma.mv(yi, V_H1, random = ~ 1 | No. / Study / effectID, data = H1_data, method = "REML")
summary(H1_model)

# checking likelihood profile plots
profile(H1_model)

# obtaining robust variance estimates
coef_test(H1_model, vcov = "CR2")
conf_int(H1_model, vcov = "CR2")

# calculating the overall amount of heterogeneity
sum(H1_model$sigma2)

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
metafor::funnel(H1_model, refline=0, level=c(95, 99), 
                shade=c("white", "grey"))

## H2 (Deterrence perspective): ----
H2_data <- MA_dat %>%
  filter(
    `Fairness/Cooperativeness A Binary` == 0, # unfair condition
    `Recipient/Player B Affliation` %in% c("Ingroup", "Outgroup")
  ) %>%
  # group_by(No., Study, `Recipient/Player B Affliation`) %>%
  select(
    No., Study, `Person A affiliation`, `Recipient/Player B Affliation`, `Sample ID`,
    `Fairness/Cooperativeness of Person A`, `Fairness/Cooperativeness of Person B`,
    `P_Average Invested Amount`, `P_SD invested amount`, `N by condition`, `Total Trials`, `Intext-Citation`,
    Country, `Punishment-Cost Ratio`, `Lab or Online`, `Group Manipulation`, Female, Male, `Game Type`
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
    `Between or Within?...39` == "Between",
    `Fairness/Cooperativeness A Binary` == 0, # unfair condition
    `Recipient/Player B Affliation` %in% c("Ingroup", "Outgroup")
  ) %>%
  select(
    No., Study, `Person A affiliation`, `Recipient/Player B Affliation`,
    `Fairness/Cooperativeness of Person A`, `Fairness/Cooperativeness of Person B`,
    `P_Average Invested Amount`, `P_SD invested amount`, `N by condition`, `Total Trials`, `Intext-Citation`,
    Country, `Punishment-Cost Ratio`, `Lab or Online`, `Group Manipulation`, Female, Male, `Game Type`
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

V_H2 <- vcalc(H2_data$vi, cluster=H2_data$Sample.ID, rho = 0.5)
H2_model <- rma.mv(yi, V_H2, random = ~ 1 | No. / Study / effectID, data = H2_data, method = "REML")
summary(H2_model)

# checking likelihood profile plots
profile(H2_model)

# obtaining robust variance estimates
coef_test(H2_model, vcov = "CR2")
conf_int(H2_model, vcov = "CR2")

# calculating the overall amount of heterogeneity
sum(H2_model$sigma2)

# forest plot
H2_f <- forest(H2_model, slab = paste(H2_data$Intext.Citation, "S", H2_data$Study, sep = " "))
# saving forest plot to a file
jpeg("fig/H2_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H2_model, slab = paste(H2_data$Intext.Citation, "S", H2_data$Study, sep = " "))
dev.off()

# funnel plot
metafor::funnel(H2_model)

# contour-enhanced funnel plot
metafor::funnel(H2_model, refline=0, level=c(95, 99), 
                shade=c("white", "grey"))

## H3 Interaction Effect between Violator’s Identity and Victim Salience ----
### The Mantel-Haenszel method ----
MH_dat <- MA_dat %>%
  group_by(No., Study) %>%
  # keep only studies that have all four A×B combinations
  filter(
    `Fairness/Cooperativeness A Binary` == 0,
    n_distinct(`Person A affiliation`) >= 2,
    n_distinct(`Recipient/Player B Affliation`) >= 2,
    all(expand.grid(
      `Person A affiliation` = c("Ingroup", "Outgroup"),
      `Recipient/Player B Affliation` = c("Ingroup", "Outgroup")
    ) %>%
      apply(1, paste, collapse = "_") %in%
      paste(`Person A affiliation`, `Recipient/Player B Affliation`, sep = "_"))
  ) %>%
  ungroup()
# 
# I_I <- MH_dat %>%
#   filter(
#     `Person A affiliation` == "Ingroup",
#     `Recipient/Player B Affliation` == "Ingroup"
#   )
# I_O <- MH_dat %>%
#   filter(
#     `Person A affiliation` == "Ingroup",
#     `Recipient/Player B Affliation` == "Outgroup"
#   )
# O_I <- MH_dat %>%
#   filter(
#     `Person A affiliation` == "Outgroup",
#     `Recipient/Player B Affliation` == "Ingroup"
#   )
# O_O <- MH_dat %>%
#   filter(
#     `Person A affiliation` == "Outgroup",
#     `Recipient/Player B Affliation` == "Outgroup"
#   )
# n1i <- I_I$`N by condition` + I_O$`N by condition`
# n2i <- O_I$`N by condition` + O_O$`N by condition`
# 
# MH_model <- rma.mh(I_I$P_Frequency, #  I O
#   I_O$P_Frequency, # I  a b
#   O_I$P_Frequency, # O  c d
#   O_O$P_Frequency,
#   n1i,
#   n2i,
#   measure = "RR"
# )
# summary(MH_model)
# # k = 8 studies included
# MH_model$beta
# forest(MH_model)
# # RR = (a/(a+b)) / (c/(c+d))
# # when person A is an ingroup member punishment probablity

### The moderator method ----
H3_model <- H1_data %>%
  rma.mv(yi, V_H1,
  mods = ~ Recipient.Player.B.Affliation,
  random = ~ 1 | No. / Study / effectID,
  data = .,
  method = "REML",
)
summary(H3_model)
# asign color for ingroup, outgroup, and neutral
H3_col <- ifelse(H1_data$Recipient.Player.B.Affliation == "Ingroup", "blue",
                 ifelse(H1_data$Recipient.Player.B.Affliation == "Outgroup", "red", "gray"))
jpeg("fig/H3_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H3_model, slab = paste(H1_data$Intext.Citation, "S", H1_data$Study, sep = " "),
       colout = H3_col)
dev.off()

#### Robustness check 1: Exclude studies with neutral recipient/player B affiliation ----
H3_data <- H1_data %>%
  filter(Recipient.Player.B.Affliation != "Neutral")

V_H3 <- vcalc(H3_data$vi, cluster=H3_data$Sample.ID, rho = 0.5)

H3_model2 <-H3_data  %>%
  rma.mv(yi, V_H3,
         mods = ~ Recipient.Player.B.Affliation,
         random = ~ 1 | No. / Study / effectID,
         data = .,
         method = "REML",
)
summary(H3_model)

#### Robustness check 2: Using the H2_data ----
H3_model3 <- H2_data %>%
  rma.mv(yi, V_H2,
         mods = ~ Person.A.affiliation,
         random = ~ 1 | No. / Study / effectID,
         data = .,
         method = "REML",
)
summary(H3_model3)

### The regression method ----
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
  )

H3_model_r <- H3_data_r %>%
  rma.mv(yi = `P_Average Invested Amount_Standardized`, V = P_SV,
  mods = ~ `Person A` * `Person B`,
  random = ~ 1 | No.,
  data = .,
  method = "REML",
)
summary(H3_model_r)

# Fair Condition Analysis ################
## H5 The Effects of Violator’s Identity on Punishment ----
H5_data <- MA_dat %>%
  filter(
    `Fairness/Cooperativeness A Binary` == 1, # fair condition
    `Person A affiliation` %in% c("Ingroup", "Outgroup")
  ) %>%
  select(
    No., Study, `Person A affiliation`, `Recipient/Player B Affliation`, `Sample ID`,
    `Fairness/Cooperativeness of Person A`, `Fairness/Cooperativeness of Person B`,
    `P_Average Invested Amount`, `P_SD invested amount`, `N by condition`, `Total Trials`, `Intext-Citation`,
    Country, `Punishment-Cost Ratio`, `Lab or Online`, `Group Manipulation`, Female, Male, `Game Type`
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
    `Between or Within?...37` == "Between",
    `Fairness/Cooperativeness A Binary` == 1, # fair condition
    `Person A affiliation` %in% c("Ingroup", "Outgroup")
  ) %>%
  select(
    No., Study, `Person A affiliation`, `Recipient/Player B Affliation`,
    `Fairness/Cooperativeness of Person A`, `Fairness/Cooperativeness of Person B`,
    `P_Average Invested Amount`, `P_SD invested amount`, `N by condition`, `Total Trials`, `Intext-Citation`,
    Country, `Punishment-Cost Ratio`, `Lab or Online`, `Group Manipulation`, Female, Male, `Game Type`
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

V_H5 <- vcalc(H5_data$vi, cluster=H5_data$Sample.ID, rho = 0.5)
H5_model <- rma.mv(yi, V_H5, random = ~ 1 | No. / Study / effectID, data = H5_data, method = "REML")
summary(H5_model)

# forest plot
H5_f <- forest(H5_model, slab = paste(H5_data$Intext.Citation, "S", H5_data$Study, sep = " "))
# saving forest plot to a file
jpeg("fig/H5_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H5_model, slab = paste(H5_data$Intext.Citation, "S", H5_data$Study, sep = " "))
dev.off()

# funnel plot
metafor::funnel(H5_model)

# contour-enhanced funnel plot
metafor::funnel(H5_model, refline=0, level=c(95, 99), 
                shade=c("white", "grey"))

## H6 The Effects of Victim’s Identity on Punishment ----
H6_data <- MA_dat %>%
  filter(
    `Fairness/Cooperativeness A Binary` == 1, # fair condition
    `Recipient/Player B Affliation` %in% c("Ingroup", "Outgroup")
  ) %>%
  # group_by(No., Study, `Recipient/Player B Affliation`) %>%
  select(
    No., Study, `Person A affiliation`, `Recipient/Player B Affliation`, `Sample ID`,
    `Fairness/Cooperativeness of Person A`, `Fairness/Cooperativeness of Person B`,
    `P_Average Invested Amount`, `P_SD invested amount`, `N by condition`, `Total Trials`, `Intext-Citation`,
    Country, `Punishment-Cost Ratio`, `Lab or Online`, `Group Manipulation`, Female, Male, `Game Type`
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
    `Between or Within?...39` == "Between",
    `Fairness/Cooperativeness A Binary` == 1, # fair condition
    `Recipient/Player B Affliation` %in% c("Ingroup", "Outgroup")
  ) %>%
  select(
    No., Study, `Person A affiliation`, `Recipient/Player B Affliation`,
    `Fairness/Cooperativeness of Person A`, `Fairness/Cooperativeness of Person B`,
    `P_Average Invested Amount`, `P_SD invested amount`, `N by condition`, `Total Trials`, `Intext-Citation`,
    Country, `Punishment-Cost Ratio`, `Lab or Online`, `Group Manipulation`, Female, Male, `Game Type`
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

V_H6 <- vcalc(H6_data$vi, cluster=H6_data$Sample.ID, rho = 0.5)
H6_model <- rma.mv(yi, V_H6, random = ~ 1 | No. / Study / effectID, data = H6_data, method = "REML")
summary(H6_model)

# forest plot
H6_f <- forest(H6_model, slab = paste(H6_data$Intext.Citation, "S", H6_data$Study, sep = " "))
# saving forest plot to a file
jpeg("fig/H6_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H6_model, slab = paste(H6_data$Intext.Citation, "S", H6_data$Study, sep = " "))
dev.off()

# funnel plot
metafor::funnel(H6_model)

# contour-enhanced funnel plot
metafor::funnel(H6_model, refline=0, level=c(95, 99), 
                shade=c("white", "grey"))

# Single moderator analyses ################
## H1 The Effects of Violator’s Identity on Punishment moderator analyses ----
### The culture ----
H1_country_dat <- H1_data %>% 
  filter(!is.na(Country))
V_H1_country <- vcalc(H1_country_dat$vi, cluster=H1_country_dat$Sample.ID, rho = 0.5)
H1_country_model <- rma.mv(yi, V_H1_country,
                           mods = ~ Country,
                           random = ~ 1 | No. / Study / effectID, data = H1_country_dat, method = "REML")
summary(H1_country_model)

### The cost ----
H1_cost_dat <- H1_data %>% 
  filter(!is.na(Punishment.Cost.Ratio)) %>% 
  mutate(Punishment.Cost.Ratio = as.numeric(Punishment.Cost.Ratio))
V_H1_cost <- vcalc(H1_cost_dat$vi, cluster=H1_cost_dat$Sample.ID, rho = 0.5)
H1_cost_model <- rma.mv(yi, V_H1_cost,
                           mods = ~ Punishment.Cost.Ratio,
                           random = ~ 1 | No. / Study / effectID, data = H1_cost_dat, method = "REML")
summary(H1_cost_model)

### Lab or online ----
H1_lab_dat <- H1_data %>% 
  filter(!is.na(Punishment.Cost.Ratio))
V_H1_lab <- vcalc(H1_lab_dat$vi, cluster=H1_lab_dat$Sample.ID, rho = 0.5)
H1_lab_model <- rma.mv(yi, V_H1_lab,
                        mods = ~ Lab.or.Online,
                        random = ~ 1 | No. / Study / effectID, data = H1_lab_dat, method = "REML")
summary(H1_lab_model)

# saving forest plot to a file
H1_lab_col <- ifelse(H1_lab_dat$Lab.or.Online == "Online", "blue", "red")
jpeg("fig/H1_lab_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H1_lab_model, slab = paste(H1_lab_dat$Intext.Citation, "S", H1_lab_dat$Study, sep = " "),
       colout = H1_lab_col)
dev.off()

### Group Manipulation ----
H1_group_dat <- H1_data %>% 
  filter(!is.na(Group.Manipulation))
V_H1_group <- vcalc(H1_group_dat$vi, cluster=H1_group_dat$Sample.ID, rho = 0.5)
H1_group_model <- rma.mv(yi, V_H1_group,
                       mods = ~ Group.Manipulation,
                       random = ~ 1 | No. / Study / effectID, data = H1_group_dat, method = "REML")
summary(H1_group_model)

# saving forest plot to a file
H1_group_col <- ifelse(H1_group_dat$Group.Manipulation == "Natural", "blue", "red")
jpeg("fig/H1_group_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H1_group_model, slab = paste(H1_group_dat$Intext.Citation, "S", H1_group_dat$Study, sep = " "),
       colout = H1_group_col)
dev.off()

### Game Type ----
H1_game_dat <- H1_data %>% 
  filter(!is.na(Game.Type))
V_H1_game <- vcalc(H1_game_dat$vi, cluster=H1_game_dat$Sample.ID, rho = 0.5)
H1_game_model <- rma.mv(yi, V_H1_game,
                         mods = ~ Game.Type,
                         random = ~ 1 | No. / Study / effectID, data = H1_game_dat, method = "REML")
summary(H1_game_model)

### Sex Ratio ----
H1_sex_dat <- H1_data %>% 
  mutate(male.ratio = Male / (Male + Female))
V_H1_sex <- vcalc(H1_sex_dat$vi, cluster=H1_sex_dat$Sample.ID, rho = 0.5)
H1_sex_model <- rma.mv(yi, V_H1_sex,
                        mods = ~ male.ratio,
                        random = ~ 1 | No. / Study / effectID, data = H1_sex_dat, method = "REML")
summary(H1_sex_model)

## H2 (Deterrence perspective): ----
### The culture ----
H2_country_dat <- H2_data %>% 
  filter(!is.na(Country))
V_H2_country <- vcalc(H2_country_dat$vi, cluster=H2_country_dat$Sample.ID, rho = 0.5)
H2_country_model <- rma.mv(yi, V_H2_country,
                           mods = ~ Country,
                           random = ~ 1 | No. / Study / effectID, data = H2_country_dat, method = "REML")
summary(H2_country_model)

### The cost ----
H2_cost_dat <- H2_data %>% 
  filter(!is.na(Punishment.Cost.Ratio)) %>% 
  mutate(Punishment.Cost.Ratio = as.numeric(Punishment.Cost.Ratio))
V_H2_cost <- vcalc(H2_cost_dat$vi, cluster=H2_cost_dat$Sample.ID, rho = 0.5)
H2_cost_model <- rma.mv(yi, V_H2_cost,
                        mods = ~ Punishment.Cost.Ratio,
                        random = ~ 1 | No. / Study / effectID, data = H2_cost_dat, method = "REML")
summary(H2_cost_model)

# plotting the cost-punishment relationship
plot(H2_cost_dat$Punishment.Cost.Ratio, H2_cost_dat$yi, 
     xlab="Punishment Cost Ratio", ylab="Effect Size", main="Scatter plot")

### Lab or online ----
H2_lab_dat <- H2_data %>% 
  filter(!is.na(Punishment.Cost.Ratio)) %>% 
  mutate(Punishment.Cost.Ratio = as.numeric(Punishment.Cost.Ratio))
V_H2_lab <- vcalc(H2_lab_dat$vi, cluster=H2_lab_dat$Sample.ID, rho = 0.5)
H2_lab_model <- rma.mv(yi, V_H2_lab,
                       mods = ~ Lab.or.Online,
                       random = ~ 1 | No. / Study / effectID, data = H2_lab_dat, method = "REML")
summary(H2_lab_model)

### Group Manipulation ----
H2_group_dat <- H2_data %>% 
  filter(!is.na(Group.Manipulation))
V_H2_group <- vcalc(H2_group_dat$vi, cluster=H2_group_dat$Sample.ID, rho = 0.5)
H2_group_model <- rma.mv(yi, V_H2_group,
                         mods = ~ Group.Manipulation,
                         random = ~ 1 | No. / Study / effectID, data = H2_group_dat, method = "REML")
summary(H2_group_model)

### Game Type ----
H2_game_dat <- H2_data %>% 
  filter(!is.na(Game.Type))
V_H2_game <- vcalc(H2_game_dat$vi, cluster=H2_game_dat$Sample.ID, rho = 0.5)
H2_game_model <- rma.mv(yi, V_H2_game,
                        mods = ~ Game.Type,
                        random = ~ 1 | No. / Study / effectID, data = H2_game_dat, method = "REML")
summary(H2_game_model)

### Sex Ratio ----
H2_sex_dat <- H2_data %>% 
  mutate(male.ratio = Male / (Male + Female))
V_H2_sex <- vcalc(H2_sex_dat$vi, cluster=H2_sex_dat$Sample.ID, rho = 0.5)
H2_sex_model <- rma.mv(yi, V_H2_sex,
                       mods = ~ male.ratio,
                       random = ~ 1 | No. / Study / effectID, data = H2_sex_dat, method = "REML")
summary(H2_sex_model)

## H5 The Effects of Violator’s Identity on Punishment moderator analyses ----
### The culture ----
H5_country_dat <- H5_data %>% 
  filter(!is.na(Country)) 
# V_H5_country <- vcalc(H5_country_dat$vi, cluster=H5_country_dat$Sample.ID, rho = 0.5)
H5_country_model <- rma.mv(yi, vi, # correlated V is close to singular
                           mods = ~ Country,
                           random = ~ 1 | No./ Study / effectID, data = H5_country_dat, method = "REML")
summary(H5_country_model)

# saving forest plot to a file
H5_country_col <- ifelse(H5_country_dat$Country == "China", "red",
                         ifelse(H5_country_dat$Country == "Turkey", "green",
                                ifelse(H5_country_dat$Country == "USA", "blue", "purple"))
                         )
jpeg("fig/H5_country_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H5_country_model, slab = paste(H5_country_dat$Intext.Citation, "S", H5_country_dat$Study, sep = " "),
       colout = H5_country_col)
dev.off()

### The cost ----
H5_cost_dat <- H5_data %>% 
  filter(!is.na(Punishment.Cost.Ratio)) %>% 
  mutate(Punishment.Cost.Ratio = as.numeric(Punishment.Cost.Ratio))
V_H5_cost <- vcalc(H5_cost_dat$vi, cluster=H5_cost_dat$Sample.ID, rho = 0.5)
H5_cost_model <- rma.mv(yi, V_H5_cost,
                        mods = ~ Punishment.Cost.Ratio,
                        random = ~ 1 | No. / Study / effectID, data = H5_cost_dat, method = "REML")
summary(H5_cost_model)

### Lab or online ----
H5_lab_dat <- H5_data %>% 
  filter(!is.na(Punishment.Cost.Ratio))
V_H5_lab <- vcalc(H5_lab_dat$vi, cluster=H5_lab_dat$Sample.ID, rho = 0.5)
H5_lab_model <- rma.mv(yi, V_H5_lab,
                       mods = ~ Lab.or.Online,
                       random = ~ 1 | No. / Study / effectID, data = H5_lab_dat, method = "REML")
summary(H5_lab_model)

### Group Manipulation ----
H5_group_dat <- H5_data %>% 
  filter(!is.na(Group.Manipulation))
V_H5_group <- vcalc(H5_group_dat$vi, cluster=H5_group_dat$Sample.ID, rho = 0.5)
H5_group_model <- rma.mv(yi, V_H5_group,
                         mods = ~ Group.Manipulation,
                         random = ~ 1 | No. / Study / effectID, data = H5_group_dat, method = "REML")
summary(H5_group_model)

# saving forest plot to a file
H5_group_col <- ifelse(H5_group_dat$Group.Manipulation == "Natural", "blue", "red")
jpeg("fig/H5_group_Forest.jpeg", height = 12, width = 8, units = "in", res = 300)
forest(H5_group_model, slab = paste(H5_group_dat$Intext.Citation, "S", H5_group_dat$Study, sep = " "),
       colout = H5_group_col)
dev.off()

### Game Type ----
H5_game_dat <- H5_data %>% 
  filter(!is.na(Game.Type))
V_H5_game <- vcalc(H5_game_dat$vi, cluster=H5_game_dat$Sample.ID, rho = 0.5)
H5_game_model <- rma.mv(yi, V_H5_game,
                        mods = ~ Game.Type,
                        random = ~ 1 | No. / Study / effectID, data = H5_game_dat, method = "REML")
summary(H5_game_model)

### Sex Ratio ----
H5_sex_dat <- H5_data %>% 
  mutate(male.ratio = Male / (Male + Female))
V_H5_sex <- vcalc(H5_sex_dat$vi, cluster=H5_sex_dat$Sample.ID, rho = 0.5)
H5_sex_model <- rma.mv(yi, V_H5_sex,
                       mods = ~ male.ratio,
                       random = ~ 1 | No. / Study / effectID, data = H5_sex_dat, method = "REML")
summary(H5_sex_model)

## H6 The Effects of Victim’s Identity on Punishment moderator analyses ----
### The culture ----
H6_country_dat <- H6_data %>% 
  filter(!is.na(Country)) 
V_H6_country <- vcalc(H6_country_dat$vi, cluster=H6_country_dat$Sample.ID, rho = 0.5)
H6_country_model <- rma.mv(yi, V_H6_country,
                           mods = ~ Country,
                           random = ~ 1 | No./ Study / effectID, data = H6_country_dat, method = "REML")
summary(H6_country_model)

### The cost ----
H6_cost_dat <- H6_data %>% 
  filter(!is.na(Punishment.Cost.Ratio)) %>% 
  mutate(Punishment.Cost.Ratio = as.numeric(Punishment.Cost.Ratio))
V_H6_cost <- vcalc(H6_cost_dat$vi, cluster=H6_cost_dat$Sample.ID, rho = 0.5)
H6_cost_model <- rma.mv(yi, V_H6_cost,
                        mods = ~ Punishment.Cost.Ratio,
                        random = ~ 1 | No. / Study / effectID, data = H6_cost_dat, method = "REML")
summary(H6_cost_model)

### Lab or online ----
H6_lab_dat <- H6_data %>% 
  filter(!is.na(Punishment.Cost.Ratio))
V_H6_lab <- vcalc(H6_lab_dat$vi, cluster=H6_lab_dat$Sample.ID, rho = 0.5)
H6_lab_model <- rma.mv(yi, V_H6_lab,
                       mods = ~ Lab.or.Online,
                       random = ~ 1 | No. / Study / effectID, data = H6_lab_dat, method = "REML")
summary(H6_lab_model)

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
# V_H6_game <- vcalc(H6_game_dat$vi, cluster=H6_game_dat$Sample.ID, rho = 0.5)
H6_game_model <- rma.mv(yi, vi, # V is close to singular
                        mods = ~ Game.Type,
                        random = ~ 1 | No. / Study / effectID, data = H6_game_dat, method = "REML")
summary(H6_game_model)

### Sex Ratio ----
H6_sex_dat <- H6_data %>% 
  mutate(male.ratio = Male / (Male + Female))
V_H6_sex <- vcalc(H6_sex_dat$vi, cluster=H6_sex_dat$Sample.ID, rho = 0.5)
H6_sex_model <- rma.mv(yi, V_H6_sex,
                       mods = ~ male.ratio,
                       random = ~ 1 | No. / Study / effectID, data = H6_sex_dat, method = "REML")
summary(H6_sex_model)