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

# H1 The Effects of Violator’s Identity on Punishment
fair_condition <- c("Fair", "Cooperate", "Neutral", "Generous")
unfair_condition <- c("Unfair", "Defect", "Selfish", "Cheat")
H1_data <- MA_dat %>%
  filter(
    `Fairness/Cooperativeness of Person A` %in% unfair_condition,
    `Person A affiliation` %in% c("Ingroup", "Outgroup")
  ) %>%
  select(
    No., Study, `Person A affiliation`, `Recipient/Player B Affliation`, `Sample ID`,
    `Fairness/Cooperativeness of Person A`, `Fairness/Cooperativeness of Person B`,
    `P_Average Invested Amount`, `P_SD invested amount`, `N by condition`, `Total Trials`, `Intext-Citation`
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

# H2 (Deterrence perspective):
H2_data <- MA_dat %>%
  filter(
    `Fairness/Cooperativeness of Person A` %in% unfair_condition,
    `Recipient/Player B Affliation` %in% c("Ingroup", "Outgroup")
  ) %>%
  # group_by(No., Study, `Recipient/Player B Affliation`) %>%
  select(
    No., Study, `Person A affiliation`, `Recipient/Player B Affliation`, `Sample ID`,
    `Fairness/Cooperativeness of Person A`, `Fairness/Cooperativeness of Person B`,
    `P_Average Invested Amount`, `P_SD invested amount`, `N by condition`, `Total Trials`, `Intext-Citation`
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

# # Interaction Effect between Violator’s Identity and Victim Salience
# # The Mantel-Haenszel method
# MH_dat <- MA_dat %>%
#   group_by(No., Study) %>%
#   # keep only studies that have all four A×B combinations
#   filter(
#     `Fairness/Cooperativeness of Person A` %in% unfair_condition,
#     n_distinct(`Person A affiliation`) >= 2,
#     n_distinct(`Recipient/Player B Affliation`) >= 2,
#     all(expand.grid(
#       `Person A affiliation` = c("Ingroup", "Outgroup"),
#       `Recipient/Player B Affliation` = c("Ingroup", "Outgroup")
#     ) %>%
#       apply(1, paste, collapse = "_") %in%
#       paste(`Person A affiliation`, `Recipient/Player B Affliation`, sep = "_"))
#   ) %>%
#   ungroup()
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

# The regression method
# H3_model <- MH_dat %>%
#   mutate(
#     `Total Trials` = as.numeric(`Total Trials`),
#     `Person A` = factor(`Person A affiliation`, levels = c("Ingroup", "Outgroup")),
#     `Person B` = factor(`Recipient/Player B Affliation`, levels = c("Ingroup", "Outgroup")),
#     P_SV = as.numeric(`P_SD invested amount`)^2 * `Total Trials` / ( `Total Trials` - 1), # Sampling variance of punishment amount
#     `P_Average Invested Amount` = as.numeric(`P_Average Invested Amount`)
#   ) %>%
#   rma.mv(yi = `P_Average Invested Amount`, V = P_SV,
#   mods = ~ `Person A` * `Person B`,
#   random = ~ 1 | No. / Study,
#   data = .,
#   method = "REML",
# )
# summary(H3_model)

################ Fair Condition Analysis ################
# H5 The Effects of Violator’s Identity on Punishment
H5_data <- MA_dat %>%
  filter(
    `Fairness/Cooperativeness of Person A` %in% fair_condition,
    `Person A affiliation` %in% c("Ingroup", "Outgroup")
  ) %>%
  select(
    No., Study, `Person A affiliation`, `Recipient/Player B Affliation`, `Sample ID`,
    `Fairness/Cooperativeness of Person A`, `Fairness/Cooperativeness of Person B`,
    `P_Average Invested Amount`, `P_SD invested amount`, `N by condition`, `Total Trials`, `Intext-Citation`
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

# H6 
H6_data <- MA_dat %>%
  filter(
    `Fairness/Cooperativeness of Person A` %in% fair_condition,
    `Recipient/Player B Affliation` %in% c("Ingroup", "Outgroup")
  ) %>%
  # group_by(No., Study, `Recipient/Player B Affliation`) %>%
  select(
    No., Study, `Person A affiliation`, `Recipient/Player B Affliation`, `Sample ID`,
    `Fairness/Cooperativeness of Person A`, `Fairness/Cooperativeness of Person B`,
    `P_Average Invested Amount`, `P_SD invested amount`, `N by condition`, `Total Trials`, `Intext-Citation`
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
V_H6 <- vcalc(H6_data$vi, cluster=H6_data$Sample.ID, rho = 0.5) # we don't have sample ID right now. assuming same sample within study
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
