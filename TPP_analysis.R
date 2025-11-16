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
MA_dat <- read_csv("Coding Sample - Coding Merged.csv", skip = 1)

# H1 The Effects of Violator’s Identity on Punishment
fair_condition <- c("Fair", "Cooperate", "Neutral", "Generous")
unfair_condition <- c("Unfair", "Defect", "Selfish", "Cheat")
H1_data <- MA_dat %>%
  filter(
    `Fairness/Cooperativeness of Person A` %in% unfair_condition,
    `Person A affiliation` %in% c("Ingroup", "Outgroup")
  ) %>%
  select(
    No., Study, `Person A affiliation`, `Recipient/Player B Affliation`,
    Levels, `Fairness/Cooperativeness of Person A`, `Fairness/Cooperativeness of Person B`,
    `P_Average Invested Amount`, `P_SD invested amount`, `N by condition`, `Total Trials`
  ) %>%
  filter(!is.na(`P_Average Invested Amount`)) %>% # Exclude mean is NA
  pivot_wider(
    names_from = `Person A affiliation`,
    values_from = c(`P_Average Invested Amount`, `P_SD invested amount`, `N by condition`,`Total Trials`)
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
V_H1 <- vcalc(H1_data$vi, cluster=H1_data$No., subgroup = H1_data$Study, rho = 0.5) # we don't have sample ID right now. assuming same sample within study
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
forest(H1_model, slab = paste(H1_data$No., H1_data$Study, sep = "."))
# Save fig in folder fig
if (!dir.exists("fig")) {
  dir.create("fig")
}
# saving forest plot to a file
tiff("fig/H1_Forest.tiff", height = 12, width = 8, units = "in", res = 300)
forest(H1_model, slab = paste(H1_data$No., H1_data$Study, sep = "."))
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
    No., Study, `Person A affiliation`, `Recipient/Player B Affliation`,
    Levels, `Fairness/Cooperativeness of Person A`, `Fairness/Cooperativeness of Person B`,
    `P_Average Invested Amount`, `P_SD invested amount`, `N by condition`, `Total Trials`
  ) %>%
  filter(!is.na(`P_Average Invested Amount`)) %>% # Exclude mean is NA
  pivot_wider(
    names_from = `Recipient/Player B Affliation`,
    values_from = c(`P_Average Invested Amount`, `P_SD invested amount`, `N by condition`,`Total Trials`)
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
V_H2 <- vcalc(H2_data$vi, cluster=H2_data$No., subgroup = H2_data$Study, rho = 0.5) # we don't have sample ID right now. assuming same sample within study
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
H2_f <- forest(H2_model, slab = paste(H2_data$No., H2_data$Study, sep = "."))
# saving forest plot to a file
tiff("fig/H2_Forest.tiff", height = 12, width = 8, units = "in", res = 300)
forest(H2_model, slab = paste(H2_data$No., H2_data$Study, sep = "."))
dev.off()

# funnel plot
metafor::funnel(H2_model)

# Interaction Effect between Violator’s Identity and Victim Salience
# The Mantel-Haenszel method
MH_dat <- MA_dat %>%
  group_by(No., Study) %>%
  # keep only studies that have all four A×B combinations
  filter(
    `Fairness/Cooperativeness of Person A` %in% unfair_condition,
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

I_I <- MH_dat %>%
  filter(
    `Person A affiliation` == "Ingroup",
    `Recipient/Player B Affliation` == "Ingroup"
  )
I_O <- MH_dat %>%
  filter(
    `Person A affiliation` == "Ingroup",
    `Recipient/Player B Affliation` == "Outgroup"
  )
O_I <- MH_dat %>%
  filter(
    `Person A affiliation` == "Outgroup",
    `Recipient/Player B Affliation` == "Ingroup"
  )
O_O <- MH_dat %>%
  filter(
    `Person A affiliation` == "Outgroup",
    `Recipient/Player B Affliation` == "Outgroup"
  )
n1i <- I_I$`N by condition` + I_O$`N by condition`
n2i <- O_I$`N by condition` + O_O$`N by condition`

MH_model <- rma.mh(I_I$P_Frequency, #  I O
  I_O$P_Frequency, # I  a b
  O_I$P_Frequency, # O  c d
  O_O$P_Frequency,
  n1i,
  n2i,
  measure = "RR"
)
summary(MH_model)
# k = 8 studies included
MH_model$beta
forest(MH_model)
# RR = (a/(a+b)) / (c/(c+d))
# when person A is an ingroup member punishment probablity

# The regression method
H3_model <- MH_dat %>%
  mutate(
    `Total Trials` = as.numeric(`Total Trials`),
    `Person A` = factor(`Person A affiliation`, levels = c("Ingroup", "Outgroup")),
    `Person B` = factor(`Recipient/Player B Affliation`, levels = c("Ingroup", "Outgroup")),
    P_SV = `P_SD invested amount`^2 * `Total Trials` / ( `Total Trials` - 1) # Sampling variance of punishment amount
  ) %>%
  rma.mv(yi = `P_Average Invested Amount`, V = P_SV,
  mods = ~ `Person A` * `Person B`,
  random = ~ 1 | No. / Study,
  data = .,
  method = "REML",
)
summary(H3_model)