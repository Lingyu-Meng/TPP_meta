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

# H1 The Effects of Aggressor’s Identity on Punishment
fair_condition <- c("Fair", "Cooperate", "Neutral", "Generous")
unfair_condition <- c("Unfair", "Defect", "Selfish", "Cheat")
H1_data <- MA_dat %>%
  filter(
    `Fairness/Cooperativeness of Person A` %in% unfair_condition,
    `Person A affiliation` %in% c("Ingroup", "Outgroup")
  ) %>%
  group_by(No., Study, `Person A affiliation`) %>% # For main effect
  summarise(
    # weighted mean across B conditions
    mean = weighted.mean(`P_Average Invested Amount`, `N by condition`),
    # pooled SD across B conditions
    sd = sqrt(sum((`N by condition` - 1) * `P_SD invested amount`^2) / sum(`N by condition` - 1)),
    n = sum(`N by condition`),
    .groups = "drop"
  ) %>%
  filter(!is.na(mean)) %>% # Exclude mean is NA
  pivot_wider(
    names_from = `Person A affiliation`,
    values_from = c(mean, sd, n)
  ) %>%
  escalc(
    measure = "SMDH", # Hedges' g
    m1i = mean_Outgroup, sd1i = sd_Outgroup, n1i = n_Outgroup,
    m2i = mean_Ingroup, sd2i = sd_Ingroup, n2i = n_Ingroup,
    data = .
  )
H1_model <- rma.mv(yi, vi, random = ~ 1 | No. / Study, data = H1_data, method = "REML")
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

# H2 (Victim Salience Effect):
H2_data <- MA_dat %>%
  filter(
    `Fairness/Cooperativeness of Person A` %in% unfair_condition,
    `Person A affiliation` %in% c("Ingroup", "Outgroup")
  ) %>%
  group_by(No., Study, `Recipient/Player B Affliation`) %>%
  summarise(
    mean = weighted.mean(`P_Average Invested Amount`, `N by condition`),
    sd = sqrt(sum((`N by condition` - 1) * `P_SD invested amount`^2) / sum(`N by condition` - 1)),
    n = sum(`N by condition`),
    .groups = "drop"
  ) %>%
  filter(!is.na(mean)) %>%
  pivot_wider(
    names_from = `Recipient/Player B Affliation`,
    values_from = c(mean, sd, n)
  ) %>%
  escalc(
    measure = "SMDH", # Hedges' g
    m1i = mean_Ingroup, sd1i = sd_Ingroup, n1i = n_Ingroup,
    m2i = mean_Outgroup, sd2i = sd_Outgroup, n2i = n_Outgroup,
    data = .
  )
H2_model <- rma.mv(yi, vi, random = ~ 1 | No. / Study, data = H2_data, method = "REML")
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