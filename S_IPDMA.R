
# Notes:       Requires 'IPD_Participant_Condition.csv' in the folder.

required_packages <- c("tidyverse", "readxl", "lme4", "clubSandwich", "cowplot", "bruceR", "ggsignif", "rempsyc")

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

if (!dir.exists("tab")) dir.create("tab")
if (!dir.exists("fig")) dir.create("fig")

# read the data
IPD_dat_raw <- readr::read_csv("IPD_Participant_Condition.csv",
                             col_types = readr::cols(Participant = readr::col_character(), Participant_ID = readr::col_character()))

# Data Preparation ----
IPD_dat <- IPD_dat_raw %>% 
  mutate(Study_ID = paste(No., Study, sep = "_"),
         P_log = car::logit(Punishment_Mean_STD), # proportions remapped to (0.025, 0.975) 
         Study_Setting = factor(Study_Setting,
                                levels = c("Lab", "Field", "Online")), # aviod set field as reference
         Group_Type    = as.factor(Group_Type),
         Game_Type     = as.factor(Game_Type),
         Real_Interaction = as.factor(Real_Interaction),
         Membership_A = as.factor(Membership_A),
         Membership_B = as.factor(Membership_B),
         Nationality  = as.factor(Nationality),
         Action_A_Binary_60 = as.factor(Action_A_Binary_60),
         Action_A_Binary_50 = as.factor(Action_A_Binary_50),
         Sex = as.factor(Sex),
         Age = scale(Age))

IPD_dat_IO <- IPD_dat %>% # data exclude neutral group
  filter(Membership_A != "Neutral") %>% 
  filter(Membership_B != "Neutral")

missing_value <- sapply(IPD_dat_IO, function(y) sum(is.na(y)))
missing_value <- missing_value[missing_value > 0]

# Modelling ----
## Unfair ----
IPD_dat_IO %>% 
  filter(Action_A_Binary_50 == "0") %>% 
  select(Study_Setting, Group_Type, Game_Type, Real_Interaction, Nationality) %>% 
  summary()
# Group_Type = Real_Interaction; Prisoner's Dilemma may = Switzerland

Full_model_unfair <- lmer(
  P_log ~ Study_Setting + Group_Type + Game_Type + Real_Interaction +
    Age + Sex + Nationality + Membership_A * Membership_B +
    (1 | Study_ID / Participant_ID),
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_50 == "0"),
  weights = Number_of_Trials
)
# Got: fixed-effect model matrix is rank deficient so dropping 5 columns / coefficients
car::vif(Full_model_unfair) # GVIF^(1/(2*Df)) of Real_Interaction ins Inf
summary(Full_model_unfair)
# If only remove Real_Interaction, will recive fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
# Remove Real_Interaction and game type, then the model is fin
model_1_unfair<- lmer(
  P_log ~ Study_Setting + Group_Type +
    Age + Sex + Membership_A * Membership_B +
    (1 | Study_ID / Participant_ID),
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_50 == "0"),
  weights = Number_of_Trials
)
car::vif(model_1_unfair)
summary(model_1_unfair)
HLM_summary(model_1_unfair)

# Remove Real_Interaction and group type, then the model is fin too
model_2_unfair<- lmer(
  P_log ~ Study_Setting + Game_Type + Age + Sex + Membership_A * Membership_B +
    (1 | Study_ID / Participant_ID),
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_50 == "0"),
  weights = Number_of_Trials
)
car::vif(model_2_unfair)
summary(model_2_unfair)

# trying remove both group type and Game_Type
model_3_unfair<- lmer(
  P_log ~ Study_Setting + Age + Sex + Membership_A * Membership_B +
    (1 | Study_ID / Participant_ID),
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_50 == "0"),
  weights = Number_of_Trials
)
car::vif(model_3_unfair)
summary(model_3_unfair)
anova(model_3_unfair, model_2_unfair) # no significant difference, suggest remove Game_Type

# trying remove Study_Setting
model_4_unfair<- lmer(
  P_log ~ Age + Sex + Membership_A * Membership_B +
    (1 | Study_ID / Participant_ID),
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_50 == "0"),
  weights = Number_of_Trials
)
car::vif(model_4_unfair)
summary(model_4_unfair)
anova(model_3_unfair, model_4_unfair) # significant difference, suggest keep Study_Setting

model_summary(list(
  model_1_unfair,
  model_2_unfair,
  model_3_unfair,
  model_4_unfair),
  file="tab/IPD_unfair.docx")

## Fair ----
IPD_dat_IO %>% 
  filter(Action_A_Binary_50 == "1") %>% 
  select(Study_Setting, Group_Type, Game_Type, Real_Interaction, Nationality) %>% 
  summary()
# Experimental = China; Prisoner's Dilemma may = Switzerland; Papua New Guinea = Field
Full_model_fair <- lmer(
  P_log ~ Study_Setting + Group_Type + Game_Type + Real_Interaction +
    Age + Sex + Nationality + Membership_A * Membership_B +
    (1 | Study_ID / Participant_ID),
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_50 == "1"),
  weights = Number_of_Trials
)
# Got: fixed-effect model matrix is rank deficient so dropping 5 columns / coefficients``
car::vif(Full_model_fair) # multiple GVIF^(1/(2*Df)) are Inf, including Real_Interaction and Nationality
# Model with study setting and group type
# adding other will lead to fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
model_1_fair <- lmer(
  P_log ~ Study_Setting + Group_Type + Age + Sex + Membership_A * Membership_B +
    (1 | Study_ID / Participant_ID),
  weights = Number_of_Trials,
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_50 == "1")
  )
car::vif(model_1_fair)
summary(model_1_fair)

# Model with study setting and game type
model_2_fair <- lmer(
  P_log ~ Study_Setting + Game_Type + Age + Sex + Membership_A * Membership_B +
    (1 | Study_ID / Participant_ID),
  weights = Number_of_Trials,
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_50 == "1")
)
car::vif(model_2_fair)
summary(model_2_fair)

# trying remove Study_Setting
model_3_fair <- lmer(
  P_log ~ Group_Type + Game_Type + Age + Sex + Membership_A * Membership_B +
    (1 | Study_ID / Participant_ID),
  weights = Number_of_Trials,
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_50 == "1")
)
car::vif(model_3_fair)
summary(model_3_fair)
anova(model_1_fair, model_3_fair) # no significant difference, suggest remove Study_Setting

# trying remove Game_Type
model_4_fair <- lmer(
  P_log ~ Group_Type + Age + Sex + Membership_A * Membership_B +
    (1 | Study_ID / Participant_ID),
  weights = Number_of_Trials,
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_50 == "1")
)
car::vif(model_4_fair)
summary(model_4_fair)
anova(model_4_fair, model_3_fair) # significant difference

model_summary(list(
  model_1_fair,
  model_2_fair,
  model_3_fair,
  model_4_fair),
  file="tab/IPD_fair.docx")

HLM_summary(model_1_unfair)
HLM_summary(model_1_fair)

lmerTest::as_lmerModLmerTest(model_1_unfair) %>%
  anova() %>%
  rownames_to_column("Effects") %>%
  nice_table(col.format.p = 7) %>% 
  flextable::save_as_docx(., path = "tab/IPDMA_unfair_anova.docx")
lmerTest::as_lmerModLmerTest(model_1_fair) %>% 
  anova() %>%
  rownames_to_column("Effects") %>%
  nice_table(col.format.p = 7) %>% 
  flextable::save_as_docx(., path = "tab/IPDMA_fair_anova.docx")


# 60 as fair ----
## Unfair ----
IPD_dat_IO %>% 
  filter(Action_A_Binary_60 == "0") %>% 
  select(Study_Setting, Group_Type, Game_Type, Real_Interaction, Nationality) %>% 
  summary()
# Group_Type = Real_Interaction?; Prisoner's Dilemma may = Switzerland

Full_model_unfair_60 <- lmer(
  P_log ~ Study_Setting + Group_Type + Game_Type + Real_Interaction +
    Age + Sex + Nationality + Membership_A * Membership_B +
    (1 | Study_ID / Participant_ID),
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_60 == "0"),
  weights = Number_of_Trials
)
# Got: fixed-effect model matrix is rank deficient so dropping 5 columns / coefficients
car::vif(Full_model_unfair_60) # GVIF^(1/(2*Df)) of Real_Interaction and Game_Type are Inf
summary(Full_model_unfair_60)
# Remove Real_Interaction and Game_Type, then the model is fin
model_1_unfair_60 <- lmer(
  P_log ~ Study_Setting + Group_Type +
    Age + Sex + Membership_A * Membership_B +
    (1 | Study_ID / Participant_ID),
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_60 == "0"),
  weights = Number_of_Trials
)
car::vif(model_1_unfair_60)
summary(model_1_unfair_60)
HLM_summary(model_1_unfair_60)

## Fair ----
IPD_dat_IO %>% 
  filter(Action_A_Binary_60 == "1") %>% 
  select(Study_Setting, Group_Type, Game_Type, Real_Interaction, Nationality) %>% 
  summary()
# Experimental = China; Prisoner's Dilemma may = Switzerland; Papua New Guinea = Field
Full_model_fair_60 <- lmer(
  P_log ~ Study_Setting + Group_Type + Game_Type + Real_Interaction +
    Age + Sex + Nationality + Membership_A * Membership_B +
    (1 | Study_ID / Participant_ID),
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_60 == "1"),
  weights = Number_of_Trials
)
# Got: fixed-effect model matrix is rank deficient so dropping 5 columns / coefficients
car::vif(Full_model_fair_60) # multiple GVIF^(1/(2*Df)) are Inf, including Real_Interaction and Nationality
# Model with study setting and group type
# adding other will lead to fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
model_1_fair_60 <- lmer(
  P_log ~ Study_Setting + Group_Type + Age + Sex + Membership_A * Membership_B +
    (1 | Study_ID / Participant_ID),
  weights = Number_of_Trials,
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_60 == "1")
)
car::vif(model_1_fair_60)
HLM_summary(model_1_fair_60)

# EMMEANS ----
significance_labels <- function(comparisons, group) {
  results <- as.data.frame(comparisons)
  labels <- case_when(
    results$p.value < 0.001 ~ "***",
    results$p.value < 0.01 ~ "**",
    results$p.value < 0.05 ~ "*",
    TRUE ~ "NS"
  )
  setNames(labels, results[[group]])
}

EMM_unfair <- emmeans(model_1_unfair,
        specs = ~ Membership_A | Membership_B)
EMM_unfair %>% as.data.frame() %>%
  mutate(
    across(c(emmean, lower.CL, upper.CL),
           ~ boot::inv.logit(.x))
  ) %>%
  mutate(
    across(c(emmean, lower.CL, upper.CL),
           ~ (.x - 0.025) / (1 - 2 * 0.025))
  ) %>% nice_table() %>% flextable::save_as_docx(., path = "tab/IPDMA_unfair_EMMEANS.docx")

unfair_pairs_by_B <- pairs(EMM_unfair)
unfair_pairs_by_A <- pairs(emmeans(model_1_unfair, ~ Membership_B | Membership_A))
unfair_pairs_by_B %>%
  nice_table() %>% flextable::save_as_docx(., path = "tab/IPDMA_unfair_EMMEANS_B.docx")
unfair_pairs_by_A %>%
  nice_table() %>% flextable::save_as_docx(., path = "tab/IPDMA_unfair_EMMEANS_A.docx")

unfair_labels_by_B <- significance_labels(unfair_pairs_by_B, "Membership_B")
unfair_labels_by_A <- significance_labels(unfair_pairs_by_A, "Membership_A")
  
slope_unfair <- ggplot(as.data.frame(EMM_unfair),
       aes(x = Membership_A,
           y = emmean,
           group = Membership_B,
           colour = Membership_B)) +
  geom_point(size = 3) +
  geom_line() +
  geom_signif(comparisons = list(c("Ingroup", "Outgroup")),
              annotation = unfair_labels_by_B[["Ingroup"]], tip_length = 0) +
  geom_signif(comparisons = list(c("Ingroup", "Outgroup")),
              annotation = unfair_labels_by_B[["Outgroup"]], tip_length = 0,
              y_position = -1.5,
              vjust = 2.2, colour = "#66CCEE") +
  annotate(geom = "text", x = 0.95, y = -1.15, label = unfair_labels_by_A[["Ingroup"]]) +
  annotate(geom = "text", x = 2.05, y = -1.0, label = unfair_labels_by_A[["Outgroup"]]) +
  labs(
    title = "Altruistic TPP",
    x = "Focal Group Membership",
    y = "Predicted Punishment (log)"
  ) +
  theme_cowplot() +
  scale_color_manual(values = c(Ingroup = "#EE6677", Outgroup = "#66CCEE")) +
  theme(legend.position = "none")

EMM_fair <- emmeans(model_1_fair,
        specs = ~ Membership_A | Membership_B)
EMM_fair %>% as.data.frame() %>%
  mutate(
    across(c(emmean, lower.CL, upper.CL),
           ~ boot::inv.logit(.x))
  ) %>%
  mutate(
    across(c(emmean, lower.CL, upper.CL),
           ~ (.x - 0.025) / (1 - 2 * 0.025))
  ) %>% nice_table() %>% flextable::save_as_docx(., path = "tab/IPDMA_fair_EMMEANS.docx")

fair_pairs_by_B <- pairs(EMM_fair)
fair_pairs_by_A <- pairs(emmeans(model_1_fair, ~ Membership_B | Membership_A))
fair_pairs_by_B %>%
  nice_table() %>% flextable::save_as_docx(., path = "tab/IPDMA_fair_EMMEANS_B.docx")
fair_pairs_by_A %>%
  nice_table() %>% flextable::save_as_docx(., path = "tab/IPDMA_fair_EMMEANS_A.docx")

fair_labels_by_B <- significance_labels(fair_pairs_by_B, "Membership_B")
fair_labels_by_A <- significance_labels(fair_pairs_by_A, "Membership_A")
  
slope_fair <- ggplot(EMM_fair %>% as.data.frame(),
                       aes(x = Membership_A,
                           y = emmean,
                           group = Membership_B,
                           colour = Membership_B)) +
  geom_point(size = 3) +
  geom_line() +
  geom_signif(comparisons = list(c("Ingroup", "Outgroup")),
              # y_position = 6.655,
              annotation = fair_labels_by_B[["Ingroup"]], tip_length = 0) +
  geom_signif(comparisons = list(c("Ingroup", "Outgroup")),
              annotation = fair_labels_by_B[["Outgroup"]], tip_length = 0,
              y_position = -2.36,
              vjust = 2.2, colour = "#66CCEE") +
  annotate(geom = "text", x = 0.95, y = -2.74, label = fair_labels_by_A[["Ingroup"]]) +
  annotate(geom = "text", x = 2.05, y = -2.82, label = fair_labels_by_A[["Outgroup"]]) +
  labs(
    title = "Antisocial TPP",
    x = "Focal Group Membership",
    y = "Predicted Punishment (log)",
    colour = 
    "Partner
Group
Membership"
  ) +
  theme_cowplot() +
  scale_color_manual(values = c(Ingroup = "#EE6677", Outgroup = "#66CCEE"))

sample_slope <- plot_grid(
  slope_unfair, slope_fair,
  labels = c('a', 'b'),
  rel_widths = c(1,1.25))

ggsave("fig/S_sample_slopes.jpeg", sample_slope, width = 12.5, height = 5)

model_1_unfair <- lmerTest::as_lmerModLmerTest(model_1_unfair)
model_1_fair <- lmerTest::as_lmerModLmerTest(model_1_fair)
model_1_unfair_60 <- lmerTest::as_lmerModLmerTest(model_1_unfair_60)
model_1_fair_60 <- lmerTest::as_lmerModLmerTest(model_1_fair_60)

IPD_screened_sample_keys <- distinct(IPD_dat_IO, No., Study, Participant, Participant_ID)
IPD_model_sample_keys <- lapply(
  list(unfair_50 = model_1_unfair, fair_50 = model_1_fair,
       unfair_60 = model_1_unfair_60, fair_60 = model_1_fair_60),
  function(model) {
    ids <- stats::model.frame(model)$Participant_ID
    filter(IPD_screened_sample_keys, Participant_ID %in% ids)
  }
)
IPD_sample_keys <- distinct(bind_rows(IPD_model_sample_keys[c("unfair_50", "fair_50")]))

save(IPD_sample_keys, IPD_screened_sample_keys, IPD_model_sample_keys,
     Full_model_unfair, Full_model_fair, Full_model_unfair_60, Full_model_fair_60,
     model_1_unfair, model_2_unfair, model_3_unfair, model_4_unfair,
     model_1_fair, model_2_fair, model_3_fair, model_4_fair,
     model_1_unfair_60, model_1_fair_60,
     EMM_unfair, EMM_fair, unfair_pairs_by_A, unfair_pairs_by_B,
     fair_pairs_by_A, fair_pairs_by_B, sample_slope,
     file = "IPD_meta.rds", compress = "xz")
