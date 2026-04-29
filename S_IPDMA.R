# Title:       R Script for TPP IPD MA
# Description: This script analyze the individual participants data meta-data 
#              from the project "Altrusitic and Antisocial Third-Party Punishment
#              in Intergroup Interaction: A Systematic Review and Meta-Analysis"
# Author:      Lingyu MENG
# Date:        2026-04-27
# Version:     1.0
# Notes:       Requires 'Individual Data.xlsx' in the folder.
#              We conducted one-stage IPDMA, by building HLM with cluster at paper,
#              study, participants level.

required_packages <- c("tidyverse", "readxl", "lme4", "clubSandwich", "cowplot", "bruceR", "brms", "ggsignif")

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
IPD_dat_raw <- read_excel("Individual Data.xlsx", sheet = "IPD_Invalid_Excluded", guess_max=Inf)

# Data Preparation ----
IPD_dat <- IPD_dat_raw %>% 
  mutate(P_log = car::logit(Punishment_STD), # proportions remapped to (0.025, 0.975) 
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
    (1 | No. / Study / Participant),
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_50 == "0"),
  weights = Number_of_Trials
)
car::vif(Full_model_unfair) # GVIF of Real_Interaction > 20
summary(Full_model_unfair)
# remove Real_Interaction and Nationality
model_1_unfair<- lmer(
  P_log ~ Study_Setting + Group_Type + Game_Type +
    Age + Sex + Membership_A * Membership_B +
    (1 | No. / Study / Participant),
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_50 == "0"),
  weights = Number_of_Trials
)
car::vif(model_1_unfair)
summary(model_1_unfair)
HLM_summary(model_1_unfair)

# trying remove Group_Type
model_2_unfair<- lmer(
  P_log ~ Study_Setting + Game_Type + Age + Sex + Membership_A * Membership_B +
    (1 | No. / Study / Participant),
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_50 == "0"),
  weights = Number_of_Trials
)
car::vif(model_2_unfair)
summary(model_2_unfair)
anova(model_1_unfair, model_2_unfair) # no significant difference

# trying remove Game_Type
model_3_unfair<- lmer(
  P_log ~ Study_Setting + Age + Sex + Membership_A * Membership_B +
    (1 | No. / Study / Participant),
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_50 == "0"),
  weights = Number_of_Trials
)
car::vif(model_3_unfair)
summary(model_3_unfair)
anova(model_3_unfair, model_2_unfair) # no significant difference

# trying remove Study_Setting
model_4_unfair<- lmer(
  P_log ~ Age + Sex + Membership_A * Membership_B +
    (1 | No. / Study / Participant),
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_50 == "0"),
  weights = Number_of_Trials
)
car::vif(model_4_unfair)
summary(model_4_unfair)
anova(model_3_unfair, model_4_unfair) # no significant difference

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
    (1 | No. / Study / Participant),
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_50 == "1"),
  weights = Number_of_Trials
)
# Model failed to converge with 1 negative eigenvalue: -1.8e-07 
# remove No.
Full_model_fair_1 <- lmer(
  P_log ~ Study_Setting + Group_Type + Game_Type + Real_Interaction +
    Age + Sex + Nationality + Membership_A * Membership_B +
    (1 | Study / Participant),
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_50 == "1"),
  weights = Number_of_Trials
)
car::vif(Full_model_fair_1) # GVIF of Real_Interaction and Nationality > 20
# remove Real_Interaction and Nationality
model_1_fair <- lmer(
  P_log ~ Study_Setting + Group_Type + Game_Type + Age + Sex + Membership_A * Membership_B +
    (1 | Study / Participant),
  weights = Number_of_Trials,
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_50 == "1")
  )
car::vif(model_1_fair)
summary(model_1_fair)

# trying remove Group_Type
model_2_fair <- lmer(
  P_log ~ Study_Setting + Game_Type + Age + Sex + Membership_A * Membership_B +
    (1 | Study / Participant),
  weights = Number_of_Trials,
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_50 == "1")
)
car::vif(model_2_fair)
summary(model_2_fair)
anova(model_1_fair, model_2_fair) # significant difference

# trying remove Study_Setting
model_3_fair <- lmer(
  P_log ~ Group_Type + Game_Type + Age + Sex + Membership_A * Membership_B +
    (1 | Study / Participant),
  weights = Number_of_Trials,
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_50 == "1")
)
car::vif(model_3_fair)
summary(model_3_fair)
anova(model_1_fair, model_3_fair) # no significant difference, p = 0.05707

# trying remove Game_Type
model_4_fair <- lmer(
  P_log ~ Group_Type + Age + Sex + Membership_A * Membership_B +
    (1 | Study / Participant),
  weights = Number_of_Trials,
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_50 == "1")
)
car::vif(model_4_fair)
summary(model_4_fair)
anova(model_4_fair, model_3_fair) # significant difference

model_summary(list(
  model_1_fair,
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
    (1 | No. / Study / Participant),
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_60 == "0"),
  weights = Number_of_Trials
)
car::vif(Full_model_unfair_60) # GVIF of Real_Interaction > 20
summary(Full_model_unfair_60)
# remove Real_Interaction and Nationality
model_1_unfair_60 <- lmer(
  P_log ~ Study_Setting + Group_Type + Game_Type +
    Age + Sex + Membership_A * Membership_B +
    (1 | No. / Study / Participant),
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
# remove No.
Full_model_fair_60 <- lmer(
  P_log ~ Study_Setting + Group_Type + Game_Type + Real_Interaction +
    Age + Sex + Nationality + Membership_A * Membership_B +
    (1 | Study / Participant),
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_60 == "1"),
  weights = Number_of_Trials
)
car::vif(Full_model_fair_60) # GVIF of Real_Interaction and Nationality > 20
# remove Real_Interaction and Nationality
model_1_fair_60 <- lmer(
  P_log ~ Study_Setting + Group_Type + Game_Type + Age + Sex + Membership_A * Membership_B +
    (1 | Study / Participant),
  weights = Number_of_Trials,
  data = IPD_dat_IO %>% 
    filter(Action_A_Binary_60 == "1")
)
car::vif(model_1_fair_60)
HLM_summary(model_1_fair_60)

# EMMEANS ----
emmeans(model_1_unfair,
        specs = ~ Membership_A | Membership_B) %>% 
  pairs()
emmeans(model_1_unfair,
        specs = ~ Membership_B | Membership_A) %>% 
  pairs()
emm_unfair <- emmeans(model_1_unfair,
                      specs = ~ Membership_A | Membership_B) %>% 
  as.data.frame()
slope_unfair <- ggplot(emm_unfair,
       aes(x = Membership_A,
           y = emmean,
           group = Membership_B,
           colour = Membership_B)) +
  geom_point(size = 3) +
  geom_line() +
  geom_signif(comparisons = list(c("Ingroup", "Outgroup")),
              annotation = c("***"), tip_length = 0) +
  geom_signif(comparisons = list(c("Ingroup", "Outgroup")),
              annotation = c("***"), tip_length = 0,
              y_position = -2.1,
              vjust = 2.2, colour = "#66CCEE") +
  annotate(geom = "text", x = 0.95, y = -1.4, label = "***") +
  annotate(geom = "text", x = 2.05, y = -1.55, label = "***") +
  labs(
    title = "Altruistic TPP",
    x = "Focal Group Membership",
    y = "Predicted Punishment (log)"
  ) +
  theme_cowplot() +
  scale_color_manual(values = c(Ingroup = "#EE6677", Outgroup = "#66CCEE")) +
  theme(legend.position = "none")

emmeans(model_1_fair,
        specs = ~ Membership_A | Membership_B) %>% 
  pairs()
emmeans(model_1_fair,
        specs = ~ Membership_B | Membership_A) %>% 
  pairs()
emm_fair <- emmeans(model_1_fair,
                      specs = ~ Membership_A | Membership_B) %>% 
  as.data.frame()
slope_fair <- ggplot(emm_fair,
                       aes(x = Membership_A,
                           y = emmean,
                           group = Membership_B,
                           colour = Membership_B)) +
  geom_point(size = 3) +
  geom_line() +
  geom_signif(comparisons = list(c("Ingroup", "Outgroup")),
              # y_position = 6.655,
              annotation = c("***"), tip_length = 0) +
  geom_signif(comparisons = list(c("Ingroup", "Outgroup")),
              annotation = c("NS."), tip_length = 0,
              y_position = -2.35,
              vjust = 2.2, colour = "#66CCEE") +
  annotate(geom = "text", x = 0.95, y = -2.76, label = "***") +
  annotate(geom = "text", x = 2.05, y = -2.85, label = "***") +
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
