suppressPackageStartupMessages(library(dplyr))

coding_path <- "Coding Sample - Final Meta Data.csv"
coding <- suppressMessages(readr::read_csv(coding_path, skip = 1, show_col_types = FALSE))

# Count each shared cohort once.
cohorts <- coding %>%
  select(No., Study, Sample_ID, Intext_Citation, N, Female, Male, Age_Mean,
         Country, Study_Setting, Group_Type, Game_Type) %>%
  mutate(
    Sample_ID = if_else(No. %in% c(35, 41, 71) | (No. == 38 & Study == "3b"),
                        NA_character_, Sample_ID),
    Age_Mean = as.numeric(na_if(tolower(Age_Mean), "not reported"))
  ) %>%
  distinct()

combine_levels <- function(values) {
  paste(sort(unique(values)), collapse = " / ")
}

study_characteristics <- cohorts %>%
  mutate(Sample = if_else(No. == 10 | (No. == 38 & Study == "1"), Sample_ID, "")) %>%
  group_by(No., Study, Sample, Intext_Citation) %>%
  summarise(
    Male_Proportion = round(sum(Male) / sum(N), 2),
    Age = na_if(round(weighted.mean(Age_Mean, N, na.rm = TRUE), 2), NaN),
    N = sum(N),
    Country = combine_levels(Country),
    Platform = combine_levels(Study_Setting),
    Group = combine_levels(Group_Type),
    Game = combine_levels(Game_Type),
    .groups = "drop"
  ) %>%
  select(Paper = No., Citation = Intext_Citation, Study, Sample,
         N, Male_Proportion, Age, Country, Platform, Group, Game)

study_ages <- cohorts %>%
  group_by(No., Study) %>%
  summarise(Age = weighted.mean(Age_Mean, N, na.rm = TRUE), .groups = "drop")

sample_summary <- cohorts %>%
  summarise(
    Papers = n_distinct(No.),
    Studies = nrow(study_ages),
    N_With_Age_Data = sum(N[!is.na(Age_Mean)]),
    Weighted_Mean_Age = round(weighted.mean(Age_Mean, N, na.rm = TRUE), 2),
    N = sum(N),
    Female = sum(Female),
    Male = sum(Male),
    Female_Percent = round(100 * Female / N, 1),
    Male_Percent = round(100 * Male / N, 1),
    Studies_Reporting_Age = sum(!is.na(study_ages$Age)),
    Min_Study_Mean_Age = round(min(study_ages$Age, na.rm = TRUE), 2),
    Max_Study_Mean_Age = round(max(study_ages$Age, na.rm = TRUE), 2)
  ) %>%
  relocate(N_With_Age_Data, Weighted_Mean_Age, .after = Studies_Reporting_Age)

cat("\nStudy characteristics\n")
study_characteristics %>%
  mutate(Male_Proportion = sprintf("%.2f", Male_Proportion),
         Age = if_else(is.na(Age), "N/A", sprintf("%.2f", Age))) %>%
  print(n = Inf, width = Inf)
cat("\nSample summary\n")
print(data.frame(Statistic = names(sample_summary),
                 Value = as.character(unlist(sample_summary))), row.names = FALSE)
