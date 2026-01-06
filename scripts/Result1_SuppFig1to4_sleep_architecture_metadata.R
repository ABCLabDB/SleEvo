############################################################
# Result 1 – Supplementary Figures 1–4
# Metadata-based characterization of sleep architecture
# across species
#
# This script generates Supplementary Figures 1–4 of Result 1,
# summarizing:
#   (1) total sleep time
#   (2) NREM sleep ratio
#   (3) sleep frequency
#   (4) sleep timing patterns
#
# Curated species-level sleep metadata is required.
# The data file is not included in this repository.
# See data/README.md for details.
############################################################

########################
# Load required packages
########################
required_packages <- c(
  "data.table",
  "dplyr",
  "ggplot2",
  "ggpubr",
  "DescTools"
)

invisible(lapply(required_packages, library, character.only = TRUE))

########################
# Load metadata
########################
DATA_DIR <- "data"

species_metadata <- fread(
  file.path(DATA_DIR, "species_sleep_metadata.txt")
) |> as.data.frame()

############################################################
# Supplementary Figure 1
# Total sleep time across species
############################################################
tst_meta <- species_metadata

longest_species <- tst_meta %>%
  filter(Total_sleep_time_per_day == max(Total_sleep_time_per_day, na.rm = TRUE)) %>%
  pull(Species_name_ensembl)

shortest_species <- tst_meta %>%
  filter(Total_sleep_time_per_day == min(Total_sleep_time_per_day, na.rm = TRUE)) %>%
  pull(Species_name_ensembl)

tst_meta <- tst_meta %>%
  mutate(
    Sleep_Time = round(Total_sleep_time_per_day, 1),
    Group = case_when(
      Species_name_ensembl %in% longest_species ~ "Longest Sleep",
      Species_name_ensembl %in% shortest_species ~ "Shortest Sleep",
      Order == "Primates" & Species_name_ensembl != "Homo_sapiens" ~ "Primates",
      Species_name_ensembl == "Homo_sapiens" ~ "Human",
      TRUE ~ "Others"
    )
  )

avg_sleep_all <- mean(tst_meta$Sleep_Time, na.rm = TRUE)
avg_sleep_primates <- tst_meta %>%
  filter(Order == "Primates") %>%
  summarize(m = mean(Sleep_Time, na.rm = TRUE)) %>%
  pull(m)

supfig1 <- ggbarplot(
  tst_meta,
  x = "Species_name_ensembl",
  y = "Sleep_Time",
  fill = "Group",
  sort.val = "desc",
  x.text.angle = 60
) +
  geom_hline(yintercept = avg_sleep_all, linetype = "dashed", color = "red") +
  geom_hline(yintercept = avg_sleep_primates, linetype = "dashed", color = "#f1b321") +
  labs(y = "Total sleep time (hours)", fill = "Group")

supfig1


############################################################
# Supplementary Figure 2
# NREM sleep ratio
############################################################
nrem_meta <- species_metadata %>%
  filter(!is.na(Percentage_of_NREM_time_per_day)) %>%
  mutate(
    NREM_ratio = round(Percentage_of_NREM_time_per_day, 2),
    Group = case_when(
      Species_name_ensembl == "Homo_sapiens" ~ "Human",
      Order == "Primates" ~ "Primates",
      Class == "Aves" ~ "Aves",
      TRUE ~ "Others"
    )
  )

avg_nrem_all <- mean(nrem_meta$NREM_ratio, na.rm = TRUE)
avg_nrem_primates <- nrem_meta %>%
  filter(Order == "Primates") %>%
  summarize(m = mean(NREM_ratio, na.rm = TRUE)) %>%
  pull(m)

avg_nrem_aves <- nrem_meta %>%
  filter(Class == "Aves") %>%
  summarize(m = mean(NREM_ratio, na.rm = TRUE)) %>%
  pull(m)

supfig2 <- ggbarplot(
  nrem_meta,
  x = "Species_name_ensembl",
  y = "NREM_ratio",
  fill = "Group",
  sort.val = "desc",
  x.text.angle = 60
) +
  geom_hline(yintercept = avg_nrem_all, linetype = "dashed", color = "red") +
  geom_hline(yintercept = avg_nrem_primates, linetype = "dashed", color = "#f1b321") +
  geom_hline(yintercept = avg_nrem_aves, linetype = "dashed", color = "#194a7a") +
  labs(y = "NREM sleep ratio", fill = "Group")

supfig2


############################################################
# Supplementary Figure 3
# Sleep frequency
############################################################
sleep_freq_meta <- species_metadata %>%
  filter(!is.na(Number_of_sleep_times_per_day)) %>%
  mutate(
    Number_of_sleep_times_per_day = factor(
      Number_of_sleep_times_per_day,
      levels = c("Once", "More than twice")
    )
  )

freq_summary <- sleep_freq_meta %>%
  count(Class, Order, Number_of_sleep_times_per_day)

supfig3 <- ggplot(
  freq_summary,
  aes(
    x = Number_of_sleep_times_per_day,
    y = n,
    fill = Number_of_sleep_times_per_day
  )
) +
  geom_col() +
  facet_wrap(~ Class, scales = "free_y") +
  labs(
    x = "Sleep frequency",
    y = "Number of species",
    fill = "Sleep frequency"
  )

supfig3


############################################################
# Supplementary Figure 4
# Sleep timing patterns
############################################################
sleep_timing_meta <- species_metadata %>%
  filter(!is.na(Sleep_timing_per_day)) %>%
  mutate(
    Sleep_timing_per_day = factor(
      Sleep_timing_per_day,
      levels = c("Sleep at night", "Sleep at daytime", "Sleep at anytime")
    )
  )

timing_summary_all <- sleep_timing_meta %>%
  count(Sleep_timing_per_day)

supfig4_all <- ggpie(
  timing_summary_all,
  x = "n",
  label = "Sleep_timing_per_day",
  title = "Sleep timing across species"
)

supfig4_all


############################################################
# End of script
############################################################
