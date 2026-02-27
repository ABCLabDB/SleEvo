############################################################
# Result 1 – Supplementary Figures 1–6
# Metadata-based characterization of sleep architecture
# across species
#
# This script generates Supplementary Figures 1–6 of Result 1,
# summarizing:
#   (1) Total sleep time
#   (2) NREM sleep ratio
#   (3) Number of sleep times per day (=sleep frequency)
#   (4) sleep timing patterns
#   (5) association between sleep timing and sleep frequency
#       in primates (Mouse lemur highlighted)
#   (6) global association between sleep timing and sleep frequency
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
  "DescTools",
  "ggmosaic",
  "scales",
  "forcats"
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

tst_meta <- species_metadata %>%
  mutate(
    Sleep_Time = round(Total_sleep_time_per_day, 1),
    Taxonomy_Class = Class  # modify if column name differs
  ) %>%
  arrange(desc(Sleep_Time)) %>%
  mutate(
    Species_name_ensembl = fct_inorder(Species_name_ensembl)
  )


avg_sleep_all <- mean(tst_meta$Sleep_Time, na.rm = TRUE)
class_colors <- c(
  "Mammalia" = "#1f4e79",
  "Aves" = "#e6b800",
  "Actinopteri" = "#bfbfbf",
  "Insecta" = "#e07a3f",
  "Lepidosauria" = "#5b7f75",
  "Sarcopterygii" = "#e8d8a9"
)

supfig1 <- ggplot(tst_meta,
                  aes(x = Species_name_ensembl,
                      y = Sleep_Time)) +
  
  # Background scaffold (fixed 24h reference)
  geom_col(aes(y = 24),
           fill = "#e6e6e6",
           width = 0.6) +
  
  # Colored bar representing actual sleep time
  geom_col(aes(y = Sleep_Time,
               fill = Taxonomy_Class),
           width = 0.6) +
  
  # Overlay colored points (taxonomy class)
  geom_point(aes(y = Sleep_Time,
                 color = Taxonomy_Class),
             size = 5,
             show.legend = FALSE) +
  
  # Global average sleep time
  geom_hline(yintercept = avg_sleep_all,
             linetype = "dashed",
             color = "red",
             linewidth = 1) +
  
  # Annotation for average line
  annotate("text",
           x = length(unique(tst_meta$Species_name_ensembl)) * 0.85,
           y = avg_sleep_all + 0.6,
           label = paste0("Average Total Sleep Time = ",
                          round(avg_sleep_all, 2), " h"),
           size = 4) +
  
  # Apply taxonomy color scale
  scale_fill_manual(values = class_colors) +
  scale_color_manual(values = class_colors) +
  
  # Axis labels and legend title
  labs(
    y = "Total Sleep Time (hours per day)",
    x = NULL,
    color = "Taxonomy Class"
  ) +
  
  # Fix y-axis to biological maximum range
  coord_cartesian(ylim = c(0, 24)) +
  
  # Clean theme for publication
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "top"
  )

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
# Supplementary Figure 5
# Association between sleep timing and sleep frequency
# in primates (Mouse lemur highlighted)
############################################################
timing_times_Meta <- species_Metadata[!is.na(species_Metadata$Sleep_timing_per_day) & !is.na(species_Metadata$Number_of_sleep_times_per_day), ]
timing_times_Meta$Sleep_Timing <- factor(timing_times_Meta$Sleep_timing_per_day)
timing_times_Meta$Number_of_sleep_times_per_day <- factor(timing_times_Meta$Number_of_sleep_times_per_day)
rownames(timing_times_Meta) <- NULL

fisher.test(timing_times_Meta$Number_of_sleep_times_per_day, timing_times_Meta$Sleep_Timing)

Primates <- timing_times_Meta[which(timing_times_Meta$Order =="Primates"),]

set.seed(42)  

df_plot <- Primates %>%
  select(Species_name_ensembl, Species_symbol_name_ensembl,
         Number_of_sleep_times_per_day, Sleep_Timing) %>%
  mutate(
    Group = ifelse(Species_symbol_name_ensembl == "Microcebus murinus", "Mouse lemur", "Other primates"),
    
    x_base = recode(Number_of_sleep_times_per_day,
                    "Once" = 1, "More than twice" = 2),
    y_base = recode(Sleep_Timing,
                    "Sleep at night" = 1, "Sleep at daytime" = 2),
    
    x_jitter = x_base + runif(n(), -0.2, 0.2),
    y_jitter = y_base + runif(n(), -0.2, 0.2)
  )

df_plot <- df_plot %>%
  mutate(
    jitter_x = as.numeric(factor(Number_of_sleep_times_per_day)) + runif(n(), -0.2, 0.2),
    jitter_y = as.numeric(factor(Sleep_Timing)) + runif(n(), -0.2, 0.2)
  )

df_plot <- df_plot %>%
  mutate(
    Number_of_sleep_times_per_day = factor(Number_of_sleep_times_per_day,
                                           levels = c("Once", "More than twice")),
    Sleep_Timing = factor(Sleep_Timing,
                          levels = c("Sleep at night", "Sleep at daytime"))
  )

primates_sleep_pattern <- ggplot(df_plot, aes(x = x_jitter, y = y_jitter)) +
  geom_point(aes(fill = Group), shape = 21, size = 5, color = "black", stroke = 1.5) +
  geom_text(aes(label = Species_name_ensembl), vjust = -1.2, size = 3) +
  
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "gray30") +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "gray30") +
  
  scale_x_continuous(breaks = c(1, 2), labels = c("Once", "More than twice"), expand = c(0.1, 0.1)) +
  scale_y_continuous(breaks = c(1, 2), labels = c("Sleep at night", "Sleep at daytime"), expand = c(0.1, 0.1)) +
  
  scale_fill_manual(values = c("Mouse lemur" = "#E64B35", "Other primates" = "gray70")) +
  
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    panel.grid = element_blank()
    # legend.position = "top"
  ) +
  
  labs(
    title = "",
    x = "Sleep Frequency",
    y = "Sleep Timing",
    fill = "Species"
  )

primates_sleep_pattern


############################################################
# Supplementary Figure 6
# Global association between sleep timing and sleep frequency
############################################################
custom_colors <- c(
  "Sleep at night" = "#1f77b4",
  "Sleep at daytime" = "#ff7f0e",
  "Sleep at anytime" = "#2ca02c"
)

p <- ggplot(data = timing_times_Meta) +
  geom_mosaic(
    aes(x = product(Number_of_sleep_times_per_day),
        fill = Sleep_Timing,
        weight = 1),
    na.rm = TRUE,
    color = "black", linewidth = 0.3
  ) +
  scale_fill_manual(values = custom_colors) +
  labs(x = "", y = "", fill = "Diurnality") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "top"
  )

gb <- ggplot_build(p)
panel_data <- gb$data[[1]]

hex_to_name <- setNames(names(custom_colors), custom_colors)

label_df <- panel_data %>%
  mutate(
    Sleep_Timing = hex_to_name[fill],
    Number_of_sleep_times_per_day = x__Number_of_sleep_times_per_day,
    x = (xmin + xmax) / 2,
    y = (ymin + ymax) / 2
  )

label_info <- timing_times_Meta %>%
  count(Number_of_sleep_times_per_day, Sleep_Timing) %>%
  group_by(Number_of_sleep_times_per_day) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

label_df <- label_df %>%
  left_join(label_info, by = c("Number_of_sleep_times_per_day", "Sleep_Timing")) %>%
  mutate(label = paste0(n, " Species\n(", round(prop * 100, 1), "%)"))

times_timing_mosicPlot <- p + geom_text(
  data = label_df %>% filter(!is.na(n)), 
  aes(x = x, y = y, label = label),
  inherit.aes = FALSE,
  size = 4.2,
  fontface = "bold",
  color = "black"
)

times_timing_mosicPlot


############################################################
# End of script
############################################################
