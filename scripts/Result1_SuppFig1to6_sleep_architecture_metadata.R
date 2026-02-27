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
    Taxonomy_Class = Class
  ) %>%
  arrange(desc(NREM_ratio)) %>%
  mutate(
    Species_name_ensembl = fct_inorder(Species_name_ensembl)
  )


avg_nrem_all <- mean(nrem_meta$NREM_ratio, na.rm = TRUE)

avg_nrem_primates <- nrem_meta %>%
  filter(Order == "Primates") %>%
  summarise(m = mean(NREM_ratio, na.rm = TRUE)) %>%
  pull(m)

avg_nrem_aves <- nrem_meta %>%
  filter(Class == "Aves") %>%
  summarise(m = mean(NREM_ratio, na.rm = TRUE)) %>%
  pull(m)

class_colors <- c(
  "Mammalia" = "#1f4e79",
  "Aves" = "#e6b800",
  "Actinopteri" = "#bfbfbf",
  "Insecta" = "#e07a3f",
  "Lepidosauria" = "#5b7f75",
  "Sarcopterygii" = "#e8d8a9"
)


supfig2 <- ggplot(nrem_meta,
                  aes(x = Species_name_ensembl)) +
  
  # 100% scaffold background
  geom_col(aes(y = 1),
           fill = "#e6e6e6",
           width = 0.6) +
  
  # Colored bars up to actual NREM ratio
  geom_col(aes(y = NREM_ratio,
               fill = Taxonomy_Class),
           width = 0.6) +
  
  # Top points (legend suppressed)
  geom_point(aes(y = NREM_ratio,
                 color = Taxonomy_Class),
             size = 6,
             show.legend = FALSE) +
  
  # Global average
  geom_hline(yintercept = avg_nrem_all,
             linetype = "dashed",
             color = "red",
             linewidth = 1) +
  
  # Primate average
  geom_hline(yintercept = avg_nrem_primates,
             linetype = "dashed",
             color = "#f1b321",
             linewidth = 1) +
  
  # Aves average
  geom_hline(yintercept = avg_nrem_aves,
             linetype = "dashed",
             color = "#194a7a",
             linewidth = 1) +
  
  # Annotations
  annotate("text",
           x = length(unique(nrem_meta$Species_name_ensembl)) * 0.85,
           y = avg_nrem_all + 0.03,
           label = paste0("All species mean = ",
                          round(avg_nrem_all, 2), "%"),
           color = "red",
           size = 4) +
  
  annotate("text",
           x = length(unique(nrem_meta$Species_name_ensembl)) * 0.85,
           y = avg_nrem_primates + 0.03,
           label = paste0("Primates mean = ",
                          round(avg_nrem_primates, 2), "%"),
           color = "#f1b321",
           size = 4) +
  
  annotate("text",
           x = length(unique(nrem_meta$Species_name_ensembl)) * 0.85,
           y = avg_nrem_aves + 0.03,
           label = paste0("Aves mean = ",
                          round(avg_nrem_aves, 2), "%"),
           color = "#194a7a",
           size = 4) +
  
  scale_fill_manual(values = class_colors) +
  scale_color_manual(values = class_colors) +
  
  labs(
    y = "NREM Sleep Ratio (%)",
    x = NULL,
    fill = "Taxonomy Class"
  ) +
  
  coord_cartesian(ylim = c(0, 1)) +
  
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "top"
  )


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

bar_data <- sleep_freq_meta %>%
  count(Class, Order, Number_of_sleep_times_per_day) %>%
  rename(Count = n)


bar_data <- sleep_freq_meta %>%
  count(Class, Order, Number_of_sleep_times_per_day) %>%
  rename(Count = n) %>%
  mutate(Species_Label = ifelse(Count == 1,
                                sleep_freq_meta$Species_name_ensembl[match(paste(Class, Order, Number_of_sleep_times_per_day), 
                                                                           paste(sleep_freq_meta$Class, sleep_freq_meta$Order, sleep_freq_meta$Number_of_sleep_times_per_day))],
                                NA))


p <- ggplot(bar_data, aes(y = reorder(Order, Count), 
                          x = ifelse(Number_of_sleep_times_per_day == "Once", -Count, Count), 
                          fill = Number_of_sleep_times_per_day)) +
  geom_bar(stat = "identity", width = 0.4) +  
  scale_fill_manual(values = c("Once" = "#0072B2", "More than twice" = "#E69F00")) +
  scale_x_continuous(labels = abs) +  # X축 값을 절대값으로 표시
  theme_minimal() +
  theme(axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)) +
  labs(title = "Distribution of Sleep Frequency Across Taxonomy Order Level",
       x = "Number of Species",
       fill = "Sleep Frequency") +
  geom_point(data = bar_data, 
             aes(x = ifelse(Number_of_sleep_times_per_day == "Once", -Count, Count), 
                 y = reorder(Order, Count), 
                 fill = Number_of_sleep_times_per_day),
             shape = 21, size = 9, stroke = 0, inherit.aes = FALSE)+
  geom_text(aes(label = Count, 
                hjust = ifelse(Number_of_sleep_times_per_day == "Once", 0.6, 0.4)),
            fontface = "bold",
            size = 4) +  # 한 번 자는 경우 왼쪽 정렬
  facet_grid(Class ~ ., scales = "free_y", space = "free_y")

bar_data_labels <- bar_data %>% filter(!is.na(Species_Label))
bar_data_labels <- bar_data_labels %>%
  mutate(
    x_pos = ifelse(Number_of_sleep_times_per_day == "Once", -Count - 0.5, Count + 0.5),
    hjust_pos = ifelse(Number_of_sleep_times_per_day == "Once", 1, 0)
  )


supfig3 <- p + geom_text(data = bar_data_labels, 
                    aes(x = x_pos,  # 🔹 X 위치 조정
                        y = reorder(Order, Count), 
                        label = Species_Label), 
                    hjust = ifelse(bar_data_labels$Number_of_sleep_times_per_day == "Once", 1, 0),  # 🔹 aes() 밖에서 설정
                    vjust = 0.5, 
                    size = 4, 
                    fontface = "bold",
                    color = "black") +
  theme(strip.text.y = element_text(size = 11, angle = 0, hjust = 0.5, vjust = 0.5)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 2)) +
  theme(strip.background = element_rect(fill = "gray90", color = "black",size = 1)) +
  geom_vline(xintercept = 0, color = "black", size = 1.2)

supfig3

############################################################
# Supplementary Figure 4A
# Sleep timing patterns across all species
############################################################
sleep_timing_meta <- species_metadata %>%
  filter(!is.na(Sleep_timing_per_day)) %>%
  mutate(
    Sleep_timing_per_day = factor(
      Sleep_timing_per_day,
      levels = c("Sleep at night",
                 "Sleep at daytime",
                 "Sleep at anytime")
    )
  )

# Count and compute percentages
timing_summary_all <- sleep_timing_meta %>%
  count(Sleep_timing_per_day) %>%
  mutate(
    percentage = n / sum(n) * 100,
    label = paste0(n, " (", round(percentage, 1), "%)")
  )

# Custom color palette (match figure style)
timing_colors <- c(
  "Sleep at night"   = "#2b9bb3",
  "Sleep at daytime" = "#67b531",
  "Sleep at anytime" = "#d9d9d9"
)

supfig4_all <- ggplot(timing_summary_all,
                      aes(x = "", y = n,
                          fill = Sleep_timing_per_day)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 5) +
  scale_fill_manual(values = timing_colors) +
  labs(
    title = paste0("Sleep timing patterns for a total of ",
                   sum(timing_summary_all$n),
                   " species"),
    fill = "Sleep Timing per day"
  ) +
  theme_void(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

supfig4_all


############################################################
# Supplementary Figure 4B
# Sleep timing patterns in primates
############################################################
timing_summary_primates <- sleep_timing_meta %>%
  filter(Order == "Primates") %>%
  count(Sleep_timing_per_day) %>%
  mutate(
    percentage = n / sum(n) * 100,
    label = paste0(n, " (", round(percentage, 1), "%)")
  )

timing_colors_primate <- c(
  "Sleep at night"   = "#7a8c7a",   # muted green
  "Sleep at daytime" = "#d8d2a8"    # pale beige
)

supfig4_primates <- ggplot(timing_summary_primates,
                           aes(x = "", y = n,
                               fill = Sleep_timing_per_day)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 5) +
  
  scale_fill_manual(values = timing_colors_primate) +
  
  labs(
    title = "Sleep timing patterns in primates",
    fill = "Sleep Timing per day"
  ) +
  
  theme_void(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

supfig4_primates



############################################################
# Supplementary Figure 5
# Association between sleep timing and sleep frequency
# in primates (Mouse lemur highlighted)
############################################################
timing_times_Meta <- species_metadata[!is.na(species_metadata$Sleep_timing_per_day) & !is.na(species_metadata$Number_of_sleep_times_per_day), ]
timing_times_Meta$Sleep_Timing <- factor(timing_times_Meta$Sleep_timing_per_day)
timing_times_Meta$Number_of_sleep_times_per_day <- factor(timing_times_Meta$Number_of_sleep_times_per_day)
rownames(timing_times_Meta) <- NULL


CochranArmitageTest(table(timing_times_Meta$Number_of_sleep_times_per_day, timing_times_Meta$Sleep_Timing))

Primates <- timing_times_Meta[which(timing_times_Meta$Order =="Primates"),]

set.seed(123)  

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
    x = "Sleep frequency",
    y = "Sleep timing",
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
  labs(x = "", y = "", fill = "Sleep timing") +
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
