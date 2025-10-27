# ------------------------------------------------------------
# Setup: Load libraries and set working directory
# ------------------------------------------------------------
setwd("~/Library/CloudStorage/OneDrive-Uppsalauniversitet/Lab LEC-8D3/ExVivo")

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(janitor)
library(gridExtra)

# ------------------------------------------------------------
# Antibody biodistribution plots ("wizard data" sheet)
# ------------------------------------------------------------
wizard_data <- read_excel("250929 Ex vivo Lec-8D3-pHrodo.xlsx",
                          sheet = "wizard data",
                          skip = 2) %>%
  janitor::clean_names()

tissues_keep <- c(
  "Term. Blood (9h)", "Pellet (9h)", "Plasma (9h)", "RH (9h)",
  "LH PFA (9h)", "Spleen (9h)", "Liver (9h)", "Thyroid (9h)",
  "Term. Blood (24h)", "Pellet (24h)", "Plasma (24h)", "RH (24h)",
  "LH PFA (24h)", "Spleen (24h)", "Liver (24h)", "Thyroid (24h)"
)

data_filt <- wizard_data %>%
  filter(tissue %in% tissues_keep) %>%
  mutate(
    timepoint = ifelse(grepl("\\(9h\\)", tissue), "9h", "24h"),
    organ = gsub(" \\(9h\\)| \\(24h\\)", "", tissue),
    antibody = gsub(".*_(ab[0-9]+).*", "\\1", tolower(sample)),
    percent_id_g = percent_id_g * 100   # multiply values by 100
  )

organ_order <- c("Term. Blood", "Pellet", "Plasma", "RH", "LH PFA", "Spleen", "Liver", "Thyroid")
data_filt$organ <- factor(data_filt$organ, levels = organ_order)
data_filt$timepoint <- factor(data_filt$timepoint, levels = c("9h", "24h"))

timepoint_colors <- c("9h" = "#FC8D62", "24h" = "#E34A33")

# ------------------------------------------------------------
# Plot each antibody with bars + SD + centered black points
# ------------------------------------------------------------
for (ab in unique(data_filt$antibody)) {
  data_ab <- data_filt %>% filter(antibody == ab)
  
  summary_data <- data_ab %>%
    group_by(organ, timepoint) %>%
    summarise(
      mean_percent_id_g = mean(percent_id_g, na.rm = TRUE),
      sd_percent_id_g = sd(percent_id_g, na.rm = TRUE),
      .groups = "drop"
    )
  
  summary_data$timepoint <- droplevels(summary_data$timepoint)
  data_ab$timepoint <- droplevels(data_ab$timepoint)
  
  active_levels <- levels(summary_data$timepoint)
  colors_fill <- timepoint_colors[active_levels]
  
  p <- ggplot(summary_data, aes(x = organ, y = mean_percent_id_g, fill = timepoint)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = mean_percent_id_g - sd_percent_id_g,
                      ymax = mean_percent_id_g + sd_percent_id_g),
                  width = 0.25, position = position_dodge(width = 0.8)) +
    geom_point(data = data_ab, aes(x = organ, y = percent_id_g, group = timepoint),
               position = position_dodge(width = 0.8),
               size = 2, alpha = 0.7, color = "black", inherit.aes = FALSE) +
    scale_fill_manual(values = colors_fill) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)) +  # fix axis
    labs(title = paste("%ID/g organ"),
         x = "", y = "%ID/g organ", fill = "timepoint") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(paste0("percentIDg_barplot_", ab, ".png"), plot = p, width = 10, height = 6, dpi = 300)
}

# ------------------------------------------------------------
# Injection and brain data plots ("Injections and brain" sheet)
# ------------------------------------------------------------
inj_data <- read_excel("250929 Ex vivo Lec-8D3-pHrodo.xlsx",
                       sheet = "Injections and brain",
                       skip = 1) %>%
  janitor::clean_names() %>%
  mutate(id = as.character(id),
         experiment = as.character(experiment))

desired_ids <- c("Bi30", "Bi32", "Bn10", "Bi31", "Bn11")
desired_experiments <- c("[125I]mAb158-scFv8D3-pHrodo", "[125I]Lec-LALA-8D3-pHrodo")

inj_filtered <- inj_data %>%
  filter(id %in% desired_ids, experiment %in% desired_experiments) %>%
  mutate(group = case_when(
    id == "Bi30" & experiment == "[125I]mAb158-scFv8D3-pHrodo" ~ "RmAB158-8D3 9h",
    id %in% c("Bi32", "Bn10") & experiment == "[125I]Lec-LALA-8D3-pHrodo" ~ "Lec8D3 9h",
    id %in% c("Bi31", "Bn11") & experiment == "[125I]Lec-LALA-8D3-pHrodo" ~ "Lec8D3 24h",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group)) %>%
  mutate(percent_id_g_rh = percent_id_g_rh * 100,
         percent_id_g_lh = percent_id_g_lh * 100)

group_colors <- c("RmAB158-8D3 9h" = "#E41A1C",
                  "Lec8D3 9h"  = "#FC4E2A",
                  "Lec8D3 24h" = "#FCA05E")
inj_filtered$group <- factor(inj_filtered$group, levels = names(group_colors))

# ------------------------------------------------------------
# Function to plot RH/LH metrics with SD + black points
# ------------------------------------------------------------
plot_metric <- function(data, value_col, y_label) {
  summary_data <- data %>%
    group_by(group) %>%
    summarise(
      mean_val = mean(.data[[value_col]], na.rm = TRUE),
      sd_val   = sd(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    right_join(data.frame(group = names(group_colors)), by = "group") %>%
    mutate(
      mean_val = ifelse(is.na(mean_val), 0, mean_val),
      sd_val = ifelse(is.na(sd_val), 0, sd_val),
      group = factor(group, levels = names(group_colors))
    )
  
  p <- ggplot(summary_data, aes(x = group, y = mean_val, fill = group)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val), width = 0.2) +
    geom_point(data = data %>% filter(!is.na(.data[[value_col]])),
               aes(x = group, y = .data[[value_col]]),
               position = position_dodge(width = 0.7),
               size = 2, alpha = 0.7, color = "black") +
    scale_fill_manual(values = group_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)) +  # fix axis
    labs(title = y_label, x = "group", y = y_label, fill = "group") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  return(p)
}

# RH plots
p_rh_idg   <- plot_metric(inj_filtered, "percent_id_g_rh", "%ID/g RH")
p_rh_ratio <- plot_metric(inj_filtered, "rh_to_blood", "RH-to-Blood Ratio")
grid_rh <- gridExtra::grid.arrange(p_rh_idg, p_rh_ratio, ncol = 2)
ggsave("RH_metrics_combined.png", grid_rh, width = 15, height = 6, dpi = 300)

# LH plots
p_lh_idg   <- plot_metric(inj_filtered, "percent_id_g_lh", "%ID/g LH")
p_lh_ratio <- plot_metric(inj_filtered, "lh_to_blood", "LH-to-Blood Ratio")
grid_lh <- gridExtra::grid.arrange(p_lh_idg, p_lh_ratio, ncol = 2)
ggsave("LH_metrics_combined.png", grid_lh, width = 15, height = 6, dpi = 300)

# ------------------------------------------------------------
# Linear plot from "Biodistribution & blood" sheet (starting AB15)
# ------------------------------------------------------------
bio_data <- read_excel(
  "250929 Ex vivo Lec-8D3-pHrodo.xlsx",
  sheet = "Biodistribution & blood",
  range = "AB15:AD1000",  # adjust columns/rows if needed
  col_names = TRUE
)

# Rename columns to preserve proper case
names(bio_data) <- gsub("rmab158", "RmAB158-8D3", names(bio_data), ignore.case = TRUE)
names(bio_data) <- gsub("lec8d3", "Lec8D3", names(bio_data), ignore.case = TRUE)


bio_long <- bio_data %>%
  pivot_longer(cols = -`Time (h)`, names_to = "antibody", values_to = "value") %>%
  mutate(`Time (h)` = as.numeric(`Time (h)`),
         value = as.numeric(value) * 100) %>%  # multiply values by 100
  filter(!is.na(value))

p_linear <- ggplot(bio_long, aes(x = `Time (h)`, y = value, color = antibody)) +
  geom_line(size = 1.2) +
  geom_point(color = "black", size = 2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)) +  # fix axis
  labs(title = "Biodistribution",
       x = "Time (h)", y = "%ID/g blood") +
  theme_minimal(base_size = 14)

ggsave("Biodistribution_linear.png", p_linear, width = 8, height = 6, dpi = 300)

# Optional: Combine RH/LH and linear plot side by side
grid_all <- gridExtra::grid.arrange(p_rh_idg, p_rh_ratio, p_linear, ncol = 3)
ggsave("Biodistribution_RH_linear.png", grid_all, width = 18, height = 6, dpi = 300)

# ------------------------------------------------------------
# Completion message
# ------------------------------------------------------------
message("✅ All plots generated successfully, only axis scaling updated.")
