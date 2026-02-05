# Combined amyloid/receptor analysis with percent area and intensity per cell
# - Reads all .xlsx files in the base directory (multiple workbooks, multiple sheets)
# - Expects columns: %Area, IntDen, Mean, Area, RawIntDen, MinThr, MaxThr, cell number
# - Computes per-cell percent area and per-cell integrated density
# - Produces summaries and plots using median/IQR + bootstrap 95% CI (recommended)

library(tidyverse)
library(readxl)

# -------- paths and files --------
setwd("~/Library/CloudStorage/OneDrive-Uppsalauniversitet/Master thesis")
# Specify which Excel files to read (adjust filenames as needed)
excel_files <- c("E2 exp. Receptor ERalpha.xlsx")  # Add more files like: c("file1.xlsx", "file2.xlsx")

# Define treatment order (extend as needed)
treatment_order <- c(
  "3d",
  "3d + Abeta",
  "3d+4",
  "3d+4 + Abeta",
  "3d 100 nM E2",
  "3d 100 nM E2 + Abeta",
  "3d+4 100 nM E2",
  "3d+4 100 nM E2 + Abeta"
)

# -------- helpers --------
median_ci <- function(x, reps = 2000, conf = 0.95) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(NA_real_, NA_real_))
  q <- (1 - conf) / 2
  boots <- replicate(reps, median(sample(x, length(x), replace = TRUE)))
  quantile(boots, c(q, 1 - q), na.rm = TRUE)
}

# Standardize cell-count column names
cell_candidates <- c("cell number", "Cell number", "cell_number", "Cells", "cells", "N_cells", "Count", "count")

# -------- load all files/sheets --------
read_one_file <- function(f) {
  sheets <- excel_sheets(f)
  # Skip sheets named "Results" (ImageJ output) and other non-data sheets
  skip_sheets <- c("Results", "Summary")
  sheets <- sheets[!(sheets %in% skip_sheets)]
  
  map_df(sheets, function(sh) {
    df <- read_excel(f, sheet = sh)
    cell_col <- intersect(cell_candidates, names(df))[1]
    if (is.na(cell_col)) {
      warning(paste("Skipping", basename(f), "sheet", sh, "- missing cell count column"))
      return(NULL)
    }

    df %>%
      rename(
        Percent.Area = `%Area`,
        Cells = all_of(cell_col)
      ) %>%
      mutate(
        Treatment = str_trim(sh),
        SourceFile = basename(f)
      )
  })
}

data_raw <- map_df(excel_files, read_one_file)

# Keep only expected treatments (warn on others)
unknown_treatments <- setdiff(unique(data_raw$Treatment), treatment_order)
if (length(unknown_treatments) > 0) {
  warning("Unexpected treatments present: ", paste(unknown_treatments, collapse = ", "))
  treatment_order <- c(treatment_order, setdiff(unknown_treatments, treatment_order))
}

data <- data_raw %>%
  mutate(Treatment = factor(Treatment, levels = treatment_order))

# -------- per-cell metrics --------
# Check if CS Date column exists
has_cs_date <- "CS Date" %in% names(data)

data_normalized <- data %>%
  mutate(
    PercentArea_per_cell = if_else(Cells > 0, Percent.Area / Cells, NA_real_),
    IntDen_per_cell = if_else(Cells > 0, IntDen / Cells, NA_real_),
    Mean_intensity = Mean
  )

# Add CS_Date column - use real column if it exists, otherwise create a dummy
if (has_cs_date) {
  data_normalized <- data_normalized %>%
    mutate(CS_Date = as.factor(`CS Date`))
} else {
  data_normalized <- data_normalized %>%
    mutate(CS_Date = as.factor(rep("All", nrow(.))))
}

# -------- filter outliers (optional) --------
# Remove extreme outliers using IQR method (points beyond 1.5*IQR from Q1/Q3)
filter_outliers <- function(x) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower_bound <- q1 - 1.5 * iqr
  upper_bound <- q3 + 1.5 * iqr
  x >= lower_bound & x <= upper_bound
}

data_normalized <- data_normalized %>%
  mutate(
    keep_IntDen = filter_outliers(IntDen_per_cell),
    keep_PercentArea = filter_outliers(PercentArea_per_cell)
  ) %>%
  # Exclude extreme outliers from all analyses
  filter(keep_IntDen & keep_PercentArea) %>%
  select(-keep_IntDen, -keep_PercentArea)

# Data for plotting (same as data_normalized, no further filtering needed)
data_plot <- data_normalized

# -------- summaries --------
pilot_summary_cells <- data_normalized %>%
  group_by(Treatment) %>%
  summarise(
    N_images = n(),
    N_cells_total = sum(Cells, na.rm = TRUE),
    Mean_PercentArea_per_cell = mean(PercentArea_per_cell, na.rm = TRUE),
    SD_PercentArea_per_cell = sd(PercentArea_per_cell, na.rm = TRUE),
    Mean_IntDen_per_cell = mean(IntDen_per_cell, na.rm = TRUE),
    SD_IntDen_per_cell = sd(IntDen_per_cell, na.rm = TRUE),
    .groups = "drop"
  )

pilot_summary_cells_median <- data_normalized %>%
  group_by(Treatment) %>%
  summarise(
    N_images = n(),
    Median_PercentArea_per_cell = median(PercentArea_per_cell, na.rm = TRUE),
    Q1_PercentArea = quantile(PercentArea_per_cell, 0.25, na.rm = TRUE),
    Q3_PercentArea = quantile(PercentArea_per_cell, 0.75, na.rm = TRUE),
    CI_PercentArea = list(median_ci(PercentArea_per_cell)),
    CI_PercentArea_low = map_dbl(CI_PercentArea, 1),
    CI_PercentArea_high = map_dbl(CI_PercentArea, 2),
    Median_IntDen_per_cell = median(IntDen_per_cell, na.rm = TRUE),
    Q1_IntDen = quantile(IntDen_per_cell, 0.25, na.rm = TRUE),
    Q3_IntDen = quantile(IntDen_per_cell, 0.75, na.rm = TRUE),
    CI_IntDen = list(median_ci(IntDen_per_cell)),
    CI_IntDen_low = map_dbl(CI_IntDen, 1),
    CI_IntDen_high = map_dbl(CI_IntDen, 2),
    .groups = "drop"
  ) %>%
  select(-CI_PercentArea, -CI_IntDen)

# -------- plots --------
# Percent area per nucleus (median/IQR + CI) with raw points (outliers excluded)
plot_percent_area <- ggplot() +
  geom_jitter(data = data_plot,
              aes(x = Treatment, y = PercentArea_per_cell, color = CS_Date),
              width = 0.15, size = 2, alpha = 0.7) +
  geom_linerange(data = pilot_summary_cells_median,
                 aes(x = Treatment, ymin = Q1_PercentArea, ymax = Q3_PercentArea),
                 size = 1.1, color = "black") +
  geom_pointrange(data = pilot_summary_cells_median,
                  aes(x = Treatment,
                      y = Median_PercentArea_per_cell,
                      ymin = CI_PercentArea_low,
                      ymax = CI_PercentArea_high),
                  shape = 21, fill = "white", size = 0.8, stroke = 0.8,
                  color = "black") +
  ylab("% area per nucleus") +
  xlab("") +
  scale_color_viridis_d(name = "CS Date") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "right")

# Intensity per nucleus (median/IQR + CI) with raw points (outliers excluded)
plot_intensity <- ggplot() +
  geom_jitter(data = data_plot,
              aes(x = Treatment, y = IntDen_per_cell, color = CS_Date),
              width = 0.15, size = 2, alpha = 0.7) +
  geom_linerange(data = pilot_summary_cells_median,
                 aes(x = Treatment, ymin = Q1_IntDen, ymax = Q3_IntDen),
                 size = 1.1, color = "black") +
  geom_pointrange(data = pilot_summary_cells_median,
                  aes(x = Treatment,
                      y = Median_IntDen_per_cell,
                      ymin = CI_IntDen_low,
                      ymax = CI_IntDen_high),
                  shape = 21, fill = "white", size = 0.8, stroke = 0.8,
                  color = "black") +
  ylab("Integrated density per nucleus") +
  xlab("") +
  scale_color_viridis_d(name = "CS Date") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "right")

# Bar charts with median and bootstrap CI
plot_percent_area_bar <- ggplot(pilot_summary_cells_median,
                                aes(x = Treatment, y = Median_PercentArea_per_cell)) +
  geom_col(fill = "steelblue", width = 0.6) +
  geom_errorbar(aes(ymin = CI_PercentArea_low, ymax = CI_PercentArea_high),
                width = 0.2, color = "black") +
  ylab("% area per nucleus (median)") +
  xlab("") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

plot_intensity_bar <- ggplot(pilot_summary_cells_median,
                             aes(x = Treatment, y = Median_IntDen_per_cell)) +
  geom_col(fill = "gray70", color = "gray40", width = 0.6) +
  geom_errorbar(aes(ymin = CI_IntDen_low, ymax = CI_IntDen_high),
                width = 0.2, color = "black") +
  ylab("Integrated density per nucleus (median)") +
  xlab("") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# Confounder check: cell count vs. intensity per cell (outliers excluded)
plot_cells_vs_intensity <- ggplot(data_plot,
       aes(x = Cells, y = IntDen_per_cell, color = Treatment)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.5) +
  facet_wrap(~ Treatment, scales = "free") +
  ylab("Integrated density per nucleus") +
  xlab("Nuclei per image") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

# Print summaries and plots when sourced
print(pilot_summary_cells)
print(pilot_summary_cells_median)
print(plot_percent_area)
print(plot_intensity)
print(plot_percent_area_bar)
print(plot_intensity_bar)
print(plot_cells_vs_intensity)
