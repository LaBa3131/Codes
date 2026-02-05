############################################
# PILOTPROJEKT: Amyloid-Belastung (%Area)
# Ziel:
# Einlesen aller Behandlungen aus einer Excel-Datei,
# Mittelung von 20 technischen Replikaten pro Behandlung
# und deskriptiver Vergleich der Amyloid-Belastung
############################################

# Pakete laden
library(tidyverse)
library(readxl)

# Excel-Datei (Dateinamen ggf. anpassen)
setwd("~/Library/CloudStorage/OneDrive-Uppsalauniversitet/Master thesis")
file <- "E2 exp. Abeta.xlsx"

# Alle Sheets einlesen und Treatment hinzufügen
sheets <- excel_sheets(file)
print("Sheet names found:")
print(sheets)

data <- map_df(sheets, function(sheet) {
  read_excel(file, sheet = sheet) %>%
    mutate(Treatment = str_trim(sheet))
})

# Enforce plotting order (keep all treatments, order the known ones)
treatment_order <- c(
  "3d + Abeta",
  "3d+4 + Abeta",
  "3d 100 nM E2 + Abeta",
  "3d+4 100 nM E2 + Abeta"
)

# Diagnostic: show which treatments are present/missing
print("Unique treatments in data:")
print(unique(data$Treatment))
print("Treatments expected but not found:")
print(setdiff(treatment_order, unique(data$Treatment)))

data <- data %>%
  mutate(Treatment = factor(Treatment, levels = treatment_order))

# Rename %Area to an R-friendly column
data <- data %>%
  rename(Percent.Area = `%Area`)

# Qualitätskontrolle
summary(data$Percent.Area)

# Pilot-Auswertung:
# 20 Bilder = technische Replikate → 1 Wert pro Behandlung
pilot_summary <- data %>%
  group_by(Treatment) %>%
  summarise(
    Mean_PercentArea = mean(Percent.Area),
    SD = sd(Percent.Area),
    N_images = n(),
    .groups = "drop"
  )

print(pilot_summary)

# Plot (descriptive)
ggplot(pilot_summary,
       aes(x = Treatment, y = Mean_PercentArea)) +
  geom_col(fill = "gray75", width = 0.6) +
  geom_errorbar(
    aes(
      ymin = Mean_PercentArea - SD,
      ymax = Mean_PercentArea + SD
    ),
    width = 0.2
  ) +
  ylab("Amyloid burden (% area)") +
  xlab("") +
  theme_classic(base_size = 14)

# OPTIONAL: show individual images (illustrative only)
ggplot(data,
       aes(x = Treatment, y = Percent.Area)) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  ylab("Amyloid burden (% area)") +
  xlab("") +
  theme_classic(base_size = 14)

# IMPORTANT:
# No statistics (pilot, n = 1 section per treatment)

############################################
# Extension: %Area per nucleus
# + Scatter cell count vs. amyloid per cell
############################################

# Ensure a cell-count column exists (tolerates common spellings)
cell_candidates <- c("cell number", "Cell number", "cell_number", "Cells", "cells", "N_cells", "Count", "count")
cell_col <- intersect(cell_candidates, names(data))[1]
if (is.na(cell_col)) {
  stop("Cell count column is missing (e.g., 'cell number').")
}

# Compute %Area per nucleus (avoid division by zero)
# Check if CS Date column exists
has_cs_date <- "CS Date" %in% names(data)

data_normalized <- data %>%
  rename(Cells = all_of(cell_col)) %>%
  mutate(
    # Percent.Area is already a percentage; divide by cell count only
    PercentArea_per_cell = if_else(Cells > 0, Percent.Area / Cells, NA_real_)
  )

# Add CS_Date column - use real column if it exists, otherwise create a dummy
if (has_cs_date) {
  data_normalized <- data_normalized %>%
    mutate(CS_Date = as.factor(`CS Date`))
} else {
  data_normalized <- data_normalized %>%
    mutate(CS_Date = as.factor(rep("All", nrow(.))))
}

# Helper: bootstrap CI for the median (returns c(low, high))
median_ci <- function(x, reps = 2000, conf = 0.95) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(NA_real_, NA_real_))
  q <- (1 - conf) / 2
  boots <- replicate(reps, median(sample(x, length(x), replace = TRUE)))
  quantile(boots, c(q, 1 - q), na.rm = TRUE)
}

# Pilot summary (per nucleus)
pilot_summary_cells <- data_normalized %>%
  group_by(Treatment) %>%
  summarise(
    N_images = n(),
    N_cells_total = sum(Cells, na.rm = TRUE),
    Mean_PercentArea_per_cell = mean(PercentArea_per_cell, na.rm = TRUE),
    SD_PercentArea_per_cell = sd(PercentArea_per_cell, na.rm = TRUE),
    .groups = "drop"
  )

print(pilot_summary_cells)

# Median + IQR + bootstrap 95% CI (recommended for bounded/skewed % data)
pilot_summary_cells_median <- data_normalized %>%
  group_by(Treatment) %>%
  summarise(
    N_images = n(),
    Median_PercentArea_per_cell = median(PercentArea_per_cell, na.rm = TRUE),
    Q1 = quantile(PercentArea_per_cell, 0.25, na.rm = TRUE),
    Q3 = quantile(PercentArea_per_cell, 0.75, na.rm = TRUE),
    CI = list(median_ci(PercentArea_per_cell)),
    CI_low = map_dbl(CI, 1),
    CI_high = map_dbl(CI, 2),
    .groups = "drop"
  ) %>%
  select(-CI)

# Plot 1: %Area per nucleus (main result)
ggplot(pilot_summary_cells,
       aes(x = Treatment, y = Mean_PercentArea_per_cell)) +
  geom_col(fill = "steelblue", width = 0.6) +
  geom_errorbar(
    aes(
      ymin = Mean_PercentArea_per_cell - SD_PercentArea_per_cell,
      ymax = Mean_PercentArea_per_cell + SD_PercentArea_per_cell
    ),
    width = 0.2
  ) +
  ylab("Amyloid (% area per nucleus)") +
  xlab("") +
  theme_classic(base_size = 14)

# Plot 1b (preferred): median + IQR + bootstrap 95% CI with raw points
ggplot() +
  geom_jitter(data = data_normalized,
              aes(x = Treatment, y = PercentArea_per_cell, color = CS_Date),
              width = 0.15, size = 2, alpha = 0.7) +
  geom_linerange(data = pilot_summary_cells_median,
                 aes(x = Treatment, ymin = Q1, ymax = Q3),
                 size = 1.1, color = "black") +
  geom_pointrange(data = pilot_summary_cells_median,
                  aes(x = Treatment,
                      y = Median_PercentArea_per_cell,
                      ymin = CI_low,
                      ymax = CI_high),
                  shape = 21, fill = "white", size = 0.8, stroke = 0.8,
                  color = "black") +
  ylab("Amyloid (% area per nucleus)") +
  xlab("") +
  scale_color_viridis_d(name = "CS Date") +
  theme_classic(base_size = 14) +
  theme(legend.position = "right")

# Plot 2: cell count vs. amyloid per cell (confounder check)
ggplot(data_normalized,
       aes(x = Cells, y = PercentArea_per_cell, color = Treatment)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.5) +
  facet_wrap(~ Treatment, scales = "free") +
  ylab("Amyloid (% area per nucleus)") +
  xlab("Nuclei per image") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

############################################
# Extension: Intensity analysis per nucleus
############################################

# Ensure an intensity column exists (tolerates common spellings)
intensity_candidates <- c("Mean", "mean", "IntDen", "RawIntDen", "Intensity", "intensity", "Mean_Intensity")
intensity_col <- intersect(intensity_candidates, names(data))[1]
if (is.na(intensity_col)) {
  stop("Intensity column is missing (e.g., 'Mean', 'IntDen').")
}

# Compute intensity per nucleus
data_normalized <- data_normalized %>%
  rename(Intensity = all_of(intensity_col)) %>%
  mutate(
    Intensity_per_cell = if_else(Cells > 0, Intensity / Cells, NA_real_)
  )

# Summary statistics for intensity per nucleus
pilot_summary_intensity <- data_normalized %>%
  group_by(Treatment) %>%
  summarise(
    N_images = n(),
    N_cells_total = sum(Cells, na.rm = TRUE),
    Mean_Intensity_per_cell = mean(Intensity_per_cell, na.rm = TRUE),
    SD_Intensity_per_cell = sd(Intensity_per_cell, na.rm = TRUE),
    .groups = "drop"
  )

print(pilot_summary_intensity)

# Median + IQR + bootstrap 95% CI for intensity
pilot_summary_intensity_median <- data_normalized %>%
  group_by(Treatment) %>%
  summarise(
    N_images = n(),
    Median_Intensity_per_cell = median(Intensity_per_cell, na.rm = TRUE),
    Q1 = quantile(Intensity_per_cell, 0.25, na.rm = TRUE),
    Q3 = quantile(Intensity_per_cell, 0.75, na.rm = TRUE),
    CI = list(median_ci(Intensity_per_cell)),
    CI_low = map_dbl(CI, 1),
    CI_high = map_dbl(CI, 2),
    .groups = "drop"
  ) %>%
  select(-CI)

# Plot: Intensity per nucleus (mean + SD)
ggplot(pilot_summary_intensity,
       aes(x = Treatment, y = Mean_Intensity_per_cell)) +
  geom_col(fill = "darkorange", width = 0.6) +
  geom_errorbar(
    aes(
      ymin = Mean_Intensity_per_cell - SD_Intensity_per_cell,
      ymax = Mean_Intensity_per_cell + SD_Intensity_per_cell
    ),
    width = 0.2
  ) +
  ylab("Intensity per nucleus") +
  xlab("") +
  theme_classic(base_size = 14)

# Plot: Intensity per nucleus (median + IQR + bootstrap 95% CI with raw points)
ggplot() +
  geom_jitter(data = data_normalized,
              aes(x = Treatment, y = Intensity_per_cell, color = CS_Date),
              width = 0.15, size = 2, alpha = 0.7) +
  geom_linerange(data = pilot_summary_intensity_median,
                 aes(x = Treatment, ymin = Q1, ymax = Q3),
                 size = 1.1, color = "black") +
  geom_pointrange(data = pilot_summary_intensity_median,
                  aes(x = Treatment,
                      y = Median_Intensity_per_cell,
                      ymin = CI_low,
                      ymax = CI_high),
                  shape = 21, fill = "white", size = 0.8, stroke = 0.8,
                  color = "black") +
  ylab("Intensity per nucleus") +
  xlab("") +
  scale_color_viridis_d(name = "CS Date") +
  theme_classic(base_size = 14) +
  theme(legend.position = "right")

# Plot: cell count vs. intensity per cell (confounder check)
ggplot(data_normalized,
       aes(x = Cells, y = Intensity_per_cell, color = Treatment)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.5) +
  facet_wrap(~ Treatment, scales = "free") +
  ylab("Intensity per nucleus") +
  xlab("Nuclei per image") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")


