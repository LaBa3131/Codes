# Colocalization analysis: Amyloid vs Estrogen Receptor overlap
# - Reads amyloid and receptor Excel files
# - Aligns data by treatment
# - Computes Pearson correlation and overlap metrics
# - Produces scatter plots and summary statistics

library(tidyverse)
library(readxl)

# -------- setup --------
setwd("~/Library/CloudStorage/OneDrive-Uppsalauniversitet/Master thesis")

# Check available files
all_files <- list.files(".", pattern = "\\.xlsx$")
print("Available Excel files in directory:")
print(all_files)

# Specify files (update these filenames if they differ)
amyloid_file <- "E2 exp. Abeta.xlsx"
receptor_file <- "E2 exp. Receptor ERalpha.xlsx"

# Verify files exist
if (!file.exists(amyloid_file)) stop(paste("Amyloid file not found:", amyloid_file))
if (!file.exists(receptor_file)) stop(paste("Receptor file not found:", receptor_file))

# Standardize cell-count column names
cell_candidates <- c("cell number", "Cell number", "cell_number", "Cells", "cells", "N_cells", "Count", "count")

# -------- load and prepare data --------
read_and_prep <- function(f, channel_name) {
  sheets <- excel_sheets(f)
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
        PercentArea = `%Area`,
        Cells = all_of(cell_col)
      ) %>%
      select(PercentArea, IntDen, Mean, Cells) %>%
      mutate(
        Treatment = str_trim(sh),
        Channel = channel_name,
        Image_ID = row_number()  # temporary ID within treatment
      )
  })
}

amyloid <- read_and_prep(amyloid_file, "Amyloid")
receptor <- read_and_prep(receptor_file, "Receptor_ERalpha")

# -------- align data --------
# Merge by treatment and image number (assumes same image order in both files)
colocalized <- amyloid %>%
  rename_with(~ paste0("Amyloid_", .), -c(Treatment, Image_ID)) %>%
  inner_join(
    receptor %>%
      rename_with(~ paste0("Receptor_", .), -c(Treatment, Image_ID)),
    by = c("Treatment", "Image_ID")
  ) %>%
  # Compute per-cell normalized metrics
  mutate(
    Amyloid_per_cell = if_else(Amyloid_Cells > 0, Amyloid_PercentArea / Amyloid_Cells, NA_real_),
    Receptor_per_cell = if_else(Receptor_Cells > 0, Receptor_PercentArea / Receptor_Cells, NA_real_),
    Amyloid_IntDen_per_cell = if_else(Amyloid_Cells > 0, Amyloid_IntDen / Amyloid_Cells, NA_real_),
    Receptor_IntDen_per_cell = if_else(Receptor_Cells > 0, Receptor_IntDen / Receptor_Cells, NA_real_)
  )

# -------- filter outliers --------
filter_outliers <- function(x) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower_bound <- q1 - 1.5 * iqr
  upper_bound <- q3 + 1.5 * iqr
  x >= lower_bound & x <= upper_bound
}

colocalized <- colocalized %>%
  mutate(
    keep_Amyloid = filter_outliers(Amyloid_IntDen_per_cell),
    keep_Receptor = filter_outliers(Receptor_IntDen_per_cell)
  ) %>%
  filter(keep_Amyloid & keep_Receptor) %>%
  select(-keep_Amyloid, -keep_Receptor)

# -------- correlation analysis --------
correlation_summary <- colocalized %>%
  group_by(Treatment) %>%
  summarise(
    N_images = n(),
    # IntDen per cell correlations
    Corr_IntDen = cor(Amyloid_IntDen_per_cell, Receptor_IntDen_per_cell, use = "complete.obs"),
    # %Area per cell correlations
    Corr_PercentArea = cor(Amyloid_per_cell, Receptor_per_cell, use = "complete.obs"),
    # Manders' overlap coefficient (fraction of signal overlap)
    # M1 = fraction of Amyloid that colocalizes with Receptor
    Manders_M1 = {
      amyl <- Amyloid_IntDen_per_cell[is.finite(Amyloid_IntDen_per_cell) & is.finite(Receptor_IntDen_per_cell)]
      recep <- Receptor_IntDen_per_cell[is.finite(Amyloid_IntDen_per_cell) & is.finite(Receptor_IntDen_per_cell)]
      sum(amyl[recep > quantile(recep, 0.25, na.rm = TRUE)]) / sum(amyl)
    },
    # M2 = fraction of Receptor that colocalizes with Amyloid
    Manders_M2 = {
      amyl <- Amyloid_IntDen_per_cell[is.finite(Amyloid_IntDen_per_cell) & is.finite(Receptor_IntDen_per_cell)]
      recep <- Receptor_IntDen_per_cell[is.finite(Amyloid_IntDen_per_cell) & is.finite(Receptor_IntDen_per_cell)]
      sum(recep[amyl > quantile(amyl, 0.25, na.rm = TRUE)]) / sum(recep)
    },
    .groups = "drop"
  )

print("=== Correlation & Manders' Overlap Summary ===")
print(correlation_summary)

# -------- scatter plots --------
# Prepare annotation data with statistics for each treatment
stats_for_plots <- colocalized %>%
  group_by(Treatment) %>%
  summarise(
    Corr_IntDen = cor(Amyloid_IntDen_per_cell, Receptor_IntDen_per_cell, use = "complete.obs"),
    Corr_PercentArea = cor(Amyloid_per_cell, Receptor_per_cell, use = "complete.obs"),
    N = n(),
    # Manders for intden
    Manders_M1_text = {
      amyl <- Amyloid_IntDen_per_cell[is.finite(Amyloid_IntDen_per_cell) & is.finite(Receptor_IntDen_per_cell)]
      recep <- Receptor_IntDen_per_cell[is.finite(Amyloid_IntDen_per_cell) & is.finite(Receptor_IntDen_per_cell)]
      round(sum(amyl[recep > quantile(recep, 0.25, na.rm = TRUE)]) / sum(amyl), 3)
    },
    Manders_M2_text = {
      amyl <- Amyloid_IntDen_per_cell[is.finite(Amyloid_IntDen_per_cell) & is.finite(Receptor_IntDen_per_cell)]
      recep <- Receptor_IntDen_per_cell[is.finite(Amyloid_IntDen_per_cell) & is.finite(Receptor_IntDen_per_cell)]
      round(sum(recep[amyl > quantile(amyl, 0.25, na.rm = TRUE)]) / sum(recep), 3)
    },
    .groups = "drop"
  )

# Prepare Manders data for visualization (long format)
manders_plot_data <- stats_for_plots %>%
  select(Treatment, Manders_M1_text, Manders_M2_text) %>%
  rename(M1 = Manders_M1_text, M2 = Manders_M2_text) %>%
  pivot_longer(cols = c(M1, M2), names_to = "Coefficient", values_to = "Value")

# Manders coefficient visualization
plot_manders <- ggplot(manders_plot_data,
       aes(x = Treatment, y = Value, fill = Coefficient)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = round(Value, 3)),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3.5) +
  ylab("Manders Coefficient (0-1)") +
  xlab("") +
  ylim(0, 1) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 20, hjust = 1))

# Plot 1: IntDen per cell correlation with Manders coefficients
plot_coloc_intden <- ggplot(colocalized,
       aes(x = Amyloid_IntDen_per_cell, y = Receptor_IntDen_per_cell, color = Treatment)) +
  geom_point(size = 2.5, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  facet_wrap(~ Treatment, scales = "free") +
  # Add statistics as text
  geom_text(data = stats_for_plots,
            aes(x = -Inf, y = Inf, label = paste0("r = ", round(Corr_IntDen, 3), 
                                                     "\nM1 = ", Manders_M1_text,
                                                     "\nM2 = ", Manders_M2_text,
                                                     "\nN = ", N)),
            vjust = 1.1, hjust = -0.05, size = 3.5, color = "black", parse = FALSE) +
  xlab("Amyloid (Integrated density per nucleus)") +
  ylab("ER-alpha (Integrated density per nucleus)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

# Plot 2: %Area per cell correlation
plot_coloc_area <- ggplot(colocalized,
       aes(x = Amyloid_per_cell, y = Receptor_per_cell, color = Treatment)) +
  geom_point(size = 2.5, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  facet_wrap(~ Treatment, scales = "free") +
  # Add correlation and N as text
  geom_text(data = stats_for_plots,
            aes(x = -Inf, y = Inf, label = paste0("r = ", round(Corr_PercentArea, 3), 
                                                     "\nN = ", N)),
            vjust = 1.1, hjust = -0.05, size = 3.5, color = "black", parse = FALSE) +
  xlab("Amyloid (% area per nucleus)") +
  ylab("ER-alpha (% area per nucleus)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

# Plot 3: Overall correlation (all treatments combined)
overall_stats <- data.frame(
  Corr_IntDen = cor(colocalized$Amyloid_IntDen_per_cell,
                    colocalized$Receptor_IntDen_per_cell,
                    use = "complete.obs"),
  N = nrow(colocalized)
)

plot_coloc_overall <- ggplot(colocalized,
       aes(x = Amyloid_IntDen_per_cell, y = Receptor_IntDen_per_cell, color = Treatment)) +
  geom_point(size = 2.5, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2, aes(group = 1), color = "black") +
  # Add overall statistics
  geom_text(data = overall_stats,
            aes(x = -Inf, y = Inf, label = paste0("Overall r = ", round(Corr_IntDen, 3),
                                                     "\nN = ", N)),
            vjust = 1.1, hjust = -0.05, size = 4, color = "black", parse = FALSE, inherit.aes = FALSE) +
  xlab("Amyloid (Integrated density per nucleus)") +
  ylab("ER-alpha (Integrated density per nucleus)") +
  theme_classic(base_size = 12)

# Print plots
print(plot_coloc_intden)
print(plot_coloc_area)
print(plot_coloc_overall)
print(plot_manders)

# -------- summary statistics --------
cat("\n=== Overall Colocalization Statistics ===\n")
overall_corr_intden <- cor(colocalized$Amyloid_IntDen_per_cell,
                           colocalized$Receptor_IntDen_per_cell,
                           use = "complete.obs")
overall_corr_area <- cor(colocalized$Amyloid_per_cell,
                         colocalized$Receptor_per_cell,
                         use = "complete.obs")
cat("Pearson r (IntDen per cell):", round(overall_corr_intden, 3), "\n")
cat("Pearson r (%Area per cell):", round(overall_corr_area, 3), "\n")
