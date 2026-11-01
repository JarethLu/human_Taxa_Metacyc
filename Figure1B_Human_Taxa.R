# =============================================================================
# Figure 1B. Integrated Fn association and CRC depletion analysis
#
# Purpose
#   Screen 392 eligible species-level genome bins (SGBs) for taxa showing:
#   (1) a negative association with Fusobacterium nucleatum (Fn), and
#   (2) depletion in colorectal cancer (CRC) compared with controls.
#
# Screening criteria (all must be met)
#   Fn association:
#     - Pooled Spearman rho < 0
#     - Nominal P < 0.05
#     - Eligible cohorts >= 7
#   CRC depletion:
#     - Pooled Hedges' g < 0
#     - BH-adjusted FDR q < 0.20
#     - Eligible cohorts >= 7
#
# Expected final set:
#   - 14 highlighted taxa in total
#   - 1 Agathobacter rectalis (red)
#   - 13 other candidates (blue)
#   - 378 background SGBs (grey)
#
# Input:
#   data/input/Figure1B_plotting_data_14candidates_corrected.xlsx
#   Sheet: All_eligible_SGBs
#
# Outputs:
#   results/figures/Figure1B_14candidates_unlabeled.tiff
#   results/figures/Figure1B_14candidates_labeled.tiff
#
# Important:
#   This script starts from a processed SGB-level summary workbook. It does not
#   perform raw-read processing, MetaPhlAn profiling, cohort-level statistics,
#   or meta-analysis.
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Check and load required packages
# -----------------------------------------------------------------------------

required_packages <- c(
  "here",
  "readxl",
  "dplyr",
  "ggplot2",
  "ggrepel"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "The following R packages are missing: ",
      paste(missing_packages, collapse = ", "),
      "\nInstall them with:\ninstall.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "),
      "))"
    ),
    call. = FALSE
  )
}

library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)


# -----------------------------------------------------------------------------
# 1. Define project paths
# -----------------------------------------------------------------------------
# Open the .Rproj file in RStudio before running this script.
# here::here() resolves paths relative to the project root.

project_root <- here::here()

input_file <- here::here(
  "data",
  "input",
  "Figure1B_plotting_data_14candidates_corrected.xlsx"
)

output_dir <- here::here("results", "figures")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

message("Project root: ", project_root)
message("Input workbook: ", input_file)
message("Output directory: ", output_dir)

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Input workbook not found:\n",
      input_file,
      "\n\nPlace 'Figure1B_plotting_data_14candidates_corrected.xlsx' ",
      "in the data/input/ folder."
    ),
    call. = FALSE
  )
}


# -----------------------------------------------------------------------------
# 2. Read the SGB-level summary table
# -----------------------------------------------------------------------------

input_sheet <- "All_eligible_SGBs"

available_sheets <- readxl::excel_sheets(input_file)

if (!input_sheet %in% available_sheets) {
  stop(
    paste0(
      "The required worksheet '", input_sheet, "' was not found.\n",
      "Available worksheets: ",
      paste(available_sheets, collapse = ", ")
    ),
    call. = FALSE
  )
}

df <- readxl::read_excel(
  path = input_file,
  sheet = input_sheet
)

message("Rows in input table: ", nrow(df))


# -----------------------------------------------------------------------------
# 3. Validate required columns and data types
# -----------------------------------------------------------------------------

required_columns <- c(
  "Species",
  "SGB",
  "Fn pooled rho",
  "Fn P",
  "Fn FDR q",
  "Fn cohorts",
  "CRC Hedges g",
  "CRC P",
  "CRC FDR q",
  "CRC cohorts"
)

missing_columns <- setdiff(required_columns, names(df))

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "Missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}

numeric_columns <- c(
  "Fn pooled rho",
  "Fn P",
  "Fn FDR q",
  "Fn cohorts",
  "CRC Hedges g",
  "CRC P",
  "CRC FDR q",
  "CRC cohorts"
)

non_numeric_columns <- numeric_columns[
  !vapply(df[numeric_columns], is.numeric, logical(1))
]

if (length(non_numeric_columns) > 0) {
  stop(
    paste0(
      "The following columns must be numeric in the Excel workbook: ",
      paste(non_numeric_columns, collapse = ", "),
      ". Check the cell formats and values in the input file."
    ),
    call. = FALSE
  )
}

if (any(is.na(df$SGB))) {
  warning("The input table contains missing SGB identifiers.")
}


# -----------------------------------------------------------------------------
# 4. Prepare plotting data
# -----------------------------------------------------------------------------

plot_df <- df %>%
  transmute(
    Species = Species,
    SGB = as.character(SGB),
    pooled_rho = `Fn pooled rho`,
    Fn_P = `Fn P`,
    Fn_FDR = `Fn FDR q`,
    Fn_cohorts = `Fn cohorts`,
    pooled_g = `CRC Hedges g`,
    CRC_P = `CRC P`,
    CRC_FDR = `CRC FDR q`,
    CRC_cohorts = `CRC cohorts`
  ) %>%
  filter(
    !is.na(pooled_rho),
    !is.na(pooled_g)
  )

message(
  "SGBs with both Fn association and CRC depletion estimates: ",
  nrow(plot_df)
)

if (nrow(plot_df) != 392) {
  warning(
    paste0(
      "Expected 392 SGBs with both estimates, but found ",
      nrow(plot_df),
      ". Check the input workbook before interpreting the figure."
    )
  )
}


# -----------------------------------------------------------------------------
# 5. Select the 14 shortlisted taxa
# -----------------------------------------------------------------------------
# Note: Fn_FDR and CRC_P are retained in the plotting data for reporting, but
# they are not part of the screening conditions used for this figure.
# The Fn criterion uses nominal P < 0.05; the CRC criterion uses BH-FDR q < 0.20.

candidate14 <- plot_df %>%
  filter(
    # Fn association
    pooled_rho < 0,
    Fn_P < 0.05,
    Fn_cohorts >= 7,

    # CRC depletion
    pooled_g < 0,
    CRC_FDR < 0.20,
    CRC_cohorts >= 7
  ) %>%
  arrange(pooled_rho)

message(
  "Number of candidates passing the specified criteria: ",
  nrow(candidate14)
)

if (nrow(candidate14) != 14) {
  stop(
    paste0(
      "The current input selected ", nrow(candidate14),
      " taxa rather than the expected 14.\n",
      "Check the input workbook, column definitions, and screening criteria. ",
      "The script has stopped to avoid producing a figure that differs from ",
      "the intended analysis."
    ),
    call. = FALSE
  )
}


# -----------------------------------------------------------------------------
# 6. Assign publication-style taxon labels
# -----------------------------------------------------------------------------

candidate_names <- c(
  "SGB4933_group" = "Agathobacter rectalis",
  "SGB4290" = "Ruminococcus sp. NSJ-71",
  "SGB4750" = "Lachnospiraceae bacterium",
  "SGB4532" = "Anaerobutyricum hallii",
  "SGB5045" = "Eubacterium ventriosum",
  "SGB4925" = "Roseburia faecis",
  "SGB4394" = "Clostridia bacterium",
  "SGB5043" = "Eubacteriaceae bacterium",
  "SGB4751" = "Clostridium SGB4751",
  "SGB4951" = "Roseburia intestinalis",
  "SGB71759" = "GGB51441_SGB71759",
  "SGB4960" = "Eubacterium ramulus",
  "SGB8095" = "Streptococcus infantis",
  "SGB15332_group" = "Faecalibacterium prausnitzii"
)

unmapped_sgbs <- setdiff(candidate14$SGB, names(candidate_names))

if (length(unmapped_sgbs) > 0) {
  stop(
    paste0(
      "No publication-style label is defined for: ",
      paste(unmapped_sgbs, collapse = ", "),
      ". Update candidate_names before plotting."
    ),
    call. = FALSE
  )
}

candidate14 <- candidate14 %>%
  mutate(
    Taxon_label = unname(candidate_names[SGB])
  )

if (any(is.na(candidate14$Taxon_label))) {
  stop("At least one selected candidate has a missing taxon label.", call. = FALSE)
}

message("\nFinal 14 candidates:")
candidate14 %>%
  select(
    Taxon_label,
    SGB,
    pooled_rho,
    Fn_P,
    pooled_g,
    CRC_P,
    CRC_FDR
  ) %>%
  print(n = Inf)


# -----------------------------------------------------------------------------
# 7. Assign plotting groups
# -----------------------------------------------------------------------------

plot_df <- plot_df %>%
  mutate(
    Group = case_when(
      SGB == "SGB4933_group" ~ "Agathobacter rectalis",
      SGB %in% candidate14$SGB ~ "Other candidates",
      TRUE ~ "Background"
    )
  )

group_counts <- table(
  factor(
    plot_df$Group,
    levels = c(
      "Background",
      "Other candidates",
      "Agathobacter rectalis"
    )
  )
)

message("\nPlot group counts:")
print(group_counts)

expected_group_counts <- c(
  "Background" = 378,
  "Other candidates" = 13,
  "Agathobacter rectalis" = 1
)

if (!identical(as.integer(group_counts), as.integer(expected_group_counts))) {
  warning(
    paste0(
      "The observed plotting-group counts differ from the expected counts ",
      "(378 background, 13 other candidates, 1 Agathobacter rectalis). ",
      "Review the input table and candidate identifiers."
    )
  )
}


# -----------------------------------------------------------------------------
# 8. Build the base quadrant plot
# -----------------------------------------------------------------------------

p_base <- ggplot(
  plot_df,
  aes(
    x = pooled_rho,
    y = pooled_g
  )
) +

  # 378 background SGBs
  geom_point(
    data = plot_df %>% filter(Group == "Background"),
    color = "grey78",
    alpha = 0.65,
    size = 2.5
  ) +

  # 13 other shortlisted candidates
  geom_point(
    data = plot_df %>% filter(Group == "Other candidates"),
    color = "#2C7FB8",
    size = 3.7
  ) +

  # Agathobacter rectalis
  geom_point(
    data = plot_df %>% filter(Group == "Agathobacter rectalis"),
    color = "#D7301F",
    size = 4.4
  ) +

  # Zero reference lines
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey75",
    linewidth = 0.55
  ) +

  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey75",
    linewidth = 0.55
  ) +

  # Axis titles
  labs(
    x = expression(
      paste(
        "Pooled Spearman ",
        rho,
        " with ",
        italic("F. nucleatum")
      )
    ),
    y = "Pooled Hedges' g (CRC vs control)"
  ) +

  # X-axis scale
  scale_x_continuous(
    breaks = c(
      -0.3, -0.2, -0.1, 0.0, 0.1,
      0.2, 0.3, 0.4, 0.5
    ),
    labels = c(
      "-0.3", "-0.2", "-0.1", "0.0", "0.1",
      "0.2", "0.3", "0.4", "0.5"
    ),
    expand = expansion(mult = c(0, 0))
  ) +

  # Y-axis scale
  scale_y_continuous(
    breaks = c(-0.5, 0.0, 0.5, 1.0),
    labels = c("-0.5", "0.0", "0.5", "1.0"),
    expand = expansion(mult = c(0, 0))
  ) +

  # Visible plotting range
  coord_cartesian(
    xlim = c(-0.3, 0.5),
    ylim = c(-0.5, 1.1),
    clip = "off"
  ) +

  # Figure theme
  theme_classic(base_size = 16) +

  theme(
    axis.title.x = element_text(
      size = 18,
      color = "black",
      margin = margin(t = 12)
    ),
    axis.title.y = element_text(
      size = 18,
      color = "black",
      margin = margin(r = 12)
    ),
    axis.text.x = element_text(
      size = 14,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 14,
      color = "black"
    ),
    axis.line = element_line(
      color = "black",
      linewidth = 0.55
    ),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.45
    ),
    axis.ticks.length = grid::unit(0.13, "cm"),
    legend.position = "none",
    plot.margin = margin(
      t = 18,
      r = 25,
      b = 18,
      l = 18
    )
  )


# -----------------------------------------------------------------------------
# 9. Figure 1B without taxon labels
# -----------------------------------------------------------------------------

p_unlabeled <- p_base

print(p_unlabeled)

unlabeled_file <- file.path(
  output_dir,
  "Figure1B_14candidates_unlabeled.tiff"
)

ggsave(
  filename = unlabeled_file,
  plot = p_unlabeled,
  device = "tiff",
  width = 7.2,
  height = 6.2,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# -----------------------------------------------------------------------------
# 10. Figure 1B with all 14 taxon labels
# -----------------------------------------------------------------------------

p_labeled <- p_base +

  geom_text_repel(
    data = candidate14,
    aes(
      x = pooled_rho,
      y = pooled_g,
      label = Taxon_label
    ),

    # Reproducible label placement
    seed = 20260814,

    # Label appearance and spacing
    size = 4.2,
    box.padding = 0.55,
    point.padding = 0.35,

    # Connector segments
    min.segment.length = 0,
    segment.color = "grey45",
    segment.size = 0.40,

    # Keep all labels and repel them in both directions
    max.overlaps = Inf,
    force = 3.5,
    direction = "both",

    # Keep labels within the plotting range
    xlim = c(-0.30, 0.50),
    ylim = c(-0.50, 1.10)
  )

print(p_labeled)

labeled_file <- file.path(
  output_dir,
  "Figure1B_14candidates_labeled.tiff"
)

ggsave(
  filename = labeled_file,
  plot = p_labeled,
  device = "tiff",
  width = 9.0,
  height = 7.2,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# -----------------------------------------------------------------------------
# 11. Final checks
# -----------------------------------------------------------------------------

message("\n====================================================")
message("Figure 1B analysis script completed.")
message("====================================================")

message("Total SGBs plotted: ", nrow(plot_df))
message(
  "Grey background SGBs: ",
  sum(plot_df$Group == "Background")
)
message(
  "Blue candidate SGBs: ",
  sum(plot_df$Group == "Other candidates")
)
message(
  "Red Agathobacter rectalis: ",
  sum(plot_df$Group == "Agathobacter rectalis")
)
message(
  "Total highlighted candidates: ",
  sum(plot_df$Group != "Background")
)

message("\nTIFF files written:")
message("1. ", unlabeled_file)
message("2. ", labeled_file)

if (!file.exists(unlabeled_file) || !file.exists(labeled_file)) {
  stop("One or more TIFF output files were not created.", call. = FALSE)
}
