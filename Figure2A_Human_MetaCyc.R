# =============================================================================
# Figure 2A. Human MetaCyc pathway scatter plot (430 pathways)
#
# Purpose
#   Classify 430 MetaCyc pathways according to the direction and FDR-adjusted
#   significance of their associations with Agathobacter rectalis (Ar) and
#   Fusobacterium nucleatum (Fn), then generate two publication-style scatter
#   plots.
#
# Classification logic
#   1. Opposite direction: Ar_rho * Fn_rho < 0
#   2. Fn-related significant pathways:
#        Ar_rho < 0, Fn_rho > 0, Ar_q < 0.20, Fn_q < 0.20
#   3. Ar-related significant pathways:
#        Ar_rho > 0, Fn_rho < 0, Ar_q < 0.20, Fn_q < 0.20
#   4. Other opposite pathways:
#        - both FDR q values >= 0.20 (expected n = 155), or
#        - exactly one FDR q value < 0.20 (expected n = 131)
#   5. Non-opposite pathways:
#        all records not assigned to one of the opposite-direction groups
#        above (expected n = 91)
#
# Expected totals
#   - Total pathways: 430
#   - Non-opposite pathways: 91
#   - Other opposite pathways: 286 (155 + 131)
#   - Significant opposite pathways: 53 (Fn-related + Ar-related)
#   - Total opposite pathways: 339
#
# Input
#   data/input/Figure2_pathway_scatter_430_R_ready.xlsx
#   Worksheet: 01_Scatter430
#
# Outputs
#   results/figures/Figure2_PathwayScatter_FDR_direction_with_text.tiff
#   results/figures/Figure2_PathwayScatter_FDR_direction_no_text.tiff
#
# Important
#   This script starts from a processed pathway-level summary workbook. It
#   does not perform HUMAnN pathway quantification, cohort-level correlations,
#   meta-analysis, or FDR adjustment.
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Check and load required packages
# -----------------------------------------------------------------------------

required_packages <- c(
  "here",
  "readxl",
  "dplyr",
  "ggplot2"
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


# -----------------------------------------------------------------------------
# 1. Define project paths
# -----------------------------------------------------------------------------
# Open Fn_Ar_CRC_Analysis.Rproj in RStudio before running this script.

project_root <- here::here()

input_file <- here::here(
  "data",
  "input",
  "Figure2_pathway_scatter_430_R_ready.xlsx"
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
      "\n\nPlace 'Figure2_pathway_scatter_430_R_ready.xlsx' ",
      "in the data/input/ folder."
    ),
    call. = FALSE
  )
}


# -----------------------------------------------------------------------------
# 2. Read pathway-level summary data
# -----------------------------------------------------------------------------

input_sheet <- "01_Scatter430"

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
# 3. Check required columns and data types
# -----------------------------------------------------------------------------

required_columns <- c(
  "Pathway",
  "Pathway_ID",
  "Ar_rho",
  "Fn_rho",
  "Ar_q",
  "Fn_q"
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
  "Ar_rho",
  "Fn_rho",
  "Ar_q",
  "Fn_q"
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

if (anyNA(df[numeric_columns])) {
  warning(
    paste0(
      "Missing values were detected in one or more statistical columns. ",
      "The original classification logic assigns records that do not meet ",
      "an opposite-direction condition to 'Non-opposite pathways'. ",
      "Review missing values before interpreting the group counts."
    ),
    call. = FALSE
  )
}

if (nrow(df) != 430) {
  warning(
    paste0(
      "Expected 430 pathways in the input worksheet, but found ",
      nrow(df),
      ". Check the workbook before interpreting the figure."
    ),
    call. = FALSE
  )
}


# -----------------------------------------------------------------------------
# 4. Define opposite Ar-Fn associations
# -----------------------------------------------------------------------------
# A pathway is defined as opposite-direction when the product of the two
# pooled Spearman correlation coefficients is negative.

df <- df %>%
  mutate(
    opposite_direction = Ar_rho * Fn_rho < 0
  )


# -----------------------------------------------------------------------------
# 5. Create display groups
# -----------------------------------------------------------------------------
# The order of conditions is retained from the supplied source code.

df <- df %>%
  mutate(
    Display_Group = case_when(

      # Significant Fn-related pathways:
      # Ar negative / Fn positive; both FDR q values < 0.20
      opposite_direction &
        Ar_q < 0.20 &
        Fn_q < 0.20 &
        Ar_rho < 0 &
        Fn_rho > 0 ~
        "Fn-related pathways",

      # Significant Ar-related pathways:
      # Ar positive / Fn negative; both FDR q values < 0.20
      opposite_direction &
        Ar_q < 0.20 &
        Fn_q < 0.20 &
        Ar_rho > 0 &
        Fn_rho < 0 ~
        "Ar-related pathways",

      # Opposite-direction pathways with exactly one FDR q value < 0.20
      opposite_direction &
        (
          (Ar_q < 0.20 & Fn_q >= 0.20) |
            (Ar_q >= 0.20 & Fn_q < 0.20)
        ) ~
        "One FDR < 0.20",

      # Opposite-direction pathways with both FDR q values >= 0.20
      opposite_direction &
        Ar_q >= 0.20 &
        Fn_q >= 0.20 ~
        "Opposite, both FDR >= 0.20",

      # Remaining records
      TRUE ~
        "Non-opposite pathways"
    )
  )


# -----------------------------------------------------------------------------
# 6. Calculate and report pathway group counts
# -----------------------------------------------------------------------------

message("\n========================================")
message("PATHWAY GROUP COUNTS")
message("========================================")

print(table(df$Display_Group))

n_nonopposite <- sum(
  df$Display_Group == "Non-opposite pathways"
)

n_fdr_none <- sum(
  df$Display_Group == "Opposite, both FDR >= 0.20"
)

n_fdr_one <- sum(
  df$Display_Group == "One FDR < 0.20"
)

n_fn_related <- sum(
  df$Display_Group == "Fn-related pathways"
)

n_ar_related <- sum(
  df$Display_Group == "Ar-related pathways"
)

n_other_opposite <- n_fdr_none + n_fdr_one
n_significant <- n_fn_related + n_ar_related

n_opposite <- sum(
  df$opposite_direction,
  na.rm = TRUE
)

message("\nExpected structure:")
message("Non-opposite = 91")
message("Other opposite = 286")
message("Fn-related + Ar-related = 53")
message("Total opposite = 339")

message("\nActual counts:")
message("Total pathways = ", nrow(df))
message("Non-opposite = ", n_nonopposite)
message("Opposite, both FDR >= 0.20 = ", n_fdr_none)
message("Opposite, exactly one FDR < 0.20 = ", n_fdr_one)
message("Other opposite = ", n_other_opposite)
message("Fn-related = ", n_fn_related)
message("Ar-related = ", n_ar_related)
message("Fn-related + Ar-related = ", n_significant)
message("Total opposite = ", n_opposite)


# -----------------------------------------------------------------------------
# 7. Define colors
# -----------------------------------------------------------------------------

col_grey <- "#E6E6E6"
col_lightblue <- "#BDD5F1"
col_fn <- "#F4B183"
col_ar <- "#A9D18E"


# -----------------------------------------------------------------------------
# 8. Split data into plotting layers
# -----------------------------------------------------------------------------
# The 155- and 131-pathway groups remain separate to preserve their original
# point sizes, although both groups use the same light-blue fill color.

df_grey <- df %>%
  filter(Display_Group == "Non-opposite pathways")

df_lightblue_155 <- df %>%
  filter(Display_Group == "Opposite, both FDR >= 0.20")

df_lightblue_131 <- df %>%
  filter(Display_Group == "One FDR < 0.20")

df_fn <- df %>%
  filter(Display_Group == "Fn-related pathways")

df_ar <- df %>%
  filter(Display_Group == "Ar-related pathways")


# -----------------------------------------------------------------------------
# 9. Calculate symmetric plot range
# -----------------------------------------------------------------------------

lim <- max(
  abs(c(df$Ar_rho, df$Fn_rho)),
  na.rm = TRUE
) + 0.04

if (!is.finite(lim) || lim <= 0) {
  stop(
    "Could not calculate a valid symmetric axis range from Ar_rho and Fn_rho.",
    call. = FALSE
  )
}


# -----------------------------------------------------------------------------
# 10. Define axis breaks
# -----------------------------------------------------------------------------

x_breaks <- pretty(
  c(-lim, lim),
  n = 5
)

y_breaks <- pretty(
  c(-lim, lim),
  n = 5
)


# -----------------------------------------------------------------------------
# 11. Build common scatter-plot layers
# -----------------------------------------------------------------------------

p_layers <- ggplot() +

  # Zero reference lines
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.50,
    colour = "#777777"
  ) +

  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.50,
    colour = "#777777"
  ) +

  # Layer 1: 91 non-opposite pathways (light grey)
  geom_point(
    data = df_grey,
    aes(
      x = Ar_rho,
      y = Fn_rho,
      fill = "Other pathways"
    ),
    shape = 21,
    size = 3.3,
    stroke = 0.18,
    colour = "white",
    alpha = 0.65
  ) +

  # Layer 2: 155 opposite pathways with both FDR q values >= 0.20
  geom_point(
    data = df_lightblue_155,
    aes(
      x = Ar_rho,
      y = Fn_rho,
      fill = "Other opposite pathways"
    ),
    shape = 21,
    size = 3.8,
    stroke = 0.22,
    colour = "white",
    alpha = 0.82
  ) +

  # Layer 3: 131 opposite pathways with exactly one FDR q value < 0.20
  geom_point(
    data = df_lightblue_131,
    aes(
      x = Ar_rho,
      y = Fn_rho,
      fill = "Other opposite pathways"
    ),
    shape = 21,
    size = 4.5,
    stroke = 0.30,
    colour = "white",
    alpha = 0.95
  ) +

  # Layer 4: significant Fn-related pathways (Ar- / Fn+; light orange)
  geom_point(
    data = df_fn,
    aes(
      x = Ar_rho,
      y = Fn_rho,
      fill = "Fn-related pathways"
    ),
    shape = 21,
    size = 5.3,
    stroke = 0.45,
    colour = "white",
    alpha = 1
  ) +

  # Layer 5: significant Ar-related pathways (Ar+ / Fn-; light green)
  geom_point(
    data = df_ar,
    aes(
      x = Ar_rho,
      y = Fn_rho,
      fill = "Ar-related pathways"
    ),
    shape = 21,
    size = 5.3,
    stroke = 0.45,
    colour = "white",
    alpha = 1
  ) +

  # Manual fill colors and legend labels
  scale_fill_manual(
    values = c(
      "Fn-related pathways" = col_fn,
      "Ar-related pathways" = col_ar,
      "Other opposite pathways" = col_lightblue,
      "Other pathways" = col_grey
    ),
    breaks = c(
      "Fn-related pathways",
      "Ar-related pathways",
      "Other opposite pathways",
      "Other pathways"
    ),
    labels = c(
      "Fn-related pathways",
      "Ar-related pathways",
      "Other opposite Ar–Fn pathways",
      "Other pathways"
    ),
    name = NULL
  ) +

  # Symmetric axes
  scale_x_continuous(
    breaks = x_breaks,
    expand = expansion(mult = c(0.02, 0.02))
  ) +

  scale_y_continuous(
    breaks = y_breaks,
    expand = expansion(mult = c(0.02, 0.02))
  ) +

  coord_cartesian(
    xlim = c(-lim, lim),
    ylim = c(-lim, lim),
    clip = "off"
  )


# -----------------------------------------------------------------------------
# 12. Version 1: annotated plot with legend
# -----------------------------------------------------------------------------

p_text <- p_layers +

  # Upper-left quadrant annotation
  annotate(
    "text",
    x = -lim * 0.93,
    y = lim * 0.90,
    label = "Ar− / Fn+",
    hjust = 0,
    vjust = 1,
    size = 4.2,
    fontface = "italic",
    colour = "#333333"
  ) +

  # Lower-right quadrant annotation
  annotate(
    "text",
    x = lim * 0.93,
    y = -lim * 0.90,
    label = "Ar+ / Fn−",
    hjust = 1,
    vjust = 0,
    size = 4.2,
    fontface = "italic",
    colour = "#333333"
  ) +

  # Main message
  annotate(
    "text",
    x = 0,
    y = lim * 0.97,
    label = "339/430 pathways showed opposite Ar–Fn associations",
    hjust = 0.5,
    vjust = 1,
    size = 4.0,
    colour = "#333333"
  ) +

  # Axis labels
  labs(
    x = expression(
      paste(
        "Pooled Spearman ",
        rho,
        " with ",
        italic("A. rectalis")
      )
    ),
    y = expression(
      paste(
        "Pooled Spearman ",
        rho,
        " with ",
        italic("F. nucleatum")
      )
    )
  ) +

  theme_classic(base_size = 12) +

  theme(
    axis.title = element_text(
      size = 14,
      colour = "black"
    ),
    axis.text = element_text(
      size = 11.5,
      colour = "black"
    ),
    axis.line = element_line(
      linewidth = 0.65,
      colour = "black"
    ),
    axis.ticks = element_line(
      linewidth = 0.50,
      colour = "black"
    ),
    legend.position = c(0.25, 0.15),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 9.0),
    legend.spacing.y = grid::unit(1.2, "mm"),
    plot.margin = margin(10, 15, 10, 10)
  ) +

  guides(
    fill = guide_legend(
      override.aes = list(
        size = 5.0,
        alpha = 1,
        shape = 21,
        colour = "white",
        stroke = 0.35
      )
    )
  )

print(p_text)


# -----------------------------------------------------------------------------
# 13. Save annotated version
# -----------------------------------------------------------------------------

with_text_file <- file.path(
  output_dir,
  "Figure2_PathwayScatter_FDR_direction_with_text.tiff"
)

ggsave(
  filename = with_text_file,
  plot = p_text,
  device = "tiff",
  width = 9.2,
  height = 7.6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# -----------------------------------------------------------------------------
# 14. Version 2: minimal plot
# -----------------------------------------------------------------------------
# Retains points, colors, axis lines, axis ticks, axis numbers, and zero lines.
# Removes axis titles, legend, quadrant annotations, and the summary sentence.

p_minimal <- p_layers +

  labs(
    x = NULL,
    y = NULL
  ) +

  theme_classic(base_size = 12) +

  theme(
    axis.title = element_blank(),

    axis.text.x = element_text(
      size = 12,
      colour = "black",
      margin = margin(t = 4)
    ),

    axis.text.y = element_text(
      size = 12,
      colour = "black",
      margin = margin(r = 4)
    ),

    axis.ticks = element_line(
      linewidth = 0.55,
      colour = "black"
    ),

    axis.line = element_line(
      linewidth = 0.70,
      colour = "black"
    ),

    legend.position = "none",

    plot.margin = margin(6, 6, 6, 6)
  )

print(p_minimal)


# -----------------------------------------------------------------------------
# 15. Save minimal version
# -----------------------------------------------------------------------------

no_text_file <- file.path(
  output_dir,
  "Figure2_PathwayScatter_FDR_direction_no_text.tiff"
)

ggsave(
  filename = no_text_file,
  plot = p_minimal,
  device = "tiff",
  width = 7.2,
  height = 7.2,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# -----------------------------------------------------------------------------
# 16. Final summary and expected-count checks
# -----------------------------------------------------------------------------

message("\n========================================")
message("FINAL SUMMARY")
message("========================================")
message("Total pathways = ", nrow(df))
message("Non-opposite pathways = ", n_nonopposite)
message("Other opposite pathways = ", n_other_opposite)
message("Fn-related pathways = ", n_fn_related)
message("Ar-related pathways = ", n_ar_related)
message("Significant pathways total = ", n_significant)
message("Total opposite = ", n_opposite)

if (n_nonopposite != 91) {
  warning(
    paste0(
      "Non-opposite pathway count is ",
      n_nonopposite,
      ", not 91."
    ),
    call. = FALSE
  )
}

if (n_fdr_none != 155) {
  warning(
    paste0(
      "Opposite pathways with both FDR q values >= 0.20: ",
      n_fdr_none,
      ", not 155."
    ),
    call. = FALSE
  )
}

if (n_fdr_one != 131) {
  warning(
    paste0(
      "Opposite pathways with exactly one FDR q value < 0.20: ",
      n_fdr_one,
      ", not 131."
    ),
    call. = FALSE
  )
}

if (n_other_opposite != 286) {
  warning(
    paste0(
      "Other opposite pathway count is ",
      n_other_opposite,
      ", not 286."
    ),
    call. = FALSE
  )
}

if (n_significant != 53) {
  warning(
    paste0(
      "Fn-related + Ar-related pathway count is ",
      n_significant,
      ", not 53."
    ),
    call. = FALSE
  )
}

if (n_opposite != 339) {
  warning(
    paste0(
      "Total opposite pathway count is ",
      n_opposite,
      ", not 339."
    ),
    call. = FALSE
  )
}


# -----------------------------------------------------------------------------
# 17. Verify output files
# -----------------------------------------------------------------------------

outputs <- c(
  with_text_file,
  no_text_file
)

output_status <- data.frame(
  File = basename(outputs),
  Exists = file.exists(outputs)
)

print(output_status)

if (!all(output_status$Exists)) {
  stop(
    "One or more expected TIFF output files were not created.",
    call. = FALSE
  )
}

message("\nFigure 2A pathway scatter analysis completed.")
