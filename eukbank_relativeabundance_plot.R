#!/usr/bin/env Rscript

###################################################
###################################################

# From ASVs to environmental
# distribution and biogeography
# visualization

###################################################
###################################################

# Plots the relative abundance of a group across environments.
#
# Generates two sets of plots from a long-format reads-per-environment
# table (as produced by the EukBank biogeography pipeline):


# Input:
#   input_file.tsv   .tsv file with group, environment, reads (no header)
#                    Recommended: the UniHolozoa_AllGroups_reads_per_env.tsv
#                    output file from relative_abundance_UniHolozoa.sh found in
#                    https://github.com/alijimv/Exploring-Unicellular-Holozoa-Diversity/
#
#   colors.tsv       .tsv file with Group, Color
#                     Defines the color and plotting order of groups.
#                     Groups in input_file.tsv with no matching entry
#                     here are excluded, with a warning.
#
# Usage:
#   Rscript eukbank_relativeabundance_plot.R <input_file.tsv> <colors.tsv>
#
# Output:
#   <basename>_Proportion_Horizontal.pdf      
#   Plot A: a single horizontal bar chart showing the proportional
#   abundance of each group within each environment.

#   <basename>_AllGroups_PageN.pdf          
#   Plot B: faceted bar charts (one panel per group, paginated)
#   showing absolute read counts per environment.

# Example:
#   Rscript eukbank_relativeabundance_plot.R UniHolozoa_AllGroups_reads_per_env.tsv colors.tsv


rm(list = ls())
options(scipen = 999)

# -------------------------------------------------
# Read command-line arguments
# -------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript eukbank_relativeabundance_plot.R input_file.tsv colors.tsv")
}

input_file <- args[1]
color_file <- args[2]

if (!file.exists(input_file)) stop(paste("File not found:", input_file))
if (!file.exists(color_file)) stop(paste("Color file not found:", color_file))

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggforce)
})

# -------------------------------------------------
# Load input data
# -------------------------------------------------
data <- read.csv(input_file, header = FALSE, sep = "\t")
names(data) <- c("Group", "Environment", "Reads")

data$Environment <- gsub("_", " ", data$Environment)

# -------------------------------------------------
# Load color mapping
# -------------------------------------------------
mapping <- read.csv(color_file,
                    sep = "\t",
                    header = TRUE,
                    stringsAsFactors = FALSE)

if (!all(c("Group", "Color") %in% colnames(mapping))) {
  stop("Color file must contain columns: Group and Color")
}

# -------------------------------------------------
# Identify excluded Groups
# -------------------------------------------------
excluded <- setdiff(unique(data$Group), mapping$Group)

if (length(excluded) > 0) {
  warning(paste(
    "Excluded groups (no color defined):",
    paste(excluded, collapse = ", ")
  ))
}

# -------------------------------------------------
# Keep only groups with defined colors
# -------------------------------------------------
data <- data[data$Group %in% mapping$Group, ]

# -------------------------------------------------
# Enforce group ordering
# -------------------------------------------------
mapping$Group <- factor(mapping$Group, levels = mapping$Group)
data$Group    <- factor(data$Group, levels = mapping$Group)

# Merge colors
data <- merge(data, mapping, by = "Group", all.x = TRUE)

color_vector <- setNames(mapping$Color, mapping$Group)

# -------------------------------------------------
# Plot A: Horizontal proportional abundance
# (fixed environment order only in plot)
# -------------------------------------------------
plotA <- ggplot(
  data,
  aes(
    y = factor(
      Environment,
      levels = c(
        "marine water",
        "marine ice",
        "marine sediment",
        "marine organism",
        "land water",
        "land freshwater",
        "land sediment",
        "land soil",
        "land organism"
      )
    ),
    x = Reads,
    fill = Group
  )
) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_vector, name = "Group") +
  scale_x_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9)
  ) +
  labs(
    x = "Proportion of Reads",
    y = NULL
  )

pdf(
  sub("\\.tsv$", "_Proportion_Horizontal.pdf", basename(input_file)),
  width = 12,
  height = 8
)

print(plotA)
dev.off()

cat("Plot A saved.\n")

# -------------------------------------------------
# Plot B: Absolute counts with faceted groups
# -------------------------------------------------
plotB <- ggplot(
  data,
  aes(
    x = Environment,
    y = Reads,
    fill = Group
  )
) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = color_vector, guide = "none") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 7
    ),
    axis.line = element_line()
  ) +
  labs(y = "Number of Reads")

# -------------------------------------------------
# Facet pagination
# -------------------------------------------------
n_groups <- length(levels(data$Group))
ncol_per_page <- 4
nrow_per_page <- 5
panels_per_page <- ncol_per_page * nrow_per_page
total_pages <- ceiling(n_groups / panels_per_page)

for (page in 1:total_pages) {
  
  output_file <- sub(
    "\\.tsv$",
    paste0("_AllGroups_Page", page, ".pdf"),
    basename(input_file)
  )
  
  pdf(output_file,
      width = 11,
      height = 14)
  
  print(
    plotB +
      ggforce::facet_wrap_paginate(
        ~Group,
        scales = "free_y",
        ncol = ncol_per_page,
        nrow = nrow_per_page,
        page = page
      )
  )
  
  dev.off()
  
  cat("Plot B page", page, "saved.\n")
}