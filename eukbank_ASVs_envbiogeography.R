#!/usr/bin/env Rscript

# Based on: http://sarahleejane.github.io/learning/r/2014/09/20/plotting-beautiful-clear-maps-with-r.html

library(maps)
library(ggplot2)
library(dplyr)
library(patchwork)
library(gridExtra)
library(scales)

###################################################
###################################################

          # From ASVs to environmental
        # distribution and biogeography
               # visualization

###################################################
###################################################

## Creates a map showing the geographic location, relative abundance/read
## count, and environment of an ASV group across EukBank samples.

## Usage:
#   Rscript eukbank_ASVs_envbiogeography.R <input_file.tsv>

## Input:
#   input_file.tsv   .tsv file with sample, reads, nreads, environment,
#                    latitude, longitude, and relative_abundance columns
#
#   The input filename must start with "ASVs_GroupX" (e.g. ASVs_Group1_...),
#   as this is used to name the output files and label the group in the map.
#
#   Recommended: the <basename>_reads_per_sample_nreads_env_lat_long_relativeabundance.tsv
#   output file from envdistribution.sh, found at:
#   https://github.com/alijimv/Exploring-Unicellular-Holozoa-Diversity/

## Output:
#   <group_name>_map_relativeabundance.pdf
#   World map with one point per sample, sized by relative abundance (%)
#   and colored by environment.
#
#   <group_name>_map_reads.pdf
#   World map with one point per sample, sized by read count and
#   colored by environment.
#
#   Files are written to the same directory as input_file, named after
#   the ASVs_GroupX pattern found in the input filename.

## Example:
#   Rscript plot_relative_abundance_map.R ASVs_Group1_reads_per_sample_nreads_env_lat_long_relativeabundance.tsv


# -----------------------------
# Arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript eukbank_ASVs_envbiogeography.R <input_file.tsv>")
}
input_file <- args[1]

# -----------------------------
# Paths and group name
# -----------------------------
input_dir  <- dirname(input_file)
group_name <- sub("(ASVs_Group[^_]+).*", "\\1", basename(input_file))

output_relative_file <- file.path(input_dir, paste0(group_name, "_map_relativeabundance.pdf"))
output_reads_file    <- file.path(input_dir, paste0(group_name, "_map_reads.pdf"))

# -----------------------------
# Load data
# -----------------------------
df <- read.delim(input_file, header = FALSE, stringsAsFactors = FALSE)

# Assign column names
colnames(df) <- c("Sample","Reads","SampleReads","Environment","Latitude","Longitude","Relative_Abundance")

# Convert numeric columns
df$Relative_Abundance <- as.numeric(gsub(",", ".", df$Relative_Abundance))
df$Latitude           <- as.numeric(df$Latitude)
df$Longitude          <- as.numeric(df$Longitude)
df$Reads              <- as.numeric(df$Reads)
df$SampleReads        <- as.numeric(df$SampleReads)

# -----------------------------
# Filter rows with coordinates
# -----------------------------
df_clean <- df %>% filter(!is.na(Latitude), !is.na(Longitude))

# -----------------------------
# Clean environment names
# -----------------------------
df_clean$Environment <- gsub("land_freshwater",  "Land freshwater",  df_clean$Environment)
df_clean$Environment <- gsub("land_organism",    "Land organism",    df_clean$Environment)
df_clean$Environment <- gsub("land_sediment",    "Land sediment",    df_clean$Environment)
df_clean$Environment <- gsub("land_soil",        "Land soil",        df_clean$Environment)
df_clean$Environment <- gsub("land_water",       "Land water",       df_clean$Environment)
df_clean$Environment <- gsub("marine_ice",       "Marine ice",       df_clean$Environment)
df_clean$Environment <- gsub("marine_organism",  "Marine organism",  df_clean$Environment)
df_clean$Environment <- gsub("marine_sediment",  "Marine sediment",  df_clean$Environment)
df_clean$Environment <- gsub("marine_water",     "Marine water",     df_clean$Environment)

# -----------------------------
# Colors for environments
# -----------------------------
colors <- c(
  "Marine water"     = "#c7eae5",
  "Marine ice"       = "#80cdc1",
  "Marine sediment"  = "#35978f",
  "Marine organism"  = "#01665e",
  "Land water"       = "#003c30",
  "Land freshwater"  = "#dfc27d",
  "Land sediment"    = "#bf812d",
  "Land soil"        = "#8c510a",
  "Land organism"    = "#543005"
)

# -----------------------------
# Base world map
# -----------------------------
world_map <- map_data("world")

base_world <- ggplot() +
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    colour = "light gray",
    fill   = "light gray"
  ) +
  coord_fixed() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom"
  ) +
  xlab("") + ylab("")

# -----------------------------
# Map: Relative Abundance
# -----------------------------
df_rel <- df_clean %>% filter(Relative_Abundance * 100 >= 0.01)
df_rel$Relative_Abundance_perc <- df_rel$Relative_Abundance * 100

# Sort so largest points are first (drawn behind)
df_rel <- df_rel %>% arrange(desc(Relative_Abundance_perc))

map_rel <- base_world +
  geom_point(
    data = df_rel,
    aes(
      x = Longitude,
      y = Latitude,
      size = Relative_Abundance_perc,
      fill = Environment
    ),
    pch = 21,
    alpha = 0.7,
    stroke = 0.1
  ) +
  scale_fill_manual(values = colors) +
  scale_size_continuous(
    name   = "Relative Abundance (%)",
    range  = c(2, 16),
    breaks = c(0.01, 0.1, 0.5, 1, 2, 3, 5, 8, 10),
    limits = c(0, 10)
  ) +
  geom_label(
    aes(x = 0, y = 120, label = group_name),
    hjust = 0.5,
    vjust = 1,
    fill  = "white",
    linewidth = 0.4
  )

if (nrow(df_rel) > 0) {
  pdf(output_relative_file, width = 15, height = 8)
  print(map_rel)
  dev.off()
} else {
  message("No data to plot relative abundance map for ", group_name)
}

# -----------------------------
# Map: Reads per sample
# -----------------------------
df_reads <- df_clean %>% filter(Reads >= 10)

# Sort so largest points are first (drawn behind)
df_reads <- df_reads %>% arrange(desc(Reads))

map_reads <- base_world +
  geom_point(
    data = df_reads,
    aes(
      x    = Longitude,
      y    = Latitude,
      size = Reads,
      fill = Environment
    ),
    pch    = 21,
    alpha  = 0.7,
    stroke = 0.1
  ) +
  scale_fill_manual(values = colors) +
  scale_size_continuous(
    name   = "Reads",
    range  = c(2, 16),
    breaks = c(10, 50, 100, 200, 1000, 10000, 15000, 20000),
    limits = c(0, 20000)
  ) +
  geom_label(
    aes(x = 0, y = 120, label = group_name),
    hjust = 0.5,
    vjust = 1,
    fill  = "white",
    linewidth = 0.4
  )

if (nrow(df_reads) > 0) {
  pdf(output_reads_file, width = 15, height = 8)
  print(map_reads)
  dev.off()
} else {
  message("No data to plot reads map for ", group_name)
}

message("Maps saved to folder: ", input_dir)



