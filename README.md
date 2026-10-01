# Exploring Unicellular Holozoa Diversity


envdistribution.sh

Bash script that summarizes the distribution of a group of ASVs from EukBank v1 across samples and environments. For each input list, it sums the reads of all ASVs in the group per sample and retrieves associated data (total sample reads, environment, and geographic coordinates) to calculate the group's relative abundance in each sample. It also generates a table of reads per ASV and environment. The resulting tables can be used downstream, e.g. to generate biogeographic maps in R, or as metadata on a phylogenetic tree to explore environmental patterns within the group.

eukbank_ASVs_envbiogeography.R

Rscript that generates world maps showing the geographic distribution of an ASV group from EukBank v1. Using the per-sample summary table produced by envdistribution.sh, it plots one map sized by relative abundance and another sized by read count, with points colored by environment.

relative_abundance_UniHolozoa.sh

Bash script that computes read abundance per environment for a set of ASV groups from EukBank v1, alongside each known unicellular Holozoa lineage. For each ASV, it resolves whether it belongs to a user-defined group or a Holozoa taxon (giving group assignments priority when both apply) and sums reads by environment accordingly. The resulting table can be used downstream to compare group relative abundance across environments, e.g. to generate plots in R.

eukbank_relativeabundance_plot.R

Rscript that plots the read abundance of groups across environments from EukBank v1. Using a reads-per-environment table (e.g. from relative_abundance_UniHolozoa.sh), it generates a horizontal bar chart showing the proportional abundance of each group per environment, plus paginated faceted bar charts showing absolute read counts per group.


