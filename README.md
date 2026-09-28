# Exploring Unicellular Holozoa Diversity


envdistribution.sh

Bash script that extracts sample-level metadata for a given list of ASVs from EukBank v1. For each ASV, it identifies the samples in which it occurs and retrieves associated data — read counts, environment, and geographic coordinates, and calculates relative abundance within the sample. The resulting tables can be used downstream, e.g. to generate biogeographic maps in R.



relative_abundance_UniHolozoa.sh

Bash script that computes read abundance per environment for a set of ASV groups from EukBank v1, alongside each known unicellular Holozoa lineage. For each ASV, it resolves whether it belongs to a user-defined group or a Holozoa taxon (giving group assignments priority when both apply) and sums reads by environment accordingly. The resulting table can be used downstream to compare group relative abundance across environments, e.g. to generate plots in R.

