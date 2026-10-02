###################################################
###################################################

 ## Expanding ASV selection to cluster amplicons ##

###################################################
###################################################

# Background:
#   Every ASV in EukBank is assigned to a hierarchical single-linkage
#   cluster, built from pairwise distances, making it possible to find, for a given ASV,
#   the largest cluster it belongs to.

# Retrieves all amplicons from a list of selected EukBank v1 ASVs,
# expanding each selected ASV to its full set of member amplicons.

# Note: not intended to be run from the terminal. Edit GROUP_NAME and the
# input filename below, then run interactively in R or with source().

# Input:
#   eukbank_18S_V4_asvs.tsv.gz           EukBank ASV table
#   eukbank_18S_V4_clusters_comp.tsv.gz  EukBank cluster composition table
#   EukBank-GROUP_NAME_selection.list
#                                        List of amplicons selected by BLAST + tree

# Output:
#   cluster_EukBank-GROUP_NAME.list
#   List of amplicons: the original selection, plus all amplicons
#   belonging to the same clusters.

# Example:
#   Set GROUP_NAME to "Choanoflagellatea", place the required EukBank
#   files and EukBank-Choanoflagellatea_selection.list
#   in the working directory, then run the script in R.

library(data.table)

EukBank_ASVs <- fread('eukbank_18S_V4_asvs.tsv.gz', quote="")
EukBank_clusters_comp <- fread('eukbank_18S_V4_clusters_comp.tsv.gz', quote="")

# Read the selected amplicons
blast_amplicons <- readLines('EukBank-GROUP_NAME_selection.list')

# Create vectors to store the amplicons
lonely_amplicons <- c()    # Amplicons that do not belong to a cluster
unknown_cluster_id <- c()  # Cluster IDs of amplicons that do belong to one
cluster_amplicons <- c()   # All amplicons belonging to those clusters
cluster_amplicons_temp <- c()

# For each selected amplicon, store its cluster ID if it belongs to one;
# otherwise, store it as a lonely amplicon
for (amplicon in blast_amplicons) {
  indices <- which(EukBank_ASVs[, "amplicon"] == amplicon)
  if (length(indices) > 0 && any(!is.na(EukBank_ASVs[indices, "unknown_cluster_id"]) && EukBank_ASVs[indices, "unknown_cluster_id"] != "")) {
    unknown_cluster_id <- c(unknown_cluster_id, as.character(EukBank_ASVs[indices, "unknown_cluster_id"]))
  } else {
    lonely_amplicons <- c(lonely_amplicons, amplicon)
  }
}

# Keep each cluster ID only once
unknown_cluster_id <- unique(sort(unknown_cluster_id))

# Retrieve all amplicons belonging to the selected clusters
for (cluster in unknown_cluster_id) {
  indices <- which(EukBank_clusters_comp$cluster == cluster)
  if (length(indices) > 0) {
    cluster_amplicons_temp <- as.vector(unlist(EukBank_clusters_comp[indices, "amplicon"]))
    cluster_amplicons <- c(cluster_amplicons, cluster_amplicons_temp)
  }
}

# Remove duplicates, in case an amplicon belongs to more than one cluster
cluster_amplicons <- unique(sort(cluster_amplicons))

# Combine lonely amplicons and cluster amplicons into the final list
final_amplicons <- c(lonely_amplicons, cluster_amplicons)

# Check that the final list matches the original selection
if (identical(sort(final_amplicons), sort(blast_amplicons))) {
  print("Same amplicons!")
} else {
  print("Mismatch: amplicon lists differ. :(")
}

write(final_amplicons, "cluster_EukBank-GROUP_NAME.list")